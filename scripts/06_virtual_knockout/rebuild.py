from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd() / "R"))
from runtime import params
import os, sys, json, hashlib, time, gc, inspect, platform
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import celloracle as co
from scipy import sparse

ROOT = Path(".")
OUT = ROOT / "results/virtual_knockout"
for d in ["objects", "tables", "sources", "protocol", "provenance", "qa", "logs"]:
    (OUT / d).mkdir(parents=True, exist_ok=True)
PRIOR = ROOT / "data/reference/promoter/hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10.parquet"
SPEC = {
    "lung": ("lung_umi", "Alveolar_macrophages", "Alveolar_macrophages_UMIs.h5ad", 30000, 56),
    "brain": ("brain_umi", "Micro-PVM_2", "three_states_UMIs.h5ad", 25056, 80),
}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def writej(p, x):
    p = Path(p)
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(json.dumps(x, indent=2, default=str))
    tmp.replace(p)


def status(side, stage, **kw):
    d = dict(side=side, stage=stage, time=time.strftime("%Y-%m-%d %H:%M:%S"), **kw)
    writej(OUT / f"status_{side}.json", d)
    print(d, flush=True)


def matsha(x):
    return hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()


def run(side):
    np.random.seed(params["celloracle_seed"])
    assert co.__version__ == "0.20.0"
    folder, state, umi, n, nd = SPEC[side]
    base = ROOT / "results/virtual_knockout/input" / folder
    checkpoint = OUT / f"objects/{side}_imputed.celloracle.oracle"
    rawlink = OUT / f"objects/{side}_raw.celloracle.links"
    ready = OUT / f"objects/{side}_ready.celloracle.oracle"
    paths = [
        base / "objects" / umi,
        PRIOR,
        ROOT / "data/reference/promoter/humanTF_v1.01.csv",
        ROOT / "config/params.yml",
        Path(__file__),
        ROOT / "results/module_overlap/results/cross/AxisA_confirmed_genes.csv",
        ROOT / "results/single_cell/provenance/Entrez_Ensembl_annotation.csv",
    ]
    manifest = pd.DataFrame(
        [dict(path=str(p.resolve()), size=p.stat().st_size, sha256=sha(p)) for p in paths]
    )
    mp = OUT / f"qa/{side}_input_manifest.tsv"
    if mp.exists():
        pd.testing.assert_frame_equal(pd.read_csv(mp, sep="\t"), manifest)
    else:
        manifest.to_csv(mp, sep="\t", index=False)
    axis_genes = pd.read_csv(paths[-2], dtype=str)
    assert len(axis_genes) == 137 and axis_genes.entrez_id.is_unique
    if not checkpoint.exists():
        status(side, "load_and_select")
        b = ad.read_h5ad(paths[0])
        a = b[b.obs.state.eq(state)].copy()
        del b
        gc.collect()
        assert a.obs_names.is_unique and a.var_names.is_unique and a.obs.donor.nunique() == nd
        assert sparse.issparse(a.X)
        assert (
            np.isfinite(a.X.data).all()
            and (a.X.data >= 0).all()
            and np.equal(a.X.data, np.floor(a.X.data)).all()
        )
        prior = pd.read_parquet(PRIOR)
        regs = set(prior.columns) - {"peak_id", "gene_short_name"}
        tf = pd.read_csv(paths[2])
        assert tf["HGNC symbol"].is_unique
        bysym = tf.set_index("HGNC symbol")
        byen = tf.set_index("Ensembl ID")
        encol = next(
            (c for c in ["gene_ids", "gene_id", "gene_identifier", "ensembl_id"] if c in a.var), None
        )
        assert encol is not None
        feature = pd.DataFrame(
            {
                "gene": a.var_names.astype(str),
                "ENSEMBL": a.var[encol].astype(str).str.split(".").str[0].values,
            }
        )
        annotation = pd.read_csv(paths[-1], dtype=str).dropna().drop_duplicates()
        pairs = annotation[annotation.ENTREZID.isin(axis_genes.entrez_id)].merge(
            feature, on="ENSEMBL", how="inner"
        )
        pairs["n_features"] = pairs.groupby("ENTREZID").gene.transform("nunique")
        pairs["n_IDs_for_feature"] = pairs.groupby("gene").ENTREZID.transform("nunique")
        unique = pairs[(pairs.n_features == 1) & (pairs.n_IDs_for_feature == 1)]
        coverage = (
            axis_genes.merge(
                unique[["ENTREZID", "ENSEMBL", "gene"]],
                left_on="entrez_id",
                right_on="ENTREZID",
                how="left",
                validate="one_to_one",
            )
            .drop(columns="ENTREZID")
            .rename(columns={"ENSEMBL": "ensembl_id"})
        )
        coverage["mapping_rule"] = (
            "unique measurable Entrez to Ensembl to original feature; no symbol/alias rescue"
        )
        coverage.to_csv(OUT / f"tables/{side}_AxisA_stable_ID_ledger.tsv", sep="\t", index=False)
        pairs.to_csv(OUT / f"provenance/{side}_AxisA_all_stable_candidates.tsv", sep="\t", index=False)
        axis = set(coverage.gene.dropna())
        assert len(axis) == 136
        rows = []
        detected = np.asarray((a.X > 0).sum(axis=0)).ravel()
        tot = np.asarray(a.X.sum(axis=0)).ravel()
        donor_detect = np.zeros(a.n_vars, dtype=int)
        for donor in a.obs.donor.unique():
            donor_detect += np.asarray(a.X[a.obs.donor.eq(donor).values].sum(axis=0)).ravel() > 0
        for i, g in enumerate(a.var_names):
            en = str(a.var.iloc[i][encol]).split(".")[0] if encol else ""
            if g in bysym.index:
                r = bysym.loc[g]
                basis = "exact_symbol"
            elif en and en in byen.index and isinstance(byen.loc[en], pd.Series):
                r = byen.loc[en]
                basis = "unique_Ensembl"
            else:
                r = pd.Series(dtype=object)
                basis = "unmapped"
            rows.append(
                dict(
                    gene=g,
                    annotation_symbol=g if basis == "exact_symbol" else r.get("HGNC symbol", ""),
                    annotation_match=basis,
                    is_TF=r.get("Is TF?", "unannotated"),
                    TF_assessment=r.get("TF assessment", ""),
                    binding_mode=r.get("Binding mode", ""),
                    prior_regulator=g in regs,
                    in_UMI=True,
                    total_UMI=float(tot[i]),
                    cells_detected=int(detected[i]),
                    donors_detected=int(donor_detect[i]),
                )
            )
        ledger = pd.DataFrame(rows).set_index("gene", drop=False)
        for g in sorted((set(bysym[bysym["Is TF?"].eq("Yes")].index) | regs) - set(ledger.index)):
            r = bysym.loc[g] if g in bysym.index else pd.Series(dtype=object)
            ledger.loc[g] = dict(
                gene=g,
                annotation_symbol=g if len(r) else "",
                annotation_match="exact_symbol" if len(r) else "unmapped",
                is_TF=r.get("Is TF?", "unannotated"),
                TF_assessment=r.get("TF assessment", ""),
                binding_mode=r.get("Binding mode", ""),
                prior_regulator=g in regs,
                in_UMI=False,
                total_UMI=0,
                cells_detected=0,
                donors_detected=0,
            )
        all_genes = set(a.var_names)
        sc.pp.filter_genes(a, min_counts=1)
        nonzero = set(a.var_names)
        sc.pp.normalize_per_cell(a, key_n_counts="n_counts_all")
        fr = sc.pp.filter_genes_dispersion(a.X, flavor="cell_ranger", n_top_genes=params["hvg_n"], log=False)
        tfkeep = ledger.loc[a.var_names, "is_TF"].eq("Yes").values
        selected = fr.gene_subset | a.var_names.isin(axis) | tfkeep
        hvg = set(a.var_names[fr.gene_subset])
        selected_genes = set(a.var_names[selected])
        ledger["nonzero_before_HVG"] = ledger.gene.isin(nonzero)
        ledger["HVG"] = ledger.gene.isin(hvg)
        ledger["selected"] = ledger.gene.isin(selected_genes)
        a = a[:, selected].copy()
        sc.pp.normalize_per_cell(a)
        if side == "lung":
            sc.pp.subsample(a, n_obs=params["lung_model_cells"], random_state=params["celloracle_seed"])
        assert a.n_obs == n and a.obs.donor.nunique() == nd
        a.obs.to_csv(OUT / f"tables/{side}_cells.tsv.gz", sep="\t")
        ledger.to_csv(OUT / f"tables/{side}_feature_ledger.tsv.gz", sep="\t", index=False)
        cov = coverage.copy()
        cov["in_UMI"] = cov.gene.isin(all_genes)
        cov["nonzero"] = cov.gene.isin(nonzero)
        cov["in_model"] = cov.gene.isin(a.var_names)
        cov["missing_reason"] = np.select(
            [cov.in_model, ~cov.in_UMI, ~cov.nonzero],
            ["included", "absent_or_ambiguous_stable_feature", "no_detected_UMI_in_state"],
            default="not_selected",
        )
        assert not cov.loc[cov.nonzero, "in_model"].eq(False).any()
        cov.to_csv(OUT / f"tables/{side}_AxisA_coverage.tsv", sep="\t", index=False)
        if side == "lung":
            e = a.copy()
            sc.pp.log1p(e)
            sc.tl.pca(e, svd_solver="arpack", random_state=params["celloracle_seed"])
            a.obsm["X_pca_2d"] = e.obsm["X_pca"][:, :2].copy()
            del e
        oracle = co.Oracle()
        oracle.import_anndata_as_raw_count(
            a, cluster_column_name="state", embedding_name="X_pca_2d" if side == "lung" else "X_umap"
        )
        oracle.import_TF_data(TF_info_matrix=prior)
        del a, prior
        gc.collect()
        status(side, "PCA", genes=oracle.adata.n_vars)
        oracle.perform_PCA()
        vr = oracle.pca.explained_variance_ratio_
        npc = min(int(np.where(np.diff(np.diff(np.cumsum(vr)) > 0.002))[0][0]), 50)
        assert npc > 0
        pd.DataFrame({"PC": np.arange(1, len(vr) + 1), "variance_ratio": vr}).to_csv(
            OUT / f"tables/{side}_PCA.tsv", sep="\t", index=False
        )
        k = int(0.025 * n)
        status(side, "knn_imputation", n_pcs=npc, k=k)
        oracle.knn_imputation(
            n_pca_dims=npc,
            k=k,
            balanced=True,
            b_sight=k * 8,
            b_maxl=k * 4,
            n_jobs=params["celloracle_model_jobs"],
        )
        assert np.isfinite(oracle.adata.layers["imputed_count"]).all()
        oracle.to_hdf5(str(checkpoint))
        writej(
            OUT / f"qa/{side}_preprocessing.json",
            dict(
                cells=n,
                donors=nd,
                genes=oracle.adata.n_vars,
                n_pcs=npc,
                k=k,
                original_UMI_valid=True,
                cell_selection_checked=True,
                HVG_method="cell_ranger",
                Ensembl_column=encol,
            ),
        )
    else:
        oracle = co.load_hdf5(str(checkpoint))
    if not rawlink.exists():
        status(side, "get_links", genes=oracle.adata.n_vars)
        links = oracle.get_links(
            cluster_name_for_GRN_unit="state",
            alpha=params["celloracle_alpha"],
            bagging_number=params["bagging_n"],
            n_jobs=params["celloracle_model_jobs"],
            verbose_level=1,
        )
        links.to_hdf5(str(rawlink))
    else:
        links = co.load_hdf5(str(rawlink))
    raw = links.links_dict[state]
    raw.to_csv(OUT / f"tables/{side}_all_edges.tsv.gz", sep="\t", index=False)
    assert np.isfinite(raw[["coef_mean", "coef_abs"]]).all().all()
    assert raw.loc[raw.p.isna(), ["coef_mean", "coef_abs"]].eq(0).all().all()
    links.filter_links()
    kept = links.filtered_links[state]
    pd.testing.assert_frame_equal(
        kept, raw.loc[raw.p <= 0.001].sort_values("coef_abs", ascending=False).iloc[:10000]
    )
    assert np.isfinite(kept[["coef_mean", "coef_abs", "p"]]).all().all()
    kept.to_csv(OUT / f"tables/{side}_filtered_edges.tsv.gz", sep="\t", index=False)
    links.to_hdf5(str(OUT / f"objects/{side}_filtered.celloracle.links"))
    status(side, "fit_simulation")
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=params["celloracle_alpha"], use_cluster_specific_TFdict=True)
    coef = oracle.coef_matrix_per_cluster[state]
    assert np.isfinite(coef).all().all() and (np.diag(coef) == 0).all()
    assert coef.index.equals(oracle.adata.var_names) and coef.columns.equals(oracle.adata.var_names)
    active = set(oracle.active_regulatory_genes)
    assert active == set(coef.index[np.any(coef.values != 0, axis=1)])
    ledger = pd.read_csv(OUT / f"tables/{side}_feature_ledger.tsv.gz", sep="\t")
    ledger["celloracle_active"] = ledger.gene.isin(active)
    ledger["outgoing_nonzero"] = ledger.gene.map((coef != 0).sum(axis=1)).fillna(0).astype(int)
    ledger["KO_eligible"] = (
        ledger.is_TF.eq("Yes") & ledger.selected & ledger.prior_regulator & ledger.celloracle_active
    )
    ledger["reason"] = np.select(
        [
            ledger.KO_eligible,
            ~ledger.is_TF.eq("Yes"),
            ~ledger.in_UMI,
            ~ledger.nonzero_before_HVG,
            ~ledger.prior_regulator,
            ~ledger.selected,
            ~ledger.celloracle_active,
        ],
        [
            "eligible",
            "not_annotated_Yes",
            "absent_assay_symbol",
            "no_detected_UMI",
            "no_promoter_prior",
            "not_selected",
            "no_active_simulation_outgoing",
        ],
        default="unresolved",
    )
    ledger.loc[ledger.is_TF.eq("Yes") | ledger.prior_regulator].sort_values("gene").to_csv(
        OUT / f"tables/{side}_TF_eligibility.tsv", sep="\t", index=False
    )
    cov = pd.read_csv(OUT / f"tables/{side}_AxisA_coverage.tsv", sep="\t", dtype={"entrez_id": str})
    cov["incoming_nonzero"] = cov.gene.map((coef != 0).sum(axis=0)).fillna(0).astype(int)
    cov.to_csv(OUT / f"tables/{side}_AxisA_coverage.tsv", sep="\t", index=False)
    oracle.to_hdf5(str(ready))
    expected = (matsha(coef.values), matsha(oracle.adata.layers["imputed_count"]), oracle.adata.obs.copy())
    del oracle, links, coef
    gc.collect()
    saved = co.load_hdf5(str(ready))
    assert (
        matsha(saved.coef_matrix_per_cluster[state].values) == expected[0]
        and matsha(saved.adata.layers["imputed_count"]) == expected[1]
    )
    assert saved.adata.obs.equals(expected[2])
    for row in manifest.itertuples():
        assert sha(row.path) == row.sha256
    writej(
        OUT / f"qa/{side}_model_verified.json",
        dict(
            save_reload_exact=True,
            input_hashes_unchanged=True,
            celloracle_filter_exact=True,
            finite_coefficients=True,
            cells=n,
            donors=nd,
            genes=saved.adata.n_vars,
            active_regulators=len(active),
            eligible_TFs=int(ledger.KO_eligible.sum()),
            axis_measured=int(cov.in_model.sum()),
            axis_incoming=int((cov.incoming_nonzero > 0).sum()),
            ready_sha256=sha(ready),
        ),
    )
    status(side, "model_ready_verified", eligible_TFs=int(ledger.KO_eligible.sum()))


if __name__ == "__main__":
    try:
        run(sys.argv[1])
    except Exception as e:
        status(sys.argv[1], "failed", error=repr(e))
        raise
