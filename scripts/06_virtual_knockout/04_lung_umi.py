from pathlib import Path
import sys, json, hashlib, time, gc, tarfile, shutil, platform, warnings
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import celloracle as co

ROOT = Path(".")
OUT = ROOT / "results/virtual_knockout/input/lung_umi"
BASE = ROOT / "results/single_cell/preprocessing/lung"
H5DIR = ROOT / "data/raw/lung_h5"
TAR = ROOT / "data/raw/single_cell/copd_lung_gse302339/GSE302339_RAW.tar"
PRIOR = ROOT / "data/reference/promoter/hg38_TFinfo_dataframe_gimmemotifsv5_fpr2_threshold_10.parquet"
STATES = {
    "Alveolar_Mph_MT_positive": "Alveolar Mph MT-positive",
    "Interstitial_Mph_perivascular": "Interstitial Mph perivascular",
    "Alveolar_Mph_proliferating": "Alveolar Mph proliferating",
    "Alveolar_Mph_CCL3": "Alveolar Mph CCL3+",
    "Monocyte_derived_Mph": "Monocyte-derived Mph",
    "Alveolar_macrophages": "Alveolar macrophages",
}
EXPECTED = [118, 820, 1285, 3010, 16648, 64542]
for d in ["sources", "objects", "tables", "qa", "logs"]:
    (OUT / d).mkdir(parents=True, exist_ok=True)


def sha(path):
    with Path(path).open("rb") as handle:
        return sha_stream(handle)


def sha_stream(handle):
    h = hashlib.sha256()
    for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
        h.update(chunk)
    return h.hexdigest()


def status(stage, **kwargs):
    r = dict(stage=stage, time=time.strftime("%Y-%m-%d %H:%M:%S"), **kwargs)
    (OUT / "status.json").write_text(json.dumps(r, indent=2))
    print(json.dumps(r), flush=True)


def prepare():
    status("input_preparation")
    source_paths = [
        BASE / "objects/celltypist_all.tsv.gz",
        BASE / "objects/cell_qc_all.tsv.gz",
        BASE / "tables/sample_ledger.tsv",
    ]
    inputs = []
    ann = pd.read_csv(source_paths[0], sep="\t")
    ann = ann[ann.celltypist_label.isin(STATES.values())].copy()
    assert len(ann) == 86423 and not ann.duplicated(["gsm", "barcode"]).any()
    qc = pd.read_csv(source_paths[1], sep="\t", usecols=["gsm", "barcode", "passes_QC"])
    ann = ann.merge(qc, on=["gsm", "barcode"], validate="one_to_one")
    assert ann.passes_QC.all() and len(ann) == 86423
    ledger = pd.read_csv(source_paths[2], sep="\t")
    fields = [
        "gsm",
        "subject_id",
        "disease_group",
        "emphysema_group",
        "age",
        "sex",
        "smoke",
        "packyear",
        "gold",
        "fev1",
        "fev1_fvc",
        "lobe",
        "lobe_laa_percent",
    ]
    ann = ann.merge(ledger[fields], on="gsm", validate="many_to_one")
    ann["donor"] = ann.subject_id
    ann["state"] = ann.celltypist_label.map({v: k for k, v in STATES.items()})
    ann["cell_id"] = ann.gsm + "|" + ann.barcode
    assert ann.cell_id.is_unique and ann.donor.nunique() == 56
    assert ann.donor.notna().all()
    for field in ["disease_group", "age", "sex"]:
        assert ann.groupby("donor")[field].nunique(dropna=False).max() == 1, field
    ann = ann.sort_values(["gsm", "barcode"]).set_index("cell_id")
    ann.to_csv(OUT / "tables/cell_metadata.tsv.gz", sep="\t")
    ann.groupby(["state", "donor", "disease_group"]).size().rename("n_cells").to_csv(
        OUT / "tables/donor_state_counts.tsv", sep="\t"
    )
    ann.groupby(["state", "disease_group"]).agg(
        n_cells=("barcode", "size"), n_donors=("donor", "nunique")
    ).to_csv(OUT / "tables/clinical_state_coverage.tsv", sep="\t")
    state_counts = ann.groupby("state").agg(n_cells=("barcode", "size"), n_donors=("donor", "nunique"))
    for (state, label), n in zip(STATES.items(), EXPECTED):
        assert state_counts.loc[state, "n_cells"] == n
    state_counts.to_csv(OUT / "tables/input_state_counts.tsv", sep="\t")
    parts = []
    verified = []
    var_ref = None
    with tarfile.open(TAR) as tar:
        members = {Path(m.name).name: m for m in tar.getmembers() if m.isfile()}
        for i, gsm in enumerate(sorted(ann.gsm.unique()), 1):
            matches = list(H5DIR.glob(gsm + "_*.h5"))
            assert len(matches) == 1
            p = matches[0]
            digest = sha(p)
            with tar.extractfile(members[p.name]) as f:
                assert digest == sha_stream(f)
            inputs.append(
                dict(path=str(p), bytes=p.stat().st_size, sha256=digest, previous_manifest_match=True)
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                a = sc.read_10x_h5(p)
            original_names = a.var_names.copy()
            a.var_names_make_unique()
            if var_ref is None:
                var_ref = a.var.copy()
                pd.DataFrame(
                    {
                        "gene_symbol": original_names,
                        "model_name": a.var_names,
                        "gene_id": a.var.gene_ids.to_numpy(),
                    }
                ).to_csv(OUT / "tables/gene_name_map.tsv", sep="\t", index=False)
            else:
                pd.testing.assert_frame_equal(a.var, var_ref)
            meta = ann[ann.gsm.eq(gsm)]
            assert meta.barcode.isin(a.obs_names).all()
            a = a[meta.barcode].copy()
            values = a.X.data
            assert (
                np.isfinite(values).all() and (values >= 0).all() and np.equal(values, np.floor(values)).all()
            )
            a.obs = meta.copy()
            parts.append(a)
            verified.append(
                dict(
                    gsm=gsm,
                    n_cells=a.n_obs,
                    tar_member_hash_match=True,
                    QC_pass=True,
                    barcode_alignment=True,
                    finite_nonnegative_integer_UMIs=True,
                )
            )
            if i % 10 == 0:
                status("reading_verified_libraries", done=i, total=ann.gsm.nunique())
    a = ad.concat(parts, join="inner", merge="same")
    del parts
    gc.collect()
    assert a.shape == (86423, len(var_ref)) and a.obs_names.is_unique and set(a.obs_names) == set(ann.index)
    for state in STATES:
        b = a[a.obs.state.eq(state)].copy()
        b.obs["state"] = pd.Categorical(b.obs.state)
        b.write_h5ad(OUT / f"objects/{state}_UMIs.h5ad", compression="gzip")
    pd.DataFrame(verified).to_csv(OUT / "qa/library_input_checks.tsv", sep="\t", index=False)
    for p in [PRIOR, Path(__file__)]:
        inputs.append(dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p), previous_manifest_match=None))
    pd.DataFrame(inputs).to_csv(OUT / "qa/input_manifest.tsv", sep="\t", index=False)
    (OUT / "qa/software_versions.json").write_text(
        json.dumps(
            dict(
                python=sys.version,
                platform=platform.platform(),
                celloracle=co.__version__,
                scanpy=sc.__version__,
                anndata=ad.__version__,
                numpy=np.__version__,
                pandas=pd.__version__,
            ),
            indent=2,
        )
    )
    status("inputs_ready", cells=86423, donors=56, states=6)


if __name__ == "__main__":
    prepare()
