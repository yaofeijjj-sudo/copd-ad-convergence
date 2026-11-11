from pathlib import Path
import sys

sys.path.insert(0, str(Path.cwd() / "R"))
from runtime import params
import sys, json, time, gc, subprocess, warnings, os, stat, errno

os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["VECLIB_MAXIMUM_THREADS"] = "8"
import numpy as np
import pandas as pd
import celloracle as co
from celloracle.trajectory.oracle_GRN import _do_simulation
from rebuild import OUT, ROOT, SPEC, sha, writej, status, matsha
from simulation_runtime import simulate_with_ood_diagnostic


def run(side):
    folder, state, _, n, nd = SPEC[side]
    assert co.__version__ == "0.20.0"
    acceleration_hash = None
    allowed_pairs = set()
    helper_hash = sha(Path(__file__).with_name("simulation_runtime.py"))
    from threadpoolctl import threadpool_info

    writej(
        OUT / f"qa/{side}_simulation_runtime.json",
        dict(
            BLAS_requested_threads=8,
            acceleration_check_sha256=acceleration_hash,
            execution_settings=dict(BLAS_threads=8, UCell_ncores=8),
            OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS"),
            threadpools=threadpool_info(),
        ),
    )
    freeze = {"protocol_sha256": sha(ROOT / "config/params.yml")}
    ready = OUT / f"objects/{side}_ready.celloracle.oracle"
    q = json.loads((OUT / f"qa/{side}_model_verified.json").read_text())
    assert sha(ready) == q["ready_sha256"]

    model = co.load_hdf5(str(ready))
    base = model.adata.layers["imputed_count"]
    genes = model.adata.var_names
    donors = model.adata.obs.donor.astype(str)
    dlist = sorted(donors.unique())
    didx = [np.flatnonzero(donors.eq(d).values) for d in dlist]
    assert len(dlist) == nd and model.adata.n_obs == n
    eligibility = pd.read_csv(OUT / f"tables/{side}_TF_eligibility.tsv", sep="\t")
    tfs = sorted(eligibility.loc[eligibility.KO_eligible, "gene"])
    assert tfs
    axis = pd.read_csv(OUT / f"tables/{side}_AxisA_coverage.tsv", sep="\t").gene.dropna().tolist()
    present = [g for g in axis if g in genes]
    aidx = genes.get_indexer(present)
    baseline_axis = base[:, aidx].copy()
    baseline_donor = np.array([baseline_axis[ix].mean(axis=0, dtype=np.float64) for ix in didx])
    od = OUT / "objects" / f"{side}_KO"
    od.mkdir(exist_ok=True)
    baseline_path = od / "baseline_AxisA.npz"
    if baseline_path.exists():
        with np.load(baseline_path) as archived:
            assert np.array_equal(archived["genes"], np.asarray(present))
            assert np.array_equal(archived["cells"], np.asarray(model.adata.obs_names, dtype=str))
            assert np.array_equal(archived["donors"], np.asarray(donors, dtype=str))
            assert np.array_equal(archived["expression"], baseline_axis)
    else:
        np.savez_compressed(
            baseline_path,
            genes=np.asarray(present),
            cells=np.asarray(model.adata.obs_names, dtype=str),
            donors=np.asarray(donors, dtype=str),
            expression=baseline_axis,
        )
    baseline_checkpoint_hash = sha(od / "baseline_AxisA.npz")
    tmp = ROOT / "results/cache" / f"celloracle_axisA_Entrez137_{side}"
    tmp.mkdir(exist_ok=True)
    pipe = tmp / "simulated.pipe"
    if not pipe.exists():
        os.mkfifo(pipe)
    assert stat.S_ISFIFO(pipe.stat().st_mode)
    pipe_keeper = os.open(pipe, os.O_RDWR | os.O_NONBLOCK)
    for name, values in [("genes", genes), ("cells", model.adata.obs_names), ("axis", axis)]:
        (tmp / f"{name}.txt").write_text("\n".join(values) + "\n")
    np.asarray(base, dtype="<f8").tofile(tmp / "baseline.bin")
    rlog = open(OUT / f"logs/{side}_ucell.log", "a")
    worker = subprocess.Popen(
        ["Rscript", "--vanilla", str(ROOT / "scripts/06_virtual_knockout/score_server.R"), str(tmp)],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=rlog,
        text=True,
        bufsize=1,
    )

    def response(expected):
        while True:
            line = worker.stdout.readline()
            if not line:
                raise RuntimeError(f"UCell worker exited {worker.poll()}; see log")
            if line.strip() == expected:
                return
            rlog.write(line)
            rlog.flush()

    response("READY")
    baseline_hash = matsha(base)
    coef = model.coef_matrix_per_cluster[state]
    model.calculate_randomized_coef_table(random_seed=params["celloracle_seed"])
    random_hash = matsha(model.coef_matrix_per_cluster_randomized[state].values)
    model.coef_matrix_per_cluster_randomized[state].to_csv(
        OUT / f"tables/{side}_randomized_coefficients.tsv.gz", sep="\t"
    )
    code_hash = sha(Path(__file__))
    scoring_hash = sha(ROOT / "scripts/06_virtual_knockout/score_server.R")
    stop_after = int(os.environ.get("CELLORACLE_STOP_AFTER_TFS", "0"))
    newly_completed = 0

    def donor_mean(x):
        return np.asarray([x[ix].mean(axis=0, dtype=np.float64) for ix in didx])

    try:
        for ti, tf in enumerate(tfs):
            if (OUT / "PAUSE_REQUEST").exists():
                status(side, "paused_at_TF_checkpoint", next_TF=tf, index=ti + 1, total=len(tfs))
                return
            dp = od / f"{tf}.complete.json"
            if dp.exists():
                record = json.loads(dp.read_text())
                assert record["protocol_sha256"] == freeze["protocol_sha256"]
                assert record["model_sha256"] == q["ready_sha256"]
                assert record["baseline_checkpoint_sha256"] == baseline_checkpoint_hash
                assert (
                    record["simulation_code_sha256"] == code_hash
                    and record["scoring_code_sha256"] == scoring_hash
                )
                if record["simulation_code_sha256"] == code_hash:
                    assert record["acceleration_check_sha256"] == acceleration_hash
                for p, h in record["outputs"].items():
                    assert sha(p) == h
                continue
            start = time.time()
            status(side, "KO", TF=tf, index=ti + 1, total=len(tfs))
            phases = {}
            phase_start = time.perf_counter()
            tfidx = genes.get_loc(tf)
            gem = pd.DataFrame(base, index=model.adata.obs_names, columns=genes)
            init = gem.copy()
            init[tf] = 0
            no_prop = _do_simulation(coef_matrix=coef, simulation_input=init, gem=gem, n_propagation=0)
            nonself = np.arange(len(genes)) != tfidx
            assert np.array_equal(no_prop.values[:, nonself], base[:, nonself])
            assert np.equal(no_prop[tf].values, 0).all()
            del init, no_prop, gem
            gc.collect()
            phases["no_propagation_seconds"] = time.perf_counter() - phase_start
            outputs = []
            info = {}
            score_baseline = None
            for mode, randomized in [("main", False), ("randomized", True)]:
                mode_start = time.perf_counter()
                with warnings.catch_warnings(record=True) as ws:
                    warnings.simplefilter("always")
                    oq = simulate_with_ood_diagnostic(
                        model,
                        perturb_condition={tf: 0.0},
                        n_propagation=params["propagation_n"],
                        ignore_warning=False,
                        use_randomized_GRN=randomized,
                        clip_delta_X=params["clip_delta_x"],
                    )
                simulation_seconds = time.perf_counter() - mode_start
                validation_start = time.perf_counter()
                sim = model.adata.layers["simulated_count"]
                delta = model.adata.layers["delta_X"]
                assert np.isfinite(sim).all() and (sim >= 0).all() and np.isfinite(delta).all()
                assert np.equal(sim[:, tfidx], 0).all()
                assert np.array_equal(delta, sim - base)
                oq.index.name = "gene"
                op = OUT / f"qa/{side}_{tf}_{mode}_range.tsv.gz"
                oq.to_csv(op, sep="\t")
                outputs.append(op)
                validation_seconds = time.perf_counter() - validation_start
                transport_start = time.perf_counter()
                worker.stdin.write(tf + "\t" + mode + "\n")
                worker.stdin.flush()
                buffer = memoryview(np.ascontiguousarray(sim, dtype="<f8")).cast("B")
                while True:
                    if worker.poll() is not None:
                        raise RuntimeError("UCell worker exited before matrix transport")
                    try:
                        fd = os.open(pipe, os.O_WRONLY | os.O_NONBLOCK)
                        os.set_blocking(fd, True)
                        break
                    except OSError as e:
                        if e.errno != errno.ENXIO:
                            raise
                        time.sleep(0.05)
                with os.fdopen(fd, "wb") as transport:
                    for start_byte in range(0, len(buffer), 8 * 1024 * 1024):
                        block = buffer[start_byte : start_byte + 8 * 1024 * 1024]
                        assert transport.write(block) == len(block)
                del buffer
                response("DONE")
                transport_scoring_seconds = time.perf_counter() - transport_start
                scoring_timing = json.loads((tmp / f"{tf}_{mode}_timing.json").read_text())
                export_start = time.perf_counter()
                if mode == "main":
                    score_baseline = np.fromfile(tmp / "baseline_score.bin", dtype="<f8")
                score_pred = np.fromfile(tmp / "simulated_score.bin", dtype="<f8")
                assert len(score_pred) == n and len(score_baseline) == n
                diff = score_pred - score_baseline
                scores = pd.DataFrame(
                    {
                        "donor": dlist,
                        "cells": [len(ix) for ix in didx],
                        "baseline_UCell": donor_mean(score_baseline),
                        "simulated_UCell": donor_mean(score_pred),
                        "delta_UCell": donor_mean(diff),
                    }
                )
                scores.insert(0, "TF", tf)
                scores.insert(0, "mode", mode)
                scores.insert(0, "side", side)
                scores["readout"] = "AxisA137_measured_target_excluded_UCell_delta"
                p = OUT / f"tables/{side}_{tf}_{mode}_donor_scores.tsv"
                scores.to_csv(p, sep="\t", index=False)
                outputs.append(p)
                sim_axis = sim[:, aidx]
                delta_axis = delta[:, aidx]
                gdelta = donor_mean(delta_axis)
                gsim = donor_mean(sim_axis)
                gtable = pd.DataFrame(
                    {
                        "donor": np.repeat(dlist, len(present)),
                        "gene": np.tile(present, len(dlist)),
                        "baseline": baseline_donor.ravel(),
                        "simulated": gsim.ravel(),
                        "delta": gdelta.ravel(),
                    }
                )
                gtable["forced_self"] = gtable.gene.eq(tf)
                gtable.insert(0, "TF", tf)
                gtable.insert(0, "mode", mode)
                gtable.insert(0, "side", side)
                p = OUT / f"tables/{side}_{tf}_{mode}_donor_genes.tsv.gz"
                gtable.to_csv(p, sep="\t", index=False)
                outputs.append(p)
                ap = od / f"{tf}_{mode}.npz"
                if mode == "main":
                    np.savez_compressed(
                        ap,
                        axis_genes=np.asarray(present),
                        delta_axis=delta_axis,
                        baseline_score=score_baseline,
                        simulated_score=score_pred,
                    )
                else:
                    np.savez_compressed(ap, baseline_score=score_baseline, simulated_score=score_pred)
                outputs.append(ap)
                downstream = oq.drop(index=tf, errors="ignore")
                affected = np.any(delta_axis[:, np.asarray(present) != tf] != 0, axis=0)
                info[mode] = dict(
                    warnings=[str(w.message) for w in ws],
                    n_axis_scored=len(present) - int(tf in present),
                    n_background=len(genes) - 1,
                    mean_donor_delta=float(scores.delta_UCell.mean()),
                    median_donor_delta=float(scores.delta_UCell.median()),
                    min_donor_delta=float(scores.delta_UCell.min()),
                    max_donor_delta=float(scores.delta_UCell.max()),
                    q25_donor_delta=float(scores.delta_UCell.quantile(0.25)),
                    q75_donor_delta=float(scores.delta_UCell.quantile(0.75)),
                    n_donors=nd,
                    axis_downstream_changed_genes=int(affected.sum()),
                    ood_gt50pct_downstream=int((downstream.Max_exceeding_ratio > 0.5).sum()),
                    range_undefined_genes=int(oq.isna().any(axis=1).sum()),
                    finite_nonnegative=True,
                    self_zero=True,
                )
                phases[mode] = dict(
                    simulation_and_OOD_seconds=simulation_seconds,
                    validation_and_diagnostic_export_seconds=validation_seconds,
                    transport_and_UCell_seconds=transport_scoring_seconds,
                    aggregation_and_checkpoint_export_seconds=time.perf_counter() - export_start,
                    UCell=scoring_timing,
                )
                del sim, delta, sim_axis, delta_axis, gtable
                for transient in ["simulated_count", "delta_X", "simulation_input"]:
                    if transient in model.adata.layers:
                        del model.adata.layers[transient]
                gc.collect()
            assert matsha(model.adata.layers["imputed_count"]) == baseline_hash
            writej(
                dp,
                dict(
                    TF=tf,
                    side=side,
                    status="complete",
                    protocol_sha256=freeze["protocol_sha256"],
                    model_sha256=q["ready_sha256"],
                    simulation_code_sha256=code_hash,
                    scoring_code_sha256=scoring_hash,
                    runtime_helper_sha256=helper_hash,
                    acceleration_check_sha256=acceleration_hash,
                    execution_settings=dict(BLAS_threads=8, UCell_ncores=8),
                    phase_seconds=phases,
                    baseline_checkpoint_sha256=baseline_checkpoint_hash,
                    randomized_coef_sha256=random_hash,
                    no_propagation_downstream_exact_zero=True,
                    no_propagation_target_excluded_score_inputs_identical=True,
                    summary=info,
                    seconds=time.time() - start,
                    outputs={str(p): sha(p) for p in outputs},
                ),
            )
            status(side, "KO_completed", TF=tf, index=ti + 1, total=len(tfs), seconds=time.time() - start)
            newly_completed += 1
            if stop_after and newly_completed >= stop_after:
                status(side, "execution_batch_checkpoint", completed_in_batch=newly_completed, total=len(tfs))
                return
    finally:
        if worker.poll() is None:
            worker.stdin.write("QUIT\n")
            worker.stdin.flush()
            worker.wait()
        os.close(pipe_keeper)
        rlog.close()
    assert len(list(od.glob("*.complete.json"))) == len(tfs)
    assert sha(ready) == q["ready_sha256"] and sha(OUT / "protocol/PROTOCOL.md") == freeze["protocol_sha256"]
    writej(
        OUT / f"qa/{side}_simulation_complete.json",
        dict(
            eligible=len(tfs),
            completed=len(tfs),
            all_main_and_randomized=True,
            model_unchanged=True,
            protocol_unchanged=True,
        ),
    )
    status(side, "all_KOs_and_controls_complete", total=len(tfs))


if __name__ == "__main__":
    try:
        run(sys.argv[1])
    except Exception as e:
        status(sys.argv[1], "failed", error=repr(e))
        raise
