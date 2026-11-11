# Recount DFC Micro-PVM_2 nuclei/donors under the CellOracle set (all DFC nuclei) and the CPS/pseudobulk filter
# (method 10Xv3.1, experiment_component_failed False, Used in analysis True, >=20 nuclei per donor).
# Run from project root: _workspace/envs/celloracle/bin/python -I 03_candidate_analysis/Methods_checks_20261019/code/02_brain_donor_recount.py
import pandas as pd
src = "cell_localization/analysis/06_single_cell_projection/results/SEAAD_DFC_axisA_localization_20260831_candidate/inputs/SEAAD_Immune_cell_metadata.tsv.gz"
m = pd.read_csv(src, sep="\t", usecols=["Donor ID", "Supertype", "method", "experiment_component_failed", "Used in analysis", "Brain Region"])
d = m[(m["Supertype"] == "Micro-PVM_2") & (m["Brain Region"] == "DFC")]
print("CellOracle set (all DFC Micro-PVM_2):", len(d), "nuclei,", d["Donor ID"].nunique(), "donors")
print(d.groupby(["method", "experiment_component_failed", "Used in analysis"]).size().to_string())
f = d[(d.method == "10Xv3.1") & (d["experiment_component_failed"].astype(str) == "False") & (d["Used in analysis"].astype(str) == "True")]
c = f.groupby("Donor ID").size(); keep = c[c >= 20]
print("CPS/pseudobulk set:", keep.sum(), "nuclei,", keep.size, "donors")
allc = d.groupby("Donor ID").size()
for k in sorted(set(allc.index) - set(keep.index)):
    print("donor not in CPS set:", k, "| all DFC Micro-PVM_2 nuclei:", allc[k], "| passing filter:", int(c.get(k, 0)))
