# Extract the 249 Patir 2019 core microglia symbols (Table S3, column "Core signature/ overlap in two dataset" == "Core").
# Run from project root: _workspace/envs/celloracle/bin/python -I 03_candidate_analysis/AxisA_auxiliary_checks_20261019/code/03b_extract_patir.py
import pandas as pd
O = "03_candidate_analysis/AxisA_auxiliary_checks_20261019/"
x = pd.read_excel(O + "source_data/glia23572-sup-0001-tables3.xlsx", header=None).iloc[6:]
core = sorted(set(x[x[13] == "Core"][1].astype(str).str.strip())); assert len(core) == 249
open(O + "source_data/patir2019_core_symbols.txt", "w").write("\n".join(core) + "\n")
