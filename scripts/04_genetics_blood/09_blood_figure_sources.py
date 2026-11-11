from pathlib import Path
import shutil

out = Path.cwd() / "results/genetics_blood"
(out / "source_data").mkdir(parents=True, exist_ok=True)
shutil.copy2(out / "results/fgsea_six_results.csv", out / "source_data/blood_six.csv")
