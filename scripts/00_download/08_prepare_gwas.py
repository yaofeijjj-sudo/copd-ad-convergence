from pathlib import Path
import os
import subprocess

root = Path.cwd()
out = root / "results/genetics_blood/results/magma_inputs"
out.mkdir(parents=True, exist_ok=True)
env = os.environ.copy()
env["LC_ALL"] = "C"
inputs = {
    "AD": root / "data/raw/gwas/ad_consensus_2026_gcst90704648_no_biobank/harmonised/GCST90704648.h.tsv.gz",
    "COPD": root / "data/raw/gwas/copd_ml_spirogram_2023_gcst90244098/harmonised/GCST90244098.h.tsv.gz",
}
for disease, raw in inputs.items():
    target = out / f"{disease}.magma.pval"
    if target.exists():
        continue
    valid = out / f"{disease}.sorted.tsv"
    extract = r"""BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$i);if($i=="rsid")s=i;if($i=="p_value")p=i;if($i=="Neff_total")ne=i;if($i=="N_cases")nc=i;if($i=="N_controls")nt=i}next} {id=$s;gsub(/\r/,"",id);pv=$p+0;if(id~/^rs[0-9]+$/&&pv>0&&pv<=1){if(disease=="AD"){if($ne+0<=0)next;total=($nc+0)+($nt+0);if(total<=0)exit 8;if(pv<1e-300)pv=1e-300;print id,pv,sprintf("%.0f",total)}else print id,$p,325027}}"""
    with valid.open("wb") as f:
        a = subprocess.Popen(["gzip", "-dc", str(raw)], stdout=subprocess.PIPE, env=env)
        b = subprocess.Popen(
            ["awk", "-v", f"disease={disease}", extract], stdin=a.stdout, stdout=subprocess.PIPE, env=env
        )
        a.stdout.close()
        c = subprocess.Popen(["sort", "-S", "1G", "-k1,1"], stdin=b.stdout, stdout=f, env=env)
        b.stdout.close()
        assert c.wait() == b.wait() == a.wait() == 0
    group = r"""function flush(){if(n==1)print saved;else print id > excluded} BEGIN{FS=OFS="\t";print "SNP","P","N"} NR==1{id=$1;saved=$0;n=1;next} $1==id{n++;next} {flush();id=$1;saved=$0;n=1} END{if(NR>0)flush()}"""
    with target.open("wb") as f:
        subprocess.run(
            ["awk", "-v", f'excluded={out / (disease + ".duplicate_rsids.txt")}', group, str(valid)],
            stdout=f,
            check=True,
            env=env,
        )
    valid.unlink()
