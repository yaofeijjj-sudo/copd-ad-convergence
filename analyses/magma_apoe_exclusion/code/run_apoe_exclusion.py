# -*- coding: utf-8 -*-
# APOE-region exclusion sensitivity for the current (SYN01, Axis A=137 / Axis B=434) AD MAGMA test.
# Region: APOE (Entrez 348) gene boundaries +/- 1 Mb, NCBI37.3 -- the definition pre-specified in the
# 2026-07-15 protocol (GWAS_and_blood/analysis/05_disease_specificity_gwas/.../preanalysis_protocol_20260715.md §8)
# and implemented in 44_run_current_axes_locus_influence_sensitivities_20260715.R.
#   V1 set_only     : region genes removed from the Axis sets; background unchanged (July implementation).
#   V2 whole_region : region genes removed from the whole analysis (MAGMA v1.10 --settings gene-exclude).
from pathlib import Path
import csv, json, os, subprocess, hashlib

R = Path('.')
O = R / '03_candidate_analysis/AxisA_MAGMA_APOE_exclusion_20261016'
SYN = R / 'GWAS_and_blood/SYN01_Entrez_downstream_20261001'
TOOL = R / '_workspace/project_tools/external_tools/magma_v1.10_official'
MAGMA = TOOL / 'magma_v1.10_mac/magma'
GENE_LOC = TOOL / 'NCBI37/NCBI37.3.gene.loc'
GENES_RAW = SYN / 'results/magma/AD_total_N.genes.raw'
BASE_SETS = SYN / 'results/magma_inputs/locked_axes.sets'
BASE_GSA = SYN / 'results/magma/AD_locked_axes.gsa.out'
WINDOW = 1_000_000

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

loc = {}
for line in GENE_LOC.open():
    e, chrom, s, t, strand, sym = line.rstrip('\n').split('\t')[:6]
    loc[e] = (chrom, int(s), int(t), sym)
apoe = loc['348']
assert apoe[0] == '19' and apoe[3] == 'APOE'
r_start, r_end = max(1, apoe[1] - WINDOW), apoe[2] + WINDOW
in_region = lambda e: e in loc and loc[e][0] == '19' and loc[e][2] >= r_start and loc[e][1] <= r_end

# genes analysed by MAGMA (genes.raw: non-comment rows, first field = Entrez ID)
analysed = [l.split()[0] for l in GENES_RAW.open() if not l.startswith('#')]
region_analysed = sorted((e for e in analysed if in_region(e)), key=lambda e: loc[e][1])

sets = {}
for line in BASE_SETS.open():
    f = line.split()
    sets[f[0]] = f[1:]

excl = O / 'inputs/APOE_1Mb_region_genes.entrez.txt'
excl.write_text(''.join(e + '\n' for e in region_analysed))
v1 = O / 'inputs/axes_no_apoe_1mb.sets'
with v1.open('w') as f:
    for k, ids in sets.items():
        f.write(k + '_no_apoe_1mb\t' + '\t'.join(e for e in ids if not in_region(e)) + '\n')

with (O / 'inputs/APOE_1Mb_region_gene_manifest.csv').open('w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['entrez_id', 'symbol', 'chr', 'start_GRCh37', 'end_GRCh37', 'in_AxisA', 'in_AxisB'])
    for e in region_analysed:
        c, s, t, sym = loc[e]
        w.writerow([e, sym, c, s, t, e in sets['AxisA'], e in sets['AxisB']])

env = os.environ.copy(); env['LC_ALL'] = 'C'
def run(args, log):
    with (O / 'logs' / log).open('w') as f:
        f.write(json.dumps([str(a) for a in args]) + '\n'); f.flush()
        subprocess.run([str(a) for a in args], stdout=f, stderr=subprocess.STDOUT, check=True, env=env)

run([MAGMA, '--gene-results', GENES_RAW, '--set-annot', v1, '--out', O / 'magma/AD_V1_set_only'], 'V1_set_only.log')
run([MAGMA, '--gene-results', GENES_RAW, '--set-annot', BASE_SETS, '--settings', f'gene-exclude={excl}',
     '--out', O / 'magma/AD_V2_whole_region'], 'V2_whole_region.log')

def read_gsa(p):
    rows = [l.split() for l in Path(p).open() if not l.startswith('#')]
    hdr = rows[0]
    return {r[0]: dict(zip(hdr, r)) for r in rows[1:]}

def holm(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i]); out = [0.0] * len(ps); run_max = 0.0
    for k, i in enumerate(order):
        run_max = max(run_max, min(1.0, (len(ps) - k) * ps[i])); out[i] = run_max
    return out

base = read_gsa(BASE_GSA)
variants = [('full_set_reference', base, lambda a: a),
            ('V1_set_only', read_gsa(O / 'magma/AD_V1_set_only.gsa.out'), lambda a: a + '_no_apoe_1mb'),
            ('V2_whole_region', read_gsa(O / 'magma/AD_V2_whole_region.gsa.out'), lambda a: a)]
out = []
for name, res, key in variants:
    rows = [res[key(a)] for a in ('AxisA', 'AxisB')]
    hp = holm([float(r['P']) for r in rows])
    for a, r, h in zip(('Axis A', 'Axis B'), rows, hp):
        out.append(dict(variant=name, axis=a, NGENES=r['NGENES'], BETA=r['BETA'], BETA_STD=r['BETA_STD'], SE=r['SE'],
                        P_one_sided=r['P'], holm_P_two_axes_within_variant=f'{h:.4g}',
                        n_region_genes_in_set=sum(in_region(e) for e in sets[a.replace(' ', '')] if e in set(analysed)),
                        background_genes=len(analysed) - (len(region_analysed) if name == 'V2_whole_region' else 0)))
with (O / 'results/AD_MAGMA_APOE_exclusion.csv').open('w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)

(O / 'results/provenance.json').write_text(json.dumps(dict(
    region=dict(build='GRCh37/NCBI37.3', chr='19', start=r_start, end=r_end, definition='APOE (Entrez 348) +/- 1 Mb',
                source='preanalysis_protocol_20260715.md §8; locus_exclusion_region_manifest_20260715.csv'),
    n_region_genes_analysed=len(region_analysed),
    inputs={str(p.relative_to(R)): sha(p) for p in (GENES_RAW, BASE_SETS, GENE_LOC)},
    magma=str(MAGMA.relative_to(R)), magma_sha256=sha(MAGMA)), indent=2, ensure_ascii=False))
print(open(O / 'results/AD_MAGMA_APOE_exclusion.csv').read())
