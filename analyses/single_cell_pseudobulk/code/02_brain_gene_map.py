# Export symbol -> Ensembl map from the H5AD whose var_names index the 2026-09-25 brain count table.
import anndata as ad, pandas as pd
R=''
a=ad.read_h5ad(R+'03_candidate_analysis/Brain_CellOracle_paper_gene_retention_20260909/objects/three_states_UMIs.h5ad',backed='r')
v=pd.DataFrame({'gene':a.var_names.astype(str),'gene_id':a.var['gene_ids'].astype(str).to_numpy()})
assert v.gene.is_unique and v.gene_id.is_unique
c=pd.read_csv(R+'虚拟敲除/AxisA_疾病标尺核对_20260925/source_data/brain_donor_pseudobulk_counts.tsv.gz',sep='\t',usecols=['gene'])
assert c.gene.tolist()==v.gene.tolist()
v.to_csv(R+'03_candidate_analysis/AxisA_pseudobulk_disease_20261014/source_data/brain_symbol_to_ensembl.tsv',sep='\t',index=False)
print(len(v))
