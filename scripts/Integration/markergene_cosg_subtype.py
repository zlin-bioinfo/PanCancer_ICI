import pandas as pd
import scanpy as sc 
import anndata as ad
import numpy as np
import cosg
import os
def read_10x_manually(dir_datset):
    adata = sc.read_mtx(filename= dir_datset + '/seu_final/matrix.mtx')
    adata_bc=pd.read_csv(filepath_or_buffer= dir_datset + '/seu_final/barcodes.tsv', header=None)
    adata_features=pd.read_csv(filepath_or_buffer= dir_datset + '/seu_final/features.tsv',header=None)
    adata_obs=pd.read_csv(filepath_or_buffer= dir_datset + '/seu_final/metadata.tsv', sep='\t')
    adata= adata.T
    adata.obs=adata_obs
    adata.obs.index=adata_bc[0].tolist()
    adata.var['gene_name']= adata_features[0].tolist()
    adata.var.index= adata.var['gene_name']
    adata = adata[adata.obs['celltype_main']=='Mono/macro', :]
    adata = adata[~adata.obs['celltype_r2'].isin(['Mast','pDC','cDC1','cDC2_CD1C','cDC2_IL1B','cDC2-ISG','cDC2_CXCL9','DC_LC-like','MoDC','mregDC','cDC2']), :]
    print(adata.obs['cohort'].unique())
    return adata
# dataset_names = ['HCC_Guo','HCC_Ma','RCC_Bi']
dataset_names = ['SKCM_Becker', 'SKCM_Plozniak', 'BCC_Yost', 
                'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
                'NSCLC_Yan', 
                'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo','HCC_Ma','RCC_Bi']
datasets = [read_10x_manually(dir_datset='/bigdata/zlin/PanCancer_ICI/data/' + name) for name in dataset_names]
adata = ad.concat(datasets, join='outer')
# Get expressed background genes
cell_types = adata.obs['celltype_r2'].unique()
# Directory to save results
output_dir = '/bigdata/zlin/PanCancer_ICI/tables/bg_genes_per_celltype/'
os.makedirs(output_dir, exist_ok=True)
# Process each cell type separately
for cell_type in cell_types:
    adata_sub = adata[adata.obs['celltype_r2'] == cell_type].copy()
    # Filter genes expressed in at least 20% of cells in this subtype
    sc.pp.filter_genes(adata_sub, min_cells=adata_sub.n_obs / 5)
    # Save expressed genes to a file
    output_path = f"{output_dir}/bg_{cell_type}.txt"
    adata_sub.var_names.to_series().to_csv(output_path, index=False, header=False)
    print(f"Saved: {output_path}")
# Get marker genes
adata.obs_names_make_unique()
adata.obs['cohort'] = adata.obs['cohort'].astype(str)
adata.obs['celltype_r2'] = adata.obs['celltype_r2'].astype(str)
# adata.obs.loc[adata.obs['cohort'].isin(['SCC_Yost','BCC_Yost']), 'cohort'] = 'BCCSCC_Yost'
# myeloids
# adata.obs.loc[adata.obs['celltype_r2'].isin(['cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'MoDC',]), 'celltype_r2'] = 'cDC2_CD1C'
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.raw = adata
groupby='celltype_r2'
cosg.cosg(
   adata,
   key_added='cosg',
   batch_key='cohort',
   batch_cell_number_threshold=50,
   mu=100,
   expressed_pct=0.1,
   remove_lowly_expressed=True,
   n_genes_user=200,
   groupby=groupby
)
# Convert to DataFrame
df = pd.DataFrame(adata.uns['cosg']['names'])
# Save as CSV
df.to_csv('/bigdata/zlin/PanCancer_ICI/tables/cosg_genes_monomac.csv', index=False)
