import pandas as pd
import scanpy as sc 
import anndata as ad
import numpy as np
import cosg

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
    print(adata.obs['cohort'].unique())
    return adata

# dataset_names = ['HCC_Guo','HCC_Ma','RCC_Bi']
dataset_names = ['SCC_Yost', 'NSCLC_Liu','SKCM_Becker', 'SKCM_Plozniak', 'BCC_Yost', 
                'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
                'NSCLC_Yan', 'NSCLC_Liu',
                'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo','HCC_Ma','RCC_Bi']
datasets = [read_10x_manually(dir_datset='/bigdata/zlin/PanCancer_ICI/data/' + name) for name in dataset_names]
adata = ad.concat(datasets, join='outer')
adata.obs_names_make_unique()
adata.obs['cohort'] = adata.obs['cohort'].astype(str)
adata.obs['celltype_r2'] = adata.obs['celltype_r2'].astype(str)
adata.obs.loc[adata.obs['cohort'].isin(['SCC_Yost','BCC_Yost']), 'cohort'] = 'BCCSCC_Yost'
adata.obs.loc[adata.obs['celltype_main'].isin(['Malignant(CNA+)','Epithelial(CNA+)']), 'celltype_r2'] = 'Epithelial(CNA+)'
adata.obs.loc[adata.obs['celltype_main']=='Melanocytes(CNA+)', 'celltype_r2'] = 'Melanocytes(CNA+)'
adata.obs.loc[adata.obs['celltype_main']=='Cycling', 'celltype_r2'] = 'Cycling'
adata.obs.loc[adata.obs['celltype_r2'].isin(['PC-cycling','GCB-cycling']), 'celltype_r2'] = 'Cycling'
adata = adata[~adata.obs['celltype_r2'].isin(['Melanocytes(CNA-)','Epithelial(CNA-)']), :]
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.raw = adata
n_gene=100
groupby='celltype_r2'
cosg.cosg(
   adata,
   key_added='cosg',
   batch_key='cohort',
   batch_cell_number_threshold=50,
   mu=100,
   expressed_pct=0.1,
   remove_lowly_expressed=True,
   n_genes_user=100,
   groupby=groupby
)
# Convert to DataFrame
df = pd.DataFrame(adata.uns['cosg']['names'])
# Save as CSV
df.to_csv('/bigdata/zlin/PanCancer_ICI/tables/cosg_genes.csv', index=False)
