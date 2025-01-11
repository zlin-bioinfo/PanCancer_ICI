#%%
import os
import diopy
import numpy as np
import pandas as pd
import scanpy as sc 
import anndata as ad
import scrublet as scr
import scanpy.external as sce
import matplotlib.pyplot as plt
from scipy import sparse
#%%
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
#%%
def run_scr(adata):
    adata.obs['sample'] = adata.obs['sample'].astype('category')
    batches = adata.obs['sample'].cat.categories.tolist()
    alldata = {}
    for batch in batches:
        tmp = adata[adata.obs['sample'] == batch,]
        print(batch, ":", tmp.shape[0], " cells")
        scrub = scr.Scrublet(tmp.X)
        out = scrub.scrub_doublets(verbose=False, n_prin_comps = 20)
        alldata[batch] = pd.DataFrame({'doublet_score':out[0],'predicted_doublets':out[1]},index = tmp.obs.index)
        print(alldata[batch].predicted_doublets.sum(), " predicted_doublets")
    # add predictions to the adata object.
    scrub_pred = pd.concat(alldata.values())
    adata.obs['doublet_scores'] = scrub_pred['doublet_score'].values
    adata.obs['predicted_doublets'] = scrub_pred['predicted_doublets'].values
    # print("Predicted doublets %d \n Removed"%sum(adata.obs['predicted_doublets']))
    adata = adata[adata.obs["predicted_doublets"] == False, :]
    return adata
def noise_removal(adata):
    # noise genes removal (mito, HSP, ribo, dissociation, and MALAT1)
    noise = pd.read_excel("../data/41586_2022_5400_MOESM3_ESM.xlsx", sheet_name="1d_Genes excluded", index_col=0, header = 2)
    noise = noise.to_numpy().flatten()
    noise = noise[~pd.isnull(noise)]
    noise = np.append(noise, "MALAT1")
    remove = np.isin( adata.var_names.values.astype(str), noise)
    adata = adata[:,~remove]
    print("%d noise genes removed"%sum(remove))
    return adata
def cal_cc(adata):
    cell_cycle_genes = [x.strip() for x in open('../data/regev_lab_cell_cycle_genes.txt')]
    print("Including %d cell cycle genes." %len(cell_cycle_genes))
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    adata.obs["cell_cycle_diff"] = adata.obs["S_score"] - adata.obs["G2M_score"]
    print("Done")
    return adata
def pp_filt(adata):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, max_genes=6000)
    sc.pp.filter_genes(adata, min_cells=20)
    adata.var['mt'] = adata.var_names.str.startswith('MT-') 
    adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
    adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs['pct_counts_mt'] < 15, :]
    # adata = adata[adata.obs['pct_counts_hb'] < 2, :]
    # adata = run_scr(adata)
    adata = noise_removal(adata)
    return adata
def pp_norm(adata, dataset, write_h5 = False):
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    adata = cal_cc(adata)
    if write_h5 == True:
        diopy.output.write_h5(adata, file = '../data/logp1_' + dataset + '.h5', save_X=False)
    return adata
def pp_dd(adata, n_high_var=2000, resolution=0.5, batch_key='sample'):
    sc.pp.highly_variable_genes(adata, n_top_genes=n_high_var, subset=True, layer='counts', flavor='seurat_v3')
    sc.pp.regress_out(adata, ['total_counts'])
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, svd_solver='arpack')
    sce.pp.bbknn(adata, batch_key=batch_key)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added='leiden_' + str(resolution))
    return adata

#%%
# Becker dataset
# reading data
results_file = '/bigdata/zlin/Melanoma_meta/data/becker_raw.h5ad'
file_list = list()
sample_list = list()
adatas = list()
trt_group = ['Combicohort','Monocohort']
for trt in trt_group:
    folder_path = "/bigdata/zlin/Melanoma/data/Melanoma_Becker/" + trt
    filenames = os.listdir(folder_path)
    sample_list = sample_list + filenames
    filenames = [os.path.join(folder_path, d) for d in filenames]
    filenames = [filename + "/raw_feature_bc_matrix/" for filename in filenames]
    file_list = file_list + filenames
del file_list[sample_list.index('E19326')]
sample_list.remove('E19326')
for i in file_list:
    adata = sc.read_10x_mtx(i,cache=True)
    adata.obs["sample"] = sample_list[file_list.index(i)]
    adatas.append(adata)
adata = ad.concat(adatas)
adata.obs_names_make_unique()
# basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print(adata.n_obs, adata.n_vars)
#%%
# add metadata
anno = pd.read_csv('/bigdata/zlin/Melanoma/data/Melanoma_Becker/Metadata.csv', sep='\t', header=0)
anno = pd.merge(pd.DataFrame(data = adata.obs['sample']), anno,left_on='sample', right_on='Probe')
anno["Time point"] = anno["Time point"].apply(lambda x: 'On' if x == 'after' else 'Pre')
adata.obs["time_point"] = anno["Time point"].values
adata.obs["treatment"] = anno["Treatment"].values
adata.obs["patient"] = anno['Patient'].str.replace("atient ","").values
adata.obs["sample"] = adata.obs["patient"] + "_" + adata.obs["time_point"]
print(adata.obs['sample'].value_counts())
#%%
clinical = pd.read_excel('/bigdata/zlin/Melanoma_meta/data/metadata_becker.xlsx')
clinical = clinical.dropna(subset='Response')
clinical = pd.merge(pd.DataFrame(data = adata.obs['patient']), clinical,left_on='patient', right_on='Patient')
adata.obs['reponse']=clinical['Response']
#%%
adata.obs['site'] = clinical['Leision type(primary site/metastasis)'].values
adata.obs['dataset'] = 'SCKM_Becker'
adata.obs['platform'] = '10X'
adata.obs['expansion'] = 'n/a'
adata.obs['subtype'] = 'n/a'
# saving
adata.write(results_file)
#%%
adata = sc.read_h5ad("/bigdata/zlin/Melanoma_meta/data/becker_raw.h5ad")
num_cell = adata.obs['sample'].value_counts()
adata = adata[adata.obs['sample'].isin(num_cell[num_cell>200].index)]
adata = pp_filt(adata)
adata_norm = adata.copy()
adata_norm = pp_norm(adata_norm, 'becker', write_h5=True)
#%%
pred = pd.read_csv('../tables/singler_becker.csv', index_col=0)
adata_norm.obs['singler'] = pred['labels'].values
adata_norm = adata_norm[~adata_norm.obs['singler'].isin(['Keratinocytes','Epithelial cells','Erythrocytes'])]
adata=adata[adata.obs_names.isin(adata_norm.obs_names)]
#%%
adata_norm = pp_dd(adata_norm)
#%%
# check marker genes
marker_genes_dict = {
    'Mast cells':['TPSAB1','TPSB2','KIT'],
    'T_NK': ['CD3D','NKG7'],
    'B_cells':['CD79A','MS4A1'],
    'Plasma':['JCHAIN','MZB1'],
    'pDC':['IL3RA','LILRA4','PLD4'],
    'Myeloids':['LYZ','C1QA','CD68'],
    'Fibroblasts':['COL1A1','COL1A2','DCN'],
    'Endothelium':['CLDN5','PLVAP'],
    'Melanoma': ['MLANA','MITF']
}
for key, marker_list in marker_genes_dict.items():
    color_param = marker_list + ['leiden_0.5','singler']
    sc.pl.umap(adata_norm, color=color_param, title = marker_list + ['Leiden 0.5','singler'], show=False, frameon=False)
    plt.savefig("../figures/Becker/Markers_" + key + '.pdf')
    plt.close()
sc.pl.dotplot(adata_norm, marker_genes_dict, 'leiden_0.5',cmap='Blues')
sc.pl.umap(adata_norm,color='leiden_0.5', frameon=False, legend_loc='on data')
#%%
adata_norm.obs['major'] = 'Unresolved'
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["0","1","5"]), "major"] = "Melanoma"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["2","3"]), "major"] = "T/NK"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["4"]), "major"] = "Myeloids"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["7"]), "major"] = "Endothelium"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["6"]), "major"] = "Fibroblasts"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["8","9"]), "major"] = "B/Plasma"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["10"]), "major"] = "pDC"
# adjustment by singler
adata_norm.obs.loc[(adata_norm.obs["leiden_0.5"].isin(["0","1","2","4","5","7"])) & (adata_norm.obs["singler"].isin(["Fibroblasts"])), "major"] =  "Fibroblasts"
adata_norm.obs.loc[(adata_norm.obs["leiden_0.5"].isin(["0","1","2","3","5","7"])) & (adata_norm.obs["singler"].isin(["Macrophages"])), "major"] =  "Myeloids"
adata_norm.obs.loc[(adata_norm.obs["leiden_0.5"].isin(["0","4","5"])) & (adata_norm.obs["singler"].isin(["CD4+ T-cells","CD8+ T-cells","NK cells"])),"major"] =  "T/NK"
adata_norm.obs.loc[adata_norm.obs["singler"].isin(["Melanocytes"]),"major"] =  "Melanoma"
adata.obs['major']=adata_norm.obs['major']
adata.obs['celltype_orig']='n/a'
sc.pl.umap(adata_norm,color='major', frameon=False)
adata.obs['batch'] = adata.obs['dataset'].astype(str) + '_' + adata.obs['sample'].astype(str)
#%%
marker_genes_dict = {
    'pDC':['LILRA4'],
    'T/NK': ['NKG7','CD3D'],
    'Myeloids':['LYZ'],
    'Melanoma': ['MLANA'],
    'Fibroblasts':['COL1A2'],
    'Endothelium':['CLDN5'],
    'B/Plasma':['CD79A']
}
marker_list = []
for v in marker_genes_dict.values():
    marker_list.extend(v)
fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(16,6), gridspec_kw={'wspace':0.6})
sc.pl.umap(adata_norm, color=['major'],title="Becker",show=False, ax=ax1, frameon=False)
sc.pl.dotplot(adata_norm, marker_list, 'major',dendrogram=False, cmap='Blues', show=False, ax=ax2)
plt.savefig('../figures/Becker/Major.pdf')
plt.show()
# %%
adata.write("../data/becker.h5ad")
#%%
# Franken
adata = sc.read_h5ad("/bigdata/zlin/Melanoma/data/additional_datasets/Franken/Franken.h5ad")
#%%
num_cell = adata.obs['sample'].value_counts()
adata = adata[adata.obs['sample'].isin(num_cell[num_cell>200].index)]
adata = pp_filt(adata)
adata_norm = adata.copy()
adata_norm = pp_norm(adata_norm, 'franken', write_h5=True)
#%%
pred = pd.read_csv('../tables/singler_franken.csv', index_col=0)
adata_norm.obs['singler'] = pred['labels'].values
adata_norm = adata_norm[~adata_norm.obs['singler'].isin(['Melanocytes','Erythrocytes'])]
adata=adata[adata.obs_names.isin(adata_norm.obs_names)]
#%%
adata_norm = pp_dd(adata_norm)
#%%
# check marker genes
marker_genes_dict = {
    'Mast cells':['TPSAB1','TPSB2','KIT'],
    'T_NK': ['CD3D','NKG7'],
    'B_cells':['CD79A','MS4A1'],
    'Plasma':['JCHAIN','MZB1'],
    'pDC':['IL3RA','LILRA4','PLD4'],
    'Myeloids':['LYZ','CD68','CSF3R','FCGR2A','CD1C'],
    'Fibroblasts':['COL1A1','COL1A2','DCN'],
    'Endothelium':['CLDN5','PLVAP'],
    'Epithelium': ['KRT19','EPCAM','KRT17']
}
for key, marker_list in marker_genes_dict.items():
    color_param = marker_list + ['leiden_0.8','singler']
    sc.pl.umap(adata_norm, color=color_param, title = marker_list + ['Leiden 0.8','singler'], show=False, frameon=False)
    plt.savefig("../figures/Franken/Markers_" + key + '.pdf')
    plt.close()
sc.pl.dotplot(adata_norm, marker_genes_dict, 'leiden_0.8',cmap='Blues')
sc.pl.umap(adata_norm,color='leiden_0.8', frameon=False, legend_loc='on data')
#%%
sc.tl.rank_genes_groups(adata_norm, 'leiden_0.8', method='wilcoxon')
sc.pl.rank_genes_groups(adata_norm, n_genes=25, sharey=False)
pd.DataFrame(adata_norm.uns['rank_genes_groups']['names']).head(5)
#%%
adata_norm.obs['major'] = 'Unresolved'
adata_norm.obs.loc[adata_norm.obs["leiden_0.8"].isin(["0","18"]), "major"] = "Epithelium"
adata_norm.obs.loc[adata_norm.obs["leiden_0.8"].isin(["2","5","6","8","15","17"]), "major"] = "T/NK"
adata_norm.obs.loc[adata_norm.obs["leiden_0.8"].isin(["3","4","7","9","13"]), "major"] = "Myeloids"
adata_norm.obs.loc[adata_norm.obs["leiden_0.8"].isin(["11","21"]), "major"] = "Endothelium"
adata_norm.obs.loc[adata_norm.obs["leiden_0.8"].isin(["1","14"]), "major"] = "Fibroblasts"
adata_norm.obs.loc[adata_norm.obs["leiden_0.8"].isin(["10","12"]), "major"] = "B/Plasma"
adata_norm.obs.loc[adata_norm.obs["leiden_0.8"].isin(["16"]), "major"] = "pDC"
adata_norm.obs.loc[adata_norm.obs["leiden_0.8"].isin(["19"]), "major"] = "Mast_cells"
#%%
# adjustment by singler
# adata_norm.obs.loc[adata_norm.obs["leiden_0.8"]=='20', "major"] = adata_norm.obs.loc[adata_norm.obs["leiden_0.8"]=='20', "singler"]
adata.obs['major']=adata_norm.obs['major']
adata.obs['celltype_orig']='n/a'
sc.pl.umap(adata_norm,color='major', frameon=False)
adata.obs['batch'] = adata.obs['dataset'].astype(str) + '_' + adata.obs['sample'].astype(str)
#%%
del adata.raw
adata.write('../data/franken.h5ad')
# %%
# Bassez cohort 1
adata = sc.read_h5ad("/bigdata/zlin/Melanoma/data/additional_datasets/Bassez/Bassez_1_raw.h5ad")
#%%
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, groupby='sample')
#%%
adata = pp_filt(adata)
#%%
adata = pp_norm(adata, 'Bassez_1')
adata = pp_dd(adata)
sc.pl.umap(adata, color=['cellType'], frameon=False)

# #%%
# # check marker genes
# marker_genes_dict = {
#     'T_NK': ['CD3D','NKG7'],
#     'B_cells':['CD79A','MS4A1'],
#     'Plasma':['JCHAIN','MZB1'],
#     'pDC':['LILRA4','PLD4'],
#     'Myeloids':['S100A9','LYZ'],
#     'Fibroblasts':['COL1A2','SRPX','COL1A1','COL1A2'],
#     'Endothelium':['CLDN5','PLVAP'],
#     'Epithelium': ['KRT19','EPCAM','KRT17'],
#     'Mast_cells':['TPSAB1','TPSB2']
# }
# for key, marker_list in marker_genes_dict.items():
#     color_param = marker_list + ['leiden_0.5','cellType']
#     sc.pl.umap(adata, color=color_param, title = marker_list + ['Leiden 0.5','celltype_orig'], show=False, frameon=False)
#     plt.savefig("../figures/Bassez_1/Markers_" + key + '.pdf')
#     plt.close()
#%%
# adjust labels
adata.obs['dataset']='BRCA_Bazzes_1'
adata.obs.rename(columns={"cellType":"celltype_orig"}, inplace=True)
adata.obs['major'] = adata.obs['celltype_orig'].astype(object)
adata.obs.loc[adata.obs["major"].isin(["T_cell"]), "major"] = "T/NK"
adata.obs.loc[adata.obs["major"].isin(["B_cell"]), "major"] = "B/Plasma"
adata.obs.loc[adata.obs["major"].isin(["Endothelial_cell"]), "major"] = "Endothelium"
adata.obs.loc[adata.obs["major"].isin(["Myeloid_cell"]), "major"] = "Myeloids"
adata.obs.loc[adata.obs["major"].isin(["Fibroblast"]), "major"] = "Fibroblasts"
adata.obs.loc[adata.obs["major"].isin(["Mast_cell"]), "major"] = "Mast_cells"
adata.obs.loc[adata.obs["major"].isin(["Cancer_cell"]), "major"] = "Malignant"
# adata.obs.loc[(adata.obs["leiden_0.5"].isin(["8"])) & (adata.obs["cellType"].isin(["B_cell"])),"major"] =  "Plasma"
# adata.obs.loc[adata.obs["major"].isin(["B_cell"]), "major"] = "B_cells"
adata.obs['batch'] = adata.obs['dataset'].astype(str) +'_'+ adata.obs['sample'].astype(str)
#%%
# marker_genes_dict = {
#     'pDC':['LILRA4'],
#     'T/NK': ['CD3D','NKG7'],
#     'Plasma':['JCHAIN'],
#     'Myeloids':['LYZ'],
#     'Mast_cells':['TPSAB1'],
#     'Malignant': ['KRT19'],
#     'Fibroblasts':['COL1A2'],
#     'Endothelium':['CLDN5'],
#     'B_cells':['MS4A1']
# }
# marker_list = []
# for v in marker_genes_dict.values():
#     marker_list.extend(v)
# fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(16,6), gridspec_kw={'wspace':0.6})
# sc.pl.umap(adata, color=['major'],title="Bassez_1",show=False, ax=ax1, frameon=False)
# sc.pl.dotplot(adata, marker_list, 'major',dendrogram=False, cmap='Blues', show=False, ax=ax2)
# plt.savefig('../figures/Bassez_1/Major.pdf')
# plt.show()
#%%
del adata.raw
adata.write("../data/bassez_1.h5ad")
# %%
# Bassez cohort 2
adata = sc.read_h5ad("/bigdata/zlin/Melanoma/data/additional_datasets/Bassez/Bassez_2_raw.h5ad")
adata.obs['dataset']='BRCA_Bazzes_2'
adata = pp_filt(adata)
#%%
adata.obs.rename(columns={"cellType":"celltype_orig"}, inplace=True)
adata.obs['major'] = adata.obs['celltype_orig'].astype(object)
adata.obs.loc[adata.obs["major"].isin(["T_cell"]), "major"] = "T/NK"
adata.obs.loc[adata.obs["major"].isin(["B_cell"]), "major"] = "B/Plasma"
adata.obs.loc[adata.obs["major"].isin(["Endothelial_cell"]), "major"] = "Endothelium"
adata.obs.loc[adata.obs["major"].isin(["Myeloid_cell"]), "major"] = "Myeloids"
adata.obs.loc[adata.obs["major"].isin(["Fibroblast"]), "major"] = "Fibroblasts"
adata.obs.loc[adata.obs["major"].isin(["Mast_cell"]), "major"] = "Mast_cells"
adata.obs.loc[adata.obs["major"].isin(["Cancer_cell"]), "major"] = "Malignant"
adata.obs['batch'] = adata.obs['dataset'] +'_'+ adata.obs['sample']
#%%
del adata.raw
adata.write("../data/bassez_2.h5ad")
#%%
# adata = sc.read_h5ad("../data/bassez_1.h5ad")
# adata.uns['log1p']["base"] = None
# bassez_1 = adata.raw.to_data()
#%%
# Yost cohort
adata_bcc = sc.read_text('/bigdata/zlin/Melanoma/data/additional_datasets/Yost/GSE123813_bcc_scRNA_counts.txt.gz').transpose()
meta_bcc = pd.read_table('/bigdata/zlin/Melanoma/data/additional_datasets/Yost/GSE123813_bcc_all_metadata.txt.gz',index_col=0).reindex(adata_bcc.obs_names.values)
adata_bcc.obs = meta_bcc
adata_bcc.obs['subtype'] = 'BCC'
# adata_scc = sc.read_text('/bigdata/zlin/Melanoma/data/additional_datasets/Yost/GSE123813_scc_scRNA_counts.txt.gz').transpose()
# meta_scc = pd.read_table('/bigdata/zlin/Melanoma/data/additional_datasets/Yost/GSE123813_scc_metadata.txt.gz',index_col=0).reindex(adata_scc.obs_names.values)
# adata_scc.obs = meta_scc
# adata_scc = adata_scc[~adata_scc.obs['patient'].isin(['su013'])]
# adata_scc.obs.loc[adata_scc.obs["patient"].isin(["su010"]), "patient"] = "su010-S"
# adata_scc.obs['subtype'] = 'SCC'
# adata = ad.concat([adata_bcc, adata_scc])
adata = adata_bcc.copy()
adata.X = sparse.csr_matrix(adata.X)
#%% 
del adata.obs['UMAP1']
del adata.obs['UMAP2']
adata.obs.rename(columns={"treatment":"time_point","cluster":"celltype_orig"}, inplace=True)
adata.obs["sample"] = adata.obs["patient"].astype(str) + "_" + adata.obs["time_point"].str.capitalize()
adata.obs["sample"] = adata.obs["sample"].astype('object')
adata.obs["time_point"] = adata.obs["time_point"].apply(lambda x: 'Post' if x == 'post' else 'Pre')
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'],
             jitter=0.4, multi_panel=True, groupby='sample')
pp_filt(adata)
#%%
adata.obs["major"] = adata.obs["celltype_orig"].astype(str)
adata.obs.loc[adata.obs["major"].str.contains("B_cells") | adata.obs["major"].str.contains("Plasma"), "major"] = "B/Plasma"
adata.obs.loc[adata.obs["major"].str.contains("Tumor") , "major"] = "Malignant"
adata.obs.loc[adata.obs["major"].str.contains("CAFs") | adata.obs["major"].str.contains("Myofibroblasts") , "major"] = "Fibroblasts"
adata.obs.loc[adata.obs["major"].str.contains("Macrophages") , "major"] = "Myeloids"
adata.obs.loc[adata.obs["major"]=="DCs" , "major"] = "Myeloids"
adata.obs.loc[adata.obs["major"]=="pDCs" , "major"] = "pDC"
adata.obs.loc[adata.obs["major"].str.contains("Endothelial") , "major"] = "Endothelium"
adata.obs.loc[adata.obs["major"].str.contains("T_cells") | adata.obs["major"].str.contains("NK_cells") | adata.obs["major"].str.contains("Naive") | adata.obs["major"].str.contains("T") | adata.obs["major"].str.contains("CD8"), "major"] = "T/NK"
adata.obs['cohort'] ='anti-PD1'
adata.obs['platform'] = '10X'
adata.obs['expansion'] = 'n/a'
adata.obs["sample"] = adata.obs["patient"].astype('string') + "_" + adata.obs["time_point"].str.capitalize()
adata.obs["sample"] = adata.obs["sample"].astype('object')
clinical = pd.read_excel("/bigdata/zlin/Melanoma/data/additional_datasets/Yost/41591_2019_522_MOESM2_ESM.xlsx", sheet_name="SuppTable1", skiprows=range(0,3))
clinical = pd.merge(pd.DataFrame(data = adata.obs['patient']), clinical,left_on='patient', right_on='Patient')
adata.obs['response'] = clinical['Response'].values
adata.obs.loc[adata.obs["response"].str.contains("Yes") , "response"] = "Yes"
adata.obs['dataset'] ='BCC_Yost'
adata.obs['batch'] = adata.obs['dataset'] + '_' + adata.obs['sample']
#%%
adata.write("../data/yost.h5ad")
#%%
adata = sc.read_h5ad("/bigdata/zlin/Melanoma/data/additional_datasets/Franken/Franken.h5ad")
adata.obs['dataset']='HNSC_Franken'
num_cell = adata.obs['sample'].value_counts()
adata = adata[adata.obs['sample'].isin(num_cell[num_cell>200].index)]
adata = pp_filt(adata)
adata_norm = adata.copy()
adata_norm = pp_norm(adata_norm, 'franken', write_h5=True)
adata_norm = pp_dd(adata_norm)
#%%
# check marker genes
marker_genes_dict = {
    'Mast cells':['TPSAB1','TPSB2'],
    'T_NK': ['CD3D','NKG7'],
    'B_cells':['CD79A','MS4A1'],
    'Plasma':['JCHAIN','MZB1'],
    'pDC':['IL3RA','LILRA4','PLD4'],
    'Myeloids':['LYZ','C1QA','CD68'],
    'Fibroblasts':['COL1A1','COL1A2','DCN'],
    'Endothelium':['CLDN5','PLVAP'],
    'Epithelium': ['EPCAM','KRT14','KRT19']
}
for key, marker_list in marker_genes_dict.items():
    color_param = marker_list + ['leiden_0.5']
    sc.pl.umap(adata_norm, color=color_param, title = marker_list + ['Leiden 0.5'], show=False, frameon=False)
    plt.savefig("../figures/Franken/Markers_" + key + '.pdf')
    plt.close()
sc.pl.dotplot(adata_norm, marker_genes_dict, 'leiden_0.5',cmap='Blues')
sc.pl.umap(adata_norm,color='leiden_0.5', frameon=False, legend_loc='on data')
# %%
sc.tl.rank_genes_groups(adata_norm, 'leiden_0.5', method='wilcoxon')
sc.pl.rank_genes_groups(adata_norm, n_genes=25, sharey=False)
#%%
adata_norm.obs['major'] = 'Unresolved'
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["2","15"]), "major"] = "Epithelium"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["1","3","10","13"]), "major"] = "T/NK"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["0","4","12"]), "major"] = "Myeloids"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["7"]), "major"] = "Endothelium"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["5","9"]), "major"] = "Fibroblasts"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["6","8"]), "major"] = "B/Plasma"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["11"]), "major"] = "pDC"
adata_norm.obs.loc[adata_norm.obs["leiden_0.5"].isin(["14"]), "major"] = "Mast_cells"
adata.obs['major']=adata_norm.obs['major']
adata.obs['celltype_orig']='n/a'
sc.pl.umap(adata_norm,color='major', frameon=False)
#%%
del adata.raw
adata.write("../data/franken.h5ad")
# %%
# Arnon
adata = sc.read_csv('/bigdata/zlin/Melanoma/data/additional_datasets/Jerby-Arnon/GSE115978_counts.csv.gz').transpose()
meta = pd.read_csv('/bigdata/zlin/Melanoma/data/additional_datasets/Jerby-Arnon/GSE115978_cell.annotations.csv',index_col=0).reindex(adata.obs_names.values)
adata.obs = meta
adata.X = sparse.csr_matrix(adata.X)
#%%
adata = adata[~adata.obs['cell.types'].isin(['?'])]
adata.obs.rename(columns={"samples": "sample","cell.types":"celltype_orig","treatment.group":"treatment"}, inplace=True)
del adata.obs['no.of.genes']
del adata.obs['no.of.reads']
adata.obs['dataset'] = 'Mel_Arnon'
adata.obs['batch'] = adata.obs['dataset'] + '_' + adata.obs['sample']
adata.obs['cohort'] ='anti-PD1'
adata.obs['platform'] = 'Smart-seq2'
adata.obs['expansion'] = 'n/a'
adata.obs['response'] = 'n/a'
#%%
adata.obs["major"] = adata.obs['celltype_orig'].astype('object')
adata.obs.loc[adata.obs["major"].isin(["B.cell"]), "major"] = "B/Plasma"
adata.obs.loc[adata.obs["major"].isin(["Mal"]), "major"] = "Malignant"
adata.obs.loc[adata.obs["major"].isin(["Endo."]), "major"] = "Endothelium"
adata.obs.loc[adata.obs["major"].isin(["CAF"]), "major"] = "Fibroblasts"
adata.obs.loc[adata.obs["major"].isin(["Macrophage"]), "major"] = "Myeloids"
adata.obs.loc[adata.obs["major"].isin(["T.CD8","T.CD4","T.cell","NK"]) , "major"] = "T/NK"
#%%
adata.write("../data/arnon.h5ad")
# %%
