# Preprocessing
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
def preprocessing(adata):
    adata.layers['counts'] = adata.X.copy()
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True, layer='counts', flavor='seurat_v3')
    sc.pp.regress_out(adata, ['total_counts','S.Score','G2M.Score'])
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, svd_solver='arpack')
    sce.pp.bbknn(adata, batch_key='dataset')
    sc.tl.umap()
    sc.tl.leiden(adata, resolution=0.8, key_added='leiden_0.8')
    return adata

# Harmony
sce.pp.harmony(adata, key='batch)
sc.tl.umap(adata, neighbors_key='X_pca_harmony')
sc.tl.leiden(adata, resolution=0.8)

# check markers
marker_genes_dict = {
    'Mast cell':['TPSAB1','TPSB2'],
    'T/NK': ['CD3D','NKG7'],
    'B_cell':['CD79A','MS4A1'],
    'Plasma':['JCHAIN','MZB1'],
    'pDC':['LILRA4','PLD4'],
    'Myeloid':['LYZ','C1QA'],
    'Fibroblast':['COL1A1','COL1A2','DCN'],
    'Endothelium':['CLDN5','PLVAP','VWF'],
    'Epithelium':['KRT14','KRT17'],
    'Melanoma': ['MLANA','MITF','S100B']
}
sc.pl.dotplot(adata, marker_genes_dict, 'leiden',cmap='Blues')

# DEG
sc.tl.rank_genes_groups(endos[2], 'leiden_0.8', method='wilcoxon')
pd.DataFrame(endos.uns['rank_genes_groups']['names']).head(10)

# COSG
n_gene=30
groupby='leiden_1'
cosg.cosg(adata,
    key_added='cosg',
    # use_raw=False, layer='log1p', ## e.g., if you want to use the log1p layer in adata
    mu=100,
    expressed_pct=0.1,
    remove_lowly_expressed=True,
    n_genes_user=100,
    groupby=groupby)
marker_gene=pd.DataFrame(adata.uns['cosg']['names'])
marker_gene.head()

# annotate cells
adata.obs['celltype_major'] = 'Unresolved'
adata.obs.loc[adata.obs["leiden"].isin(["0","1","2","10"]) & adata.obs["celltype_sc"].str.contains("Mal"), "celltype_major"] = "Melanoma"
adata.obs.loc[adata.obs["leiden"].isin(["3","4"]) & adata.obs["celltype_bped_main"].isin(["CD4+ T-cells","CD8+ T-cells"]), "celltype_major"] = "T_cells"
adata.obs.loc[adata.obs["leiden"].isin(["3"]) & adata.obs["celltype_bped_main"].isin(["NK cells"]), "celltype_major"] = "NK_cells"
adata.obs.loc[adata.obs["leiden"].isin(["5"]) & adata.obs["celltype_sc"].isin(["Macrophage"]), "celltype_major"] = "Myeloids"
adata.obs.loc[adata.obs["leiden"].isin(["6"]) & adata.obs["celltype_sc"].str.contains("Endo."), "celltype_major"] = "Endothelium"
adata.obs.loc[adata.obs["leiden"].isin(["7"]) & adata.obs["celltype_sc"].str.contains("CAF"), "celltype_major"] = "Fibroblasts"
adata.obs.loc[adata.obs["leiden"].isin(["8"]) & adata.obs["celltype_sc"].str.contains("B.cell"), "celltype_major"] = "B_cells"
adata.obs.loc[adata.obs["leiden"].isin(["9"])& adata.obs["celltype_bped_main"].isin(["Epithelial cells","Keratinocytes"]), "celltype_major"] = "Epithelium"
adata.obs.loc[adata.obs["leiden"].isin(["11"]), "celltype_major"] = "Plasma"
adata.obs.loc[adata.obs["leiden"].isin(["12"]), "celltype_major"] = "pDC"

# remove ambiguous cells
adata = adata[~adata.obs['celltype_major'].isin(['Unresolved'])]
adata.obs['celltype_major'].value_counts()

# rearrange cell types
new_categories = ['T_cells', 'NK_cells', 'B_cells', 'Plasma', 'pDC', 'Myeloids','Fibroblasts','Endothelium', 'Epithelium', 'Melanoma']
new_dtype = pd.CategoricalDtype(categories=new_categories, ordered=True)
adata.obs['celltype_major'] = adata.obs['celltype_major'].astype(new_dtype)

marker_genes_dict = {
    'T/NK': ['CD3D','NKG7','GNLY'],
    'B':['CD79A','MS4A1'],
    'Plasma':['JCHAIN'],
    'pDC':['LILRA4'],
    'Myeloid':['LYZ','C1QA'],
    'Fibroblast':['COL1A1','COL1A2'],
    'Endothelium':['PLVAP','VWF'],
    'Epithelium':['KRT14','KRT17'],
    'Melanoma': ['MLANA','S100B']
}
sc.pl.dotplot(adata, marker_genes_dict,'celltype_major',cmap='Blues', var_group_rotation=60)

# save object
diopy.output.write_h5(adata_counts, file = '/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/processed_major_after.h5')

# remove unmatched samples
def filter_patients(adata):
    # Create the crosstab
    ct = pd.crosstab(adata.obs['patient'], adata.obs['time_point'])
    filtered = ct[(ct.iloc[:, 0] < 100) | (ct.iloc[:, 1] < 100)]
    return filtered.index.tolist()

rm_pt = {}
for i in range(len(datasets)):
    rm = filter_patients(datasets[i])
    rm_pt[datasets_names[i]] = rm
rm_pt['BCC_Yost'] = rm_pt['BCC_Yost'] 
print(rm_pt)

# set_figure_params for a single figure
with plt.rc_context({"figure.figsize": (7, 4), "figure.dpi": (80)}):
    sc.pl.violin(endo, groupby = 'leiden_1.5', keys=['lymphatic_score','EndMT_score'])

# shuffle the order of cells
np.random.seed(100)
random_order = np.random.permutation(list(range(nk.n_obs)))
sc.pl.umap(nk[random_order,:], color='celltype_r2')

# Endothelial cells
# Arterial Endothelial Cells
arterial_genes = ["IGFBP3","FN1","DLL4", "HEY2", "EFNB2", "NOTCH1", "ROBO4","SEMA3G"]
arterial_score_name = 'arterial_score'
sc.tl.score_genes(endo, gene_list=arterial_genes, score_name=arterial_score_name)

# Capillary Endothelial Cells
capillary_genes = ["PECAM1", "CD34", "VWF", "CLDN5","COL15A1"]
capillary_score_name = 'capillary_score'
sc.tl.score_genes(endo, gene_list=capillary_genes, score_name=capillary_score_name)

# Venous Endothelial Cells
venous_genes = ["FLT4", "ACKR1", "KDR", "NR2F2","SELE","CYP1B1"]
venous_score_name = 'venous_score'
sc.tl.score_genes(endo, gene_list=venous_genes, score_name=venous_score_name)

# Lymphatic Endothelial Cells
lymphatic_genes = ["PROX1", "LYVE1", "PDPN","FLT4","VEGFR3","CCL21"]
lymphatic_score_name = 'lymphatic_score'
sc.tl.score_genes(endo, gene_list=lymphatic_genes, score_name=lymphatic_score_name)

# EndMT Endothelial Cells
endmt_genes = ['TWIST1','SNAI1','SNAI2','ZEB2','COL1A1','COL1A2','COL3A1','COL6A1','MMP9','ACTA2','S100A4','CD44','CDH2', 'CDH11']
endmt_score_name = 'EndMT_score'
sc.tl.score_genes(endo, gene_list=endmt_genes, score_name=endmt_score_name)

# CD8
marker_genes_dict = {
    'Naive':['TCF7','LEF1','CCR7','SELL','MAL'],
    'Central memory':['IL7R','GPR183','ZFP36L2','CXCR4'],
    'Resident memory':['ZNF683','CD52','HOPX','ID2','CXCR6','XCL1'],
    'Terminal differentiation':['TBX21','ASCL2','CX3CR1','KLRG1'],
    'NK-like':['KLRD1','TYROBP','KIR2DL3','KIR2DL1','KIR3DL1','KIR3DL2','CD160','EOMES','TXK','KLRC1','KIR2DL4'],
    'Effector':['GZMK','CXCR5','CCR4','CD28','CXCR3','GZMH','CD27','HLA-DRB1'],
    'Exhuation':['PDCD1','CXCL13','LAYN'],
    'Interferon response':['STAT1','IFIT1','ISG15','CCR1'],
    'IL-17 response':['SLC4A10','KLRB1','TMIGD2','RORA','RORC','ZBTB16','IL26','IL17A','IL23R']
    'NME1 related':['NME1','NME2','MND1','SPC24','MYB']
}

# CD4
marker_genes_dict = {
    'Activation':['CD40LG'],
    'Naive':['TCF7','LEF1','CCR7','SELL','MAL','CXCR5','ADSL','IL16','IL7R'],
    'Effector memory':['TNF','AREG','TIMP1','CREM','CCL5','CAPG','GZMK','KLRG1','CX3CR1','TBX21'],
    'IL-17 response':['RORA','RORC','CCR6','IL23R','IL22','IL17A','IL17F','IL26'],
    'Follicular helper':['TOX','TOX2','IL21','GNG4','CD200','BCL6','ZBED2','CCL3','CCL4','IFNG','GZMB','LAG3','HAVCR2'],
    'Regulatory':['FOXP3','RTKN2','IL2RA','S1PR1','TNFRSF9','CTLA4','LAYN'],
    'Interferon response':['STAT1','IFIT1','ISG15','CCR1'],
    'NME1 related':['NME1','NME2','MND1','SPC24','MYB']
}

# NK
marker_genes_dict = {
    'Pro_inflammatory_markers':['XCL1', 'XCL2'],
    'Cytotoxic_markers':['FCGR3A', 'KLRK1', 'NCR3', 'NCR1', 'PRF1', 'GZMB'],
    'ILC3_like_markers':['KIT', 'IL7R', 'KLRB1']
}

# Myeloids
marker_genes_dict = {
    'Mast':['KIT'],
    'pDC_LILRA4':['GZMB','IL3RA'],
    'cDC1':['BATF3','IDO1','CADM1'],
    'cDC2_CD1C':['CD1C','FCER1A','CLEC10A'],# 'HLA_DQA1'
    'cDC3_LAMP3':['LAMP3','CCR7','FSCN1'],
    'Monocytes':['FCN1'],
    'Mono_CD14': ['CD14','S100A9'],
    'Mono_CD16':['FCGR3A'],
    'Macro_intermediate':['CD163','CD68'],
    'Macrophages':['FCGR2A','CSF1R'],
    'Macro_LYVE1':['PLTP'],
    'Macro_NLRP3':['NLRP3','EREG','IL1B','CXCL2'],
    'Macro_PPARG':['PPARG','MRC1','MSR1'],
    'Macro_ISG15':['ISG15','CXCL10','GBP1'],# 'IFITM3',
    'Macro_C1QC':['C1QA','C1QB','C1QC','APOE'],
    'Macro_SPP1':['SPP1','VEGFA','GPNMB']
}


# sequential numbers
def seq(start, end):
    start_num = int(start)
    end_num = int(end)
    return [str(num) for num in range(start_num, end_num + 1)]
seq('20','34') 

# maximum number of columns
pd.options.display.max_columns = None
pd.options.display.max_rows = None

# retain genes for dimentioanlity reduction
genes_to_retain = ["CD3D","CD3E","CD4","CD8A","GNLY","KLRF1","NCAM1",
                   "CD79A","MZB1","MS4A1","JCHAIN","TPSAB1",
                   "LYZ","PTPRC","CD163","CD68","TYROBP",
                   "COL1A1","DCN","FAP","PLVAP","CLDN5","FTL1","VWF","MKI67"]
for gene in genes_to_retain:
    adata.var.loc[gene, "highly_variable"] = True
adata = adata[:, adata.var["highly_variable"] == True]
adata.layers['counts'] = adata.obsm['counts'][:, adata.var["highly_variable"] == True].copy()

adata.obs['celltype_major'] = 'Unassigned'
adata.obs.loc[adata.obs["celltype_sc"].str.contains("Mal"), "celltype_major"] = "Melanoma"
adata.obs.loc[adata.obs["celltype_bped_main"].isin(["CD4+ T-cells","CD8+ T-cells"]), "celltype_major"] = "T_cells"
adata.obs.loc[adata.obs["celltype_sc"].isin(['NK']) & adata.obs["celltype_bped_main"].isin(["NK cells"]), "celltype_major"] = "NK_cells"
adata.obs.loc[adata.obs["celltype_sc"].isin(['Macrophage']) , "celltype_major"] = "Myeloids"
adata.obs.loc[adata.obs["celltype_sc"].str.contains("Endo."), "celltype_major"] = "Endothelium"
adata.obs.loc[adata.obs["celltype_bped_main"].str.contains("B-cells") & adata.obs["leiden_0.5"].isin(['11']), "celltype_major"] = "B_cells"
adata.obs.loc[adata.obs["leiden_0.5"].isin(['12']), "celltype_major"] = "Epithelium"
adata.obs.loc[adata.obs["leiden_0.5"].isin(['10','25']), "celltype_major"] = "Fibroblasts"
adata.obs.loc[adata.obs["leiden_0.5"].isin(['17']), "celltype_major"] = "Plasma"
adata.obs.loc[adata.obs["leiden_0.5"].isin(['22']), "celltype_major"] = "pDC"

sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes= 2000, layer="counts", batch_key="dataset", subset=True)
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="dataset")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="celltype_r2",
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)
sc.tl.leiden(adata)


#example script
#to run:
#mkdir <out_dir>
#python name_of_this_script.py <adata_path> <out_dir>
#gzip <out_dir>/*
import scanpy as sc
from scipy import io
import sys

adata = sc.read_h5ad(sys.argv[1])
out_dir = sys.argv[2]

adata = adata.raw.to_adata() #only if adata has RAW saved and thats what you want!!

with open(out_dir + '/barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')
        
with open(out_dir + '/features.tsv', 'w') as f:
    for item in ['\t'.join([x,x,'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')
        
io.mmwrite(out_dir +'/matrix', adata.X.T)

adata.obs.to_csv(sys.argv[1] + '.metadata.csv')

# cosg
n_gene=30
groupby='type'
cosg.cosg(adata,
    key_added='cosg',
    # use_raw=False, layer='log1p', ## e.g., if you want to use the log1p layer in adata
    mu=100,
    expressed_pct=0.1,
    remove_lowly_expressed=True,
    n_genes_user=5000,
    groupby=groupby)
df_tmp=pd.DataFrame(adata.uns['cosg']['names'][:10,]).T
marker_genes_list={idx: list(row.values) for idx, row in df_tmp.iterrows()}
marker_genes_list = {k: v for k, v in marker_genes_list.items() if not any(isinstance(x, float) for x in v)}
sc.pl.dotplot(adata, marker_genes_list,
             groupby=groupby,              
              swap_axes=False,
             standard_scale='var',
             cmap='Spectral_r')
[genename[0] for genename in adata.uns['cosg']['names']]
[score[0] for score in adata.uns['cosg']['scores']]
rnk=pd.DataFrame({'0':[score[0] for score in adata.uns['cosg']['names'][:200]], '1':[score[0] for score in adata.uns['cosg']['scores'][:200]]})
rnk.set_index('0', inplace=True)
msig.list_category(dbver="2023.2.Hs")
gp.get_library_name(organism='human')
pre_res = gp.prerank(rnk=rnk, # or rnk = rnk,
                     gene_sets='WikiPathway_2023_Human',
                     threads=4,
                     min_size=5,
                     max_size=1000,
                     permutation_num=1000, # reduce number to speed up testing
                     outdir=None, # don't write to disk
                     seed=6,
                     verbose=True, # see what's going on behind the scenes
                    )
term = pre_res.res2d.Term
# gp.gseaplot(res.ranking, term=term[i], **res.results[term[i]])
axs = pre_res.plot(terms=term[:10])
# filtering results
def convert_to_numeric(value):
    numerator, denominator = map(int, value.split('/'))
    return numerator / denominator
enr_res[enr_res['Overlap'].apply(convert_to_numeric)>0.5]