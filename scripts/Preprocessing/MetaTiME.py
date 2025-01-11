import diopy
import numpy as np
import pandas as pd
import scanpy as sc 
import anndata as ad
import scanpy.external as sce
import matplotlib.pyplot as plt
# import metatime
# from metatime import config
# from metatime import mecmapper
# from metatime import mecs
# from metatime import annotator
import pickle

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

becker = diopy.input.read_h5(file='/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/processed_major_after.h5')
franken = diopy.input.read_h5(file='/bigdata/zlin/Melanoma_meta/data/HNSC_Franken/processed_major_after.h5')
bassez1 = diopy.input.read_h5(file='/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/filtered.h5')
bassez2 = diopy.input.read_h5(file='/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/filtered.h5')
yost = diopy.input.read_h5(file='/bigdata/zlin/Melanoma_meta/data/BCC_Yost/filtered.h5')
crc_li = diopy.input.read_h5(file='/bigdata/zlin/Melanoma_meta/data/CRC_Li/processed_major_after.h5')
tnbc_li = diopy.input.read_h5(file = '/bigdata/zlin/Melanoma_meta/data/TNBC_Li/processed_major_after.h5')

def counts_reverse(adata):
    adata.X = adata.layers['counts'].copy()
    return adata

counts_reverse(bassez1)
counts_reverse(bassez2)
counts_reverse(yost)

becker = becker[~becker.obs['celltype_major'].isin(['Melanoma','Epithelium'])]
franken = franken[~franken.obs['celltype_major'].isin(['Epi/Malignant'])]
bassez1 = bassez1[~bassez1.obs['celltype_major'].isin(['Malignant'])]
bassez2 = bassez2[~bassez2.obs['celltype_major'].isin(['Malignant'])]
yost = yost[~yost.obs['celltype_major'].isin(['Malignant','Melanocytes'])]
crc_li = crc_li[~crc_li.obs['celltype_major'].isin(['Epi/Malignant'])]

# Load the pre-trained MeCs
mecmodel = mecs.MetatimeMecs.load_mec_precomputed()
# Load functional annotation for MetaTiME-TME
mectable = mecs.load_mecname(mecDIR = config.SCMECDIR, mode ='table' )
mecnamedict = mecs.getmecnamedict_ct(mectable) 

datasets = [yost[:,mecmodel.mec_score.index.intersection(yost.var_names)], 
            crc_li[:,mecmodel.mec_score.index.intersection(crc_li.var_names)], 
            tnbc_li[:,mecmodel.mec_score.index.intersection(tnbc_li.var_names)], 
            becker[:,mecmodel.mec_score.index.intersection(becker.var_names)], 
            bassez1[:,mecmodel.mec_score.index.intersection(bassez1.var_names)], 
            bassez2[:,mecmodel.mec_score.index.intersection(bassez2.var_names)], 
            franken[:,mecmodel.mec_score.index.intersection(franken.var_names)]]

print('Get started!')
pred_list = list()
for adata in datasets:
    adata.layers['counts'] = adata.X.copy()
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset=False, layer='counts', flavor='seurat_v3')
    sc.pp.regress_out(adata, ['total_counts'])
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, svd_solver='arpack')
    sce.pp.bbknn(adata, batch_key='sample')
    sc.tl.umap(adata)
    adata = annotator.overcluster(adata)
    pdata = mecmapper.projectMecAnn(adata, mecmodel.mec_score)
    projmat, mecscores = annotator.pdataToTable(pdata, mectable, gcol = 'overcluster')
    projmat, gpred, gpreddict = annotator.annotator(projmat, mecnamedict, gcol = 'overcluster')
    pred_list.append(pd.concat([projmat['MetaTiME_overcluster'].str.split(': ').str.get(1).str.split('_').str.get(1), mecscores], axis=1))

# Save the object to a file
with open('/bigdata/zlin/Melanoma_meta/data/pred_MetaTiME.pkl', 'wb') as file:
    pickle.dump(pred_list, file)

print('Done!')

with open('/bigdata/zlin/Melanoma_meta/data/pred_MetaTiME.pkl', 'rb') as file:
    pred_list = pickle.load(file)

datasets = [yost, crc_li, tnbc_li, becker, bassez1, bassez2, franken]
dataset_names = ['BCC_Yost', 'CRC_Li', 'TNBC_Li', 'SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'HNSC_Franken']
for i in range(len(datasets)):
    if datasets[i].obs.shape[0] == pred_list[i].shape[0]:
        datasets[i].obs = pd.concat([datasets[i].obs, pred_list[i]], axis=1)
        output_file = '/bigdata/zlin/Melanoma_meta/data/' + dataset_names[i] + '/non_malignant.h5'
        diopy.output.write_h5(datasets[i], file = output_file)
        print(dataset_names[i] + ' Done!')
    else: 
        break

