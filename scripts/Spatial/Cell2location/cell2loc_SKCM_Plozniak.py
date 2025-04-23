# CUDA_VISIBLE_DEVICES=1 python /bigdata/zlin/PanCancer_ICI/scripts/Spatial/cell2loc_SKCM_Plozniak.py
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import squidpy as sq
import pandas as pd
import cell2location
import anndata as ad
from matplotlib import rcParams
results_folder = '/bigdata/zlin/PanCancer_ICI/data/SKCM_Plozniak/cell2loc_results'
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
# load reference data
dir_datset='/bigdata/zlin/PanCancer_ICI/data/SKCM_Plozniak'
adata_ref=sc.read_mtx(filename= dir_datset + '/seu_final/matrix.mtx')
adata_bc=pd.read_csv(filepath_or_buffer= dir_datset + '/seu_final/barcodes.tsv', header=None)
adata_features=pd.read_csv(filepath_or_buffer= dir_datset + '/seu_final/features.tsv',header=None)
adata_obs=pd.read_csv(filepath_or_buffer= dir_datset + '/seu_final/metadata.tsv', sep='\t')
adata_ref= adata_ref.T
adata_ref.obs=adata_obs
adata_ref.obs.index=adata_bc[0].tolist()
adata_ref.var['gene_name']= adata_features[0].tolist()
adata_ref.var.index= adata_ref.var['gene_name']
adata_ref=adata_ref[~adata_ref.obs['celltype_main'].isin(['Cycling T/NK','Cycling','Melanocytes(CNA-)'])]
# QC
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# filter the object
adata_ref = adata_ref[:, selected].copy()
# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='sample',
                        # cell type, covariate used for constructing signatures
                        labels_key='celltype_r2',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['patient']
                       )
# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
# view anndata_setup as a sanity check
mod.view_anndata_setup()
mod.train(max_epochs=200)
# export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500}
)
# Save model
mod.save(f"{ref_run_name}", overwrite=True)
# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

# # reload model and output h5ad:
# adata_file = f"{ref_run_name}/sc.h5ad"
# adata_ref = sc.read_h5ad(adata_file)
# mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

# load visium data
samples = ['Pt10_K28284','Pt11_K28439','Melanoma2_3','Pt12_K29368','Pt09_K3000','Melanoma2_1']
dict_adata_vis = {
    sample: sq.read.visium(path=f'/bigdata/zlin/PanCancer_ICI/data/SKCM_Plozniak/{sample}/outs', library_id=sample)
    for sample in samples
}
# Iterate through dictionary and process each sample
for sample, adata_vis in dict_adata_vis.items():
    print(sample)
    run_name = f'{results_folder}/cell2location_map/{sample}'
    adata_vis.var_names_make_unique()
    # QC
    from scipy.sparse import csr_matrix
    adata_vis.X = adata_vis.X.toarray()
    sc.pp.calculate_qc_metrics(adata_vis, inplace=True)
    adata_vis.X = csr_matrix(adata_vis.X)
    adata_vis.var['mt'] = [gene.startswith('mt-') for gene in adata_vis.var_names]
    adata_vis.obs['mt_frac'] = adata_vis[:, adata_vis.var['mt'].tolist()].X.sum(1).A.squeeze()/adata_vis.obs['total_counts']
    # add sample name to obs names
    adata_vis.obs["sample"] = sample
    adata_vis.obs_names = adata_vis.obs["sample"] + '_' + adata_vis.obs_names
    adata_vis.obs.index.name = 'spot_id'
    # mitochondria-encoded (MT) genes should be removed for spatial mapping
    adata_vis.obsm['mt'] = adata_vis[:, adata_vis.var['mt'].values].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var['mt'].values]
    sc.pp.calculate_qc_metrics(adata_vis, inplace=True)
    sc.pp.filter_cells(adata_vis, min_counts=200)
    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()
    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis)
    # create and train the model
    mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=4,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20)
    mod.view_anndata_setup()
    mod.train(max_epochs=15000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1)
    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
    )
    # Save model
    mod.save(f"{run_name}", overwrite=True)
    # Save anndata object with results
    adata_file = f"{run_name}/sp.h5ad"
    adata_vis.write(adata_file)
    adata_file


