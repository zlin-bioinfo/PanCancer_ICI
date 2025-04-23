import scanpy as sc
import decoupler as dc
import plotnine as p9
import liana as li
import os
from liana.method import MistyData, genericMistyData, lrMistyData
from liana.method.sp import RandomForestModel, LinearModel, RobustLinearModel
# Start overall timing
import time
start_time = time.time()
# Create directory 
directory_name = "/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/MISTy/"
if not os.path.exists(directory_name):
    os.makedirs(directory_name)
# Load and normalize data
samples = ['HCC7NR', 'HCC6NR', 'HCC4R', 'HCC3R', 'HCC2R', 'HCC5NR', 'HCC1R']
dict_adata = {
    sample: sc.read(filename=f'/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/cell2loc_results/cell2location_map/{sample}/sp.h5ad', library_id=sample)
    for sample in samples
}
# Iterate through dictionary and process each sample
for sample, adata in dict_adata.items():
    sample_start_time = time.time()
    print(sample)
    sc.pp.filter_cells(adata, min_counts=300)
    sc.pp.filter_cells(adata, min_genes=200)
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Use 5% quantile of the posterior distribution, representing the value of cell abundance that the model has high confidence
    adata.obsm['q05_cell_abundance_w_sf'].columns=adata.uns['mod']['factor_names'].copy()
    comps = li.ut.obsm_to_adata(adata, 'q05_cell_abundance_w_sf')
    # Colocalization analysis(intra view only)
    directory_name = "/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/MISTy/colocalization/"
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    misty = genericMistyData(intra=comps)
    misty(model=RandomForestModel, n_jobs=-1, verbose = True)
    misty.uns['target_metrics'].to_csv(f'/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/MISTy/colocalization/{sample}_target_metrics_RF.csv',index=False)
    misty.uns['interactions'].to_csv(f'/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/MISTy/colocalization/{sample}_interactions_RF.csv',index=False)
    misty(model=LinearModel, k_cv=10, seed=1000,verbose = True)
    misty.uns['target_metrics'].to_csv(f'/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/MISTy/colocalization/{sample}_target_metrics_LM.csv',index=False)
    misty.uns['interactions'].to_csv(f'/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/MISTy/colocalization/{sample}_interactions_LM.csv',index=False)
    # Custom view(PROGENy and TFs)
    directory_name = "/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/MISTy/PROGENY_TFs/"
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    # PROGENy activity estimation
    progeny = dc.get_progeny(organism='human', top=500)
    # use multivariate linear model to estimate activity
    dc.run_mlm(
        mat=adata,
        net=progeny,
        source='source',
        target='target',
        weight='weight',
        verbose=True,
        use_raw=False)
    # extract progeny activities as an AnnData object
    acts_progeny = li.ut.obsm_to_adata(adata, 'mlm_estimate')   
    # get TFs prior knowledge
    net = dc.get_collectri()
    # Estimate activities
    dc.run_ulm(
        mat=adata,
        net=net,
        verbose=True,
        use_raw=False)
    # extract activities
    acts_tfs = li.ut.obsm_to_adata(adata, 'ulm_estimate')
    # Calculate spatial neighbors
    li.ut.spatial_neighbors(acts_tfs, cutoff=0.1, bandwidth=200, set_diag=False)
    # transfer spatial information to progeny activities
    # NOTE: spatial connectivities can differ between views, but in this case we will use the same
    acts_progeny.obsm['spatial'] = acts_tfs.obsm['spatial']
    acts_progeny.obsp['spatial_connectivities'] = acts_tfs.obsp['spatial_connectivities']
    misty = MistyData(data={"intra": comps, "TFs": acts_tfs, "Pathways": acts_progeny})
    misty
    misty(model=LinearModel, verbose=True, bypass_intra=True)
    misty.uns['target_metrics'].to_csv(f'/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/MISTy/PROGENY_TFs/{sample}_target_metrics.csv',index=False)
    misty.uns['interactions'].to_csv(f'/bigdata/zlin/PanCancer_ICI/data/GSE238264_RAW/MISTy/PROGENY_TFs/{sample}_interactions.csv',index=False)
    sample_end_time = time.time()
    print(f"Processing time for {sample}: {(sample_end_time - sample_start_time) / 60:.2f} minutes")
# End overall timing
end_time = time.time()
print(f"Total execution time: {(end_time - start_time)/3600 :.2f} hours")