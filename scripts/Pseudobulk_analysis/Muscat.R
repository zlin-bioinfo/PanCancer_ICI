rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','ggsci','dior', 'qs','BiocParallel','muscat','scater','muscatWrapper')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
options(future.globals.maxSize = 1e9)

datasets <- c('SKCM_Becker', 'BRCA_Bassez2')
list_seu <- lapply(datasets, function(dataset){
  seu <- read_h5(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/R2_annotated.h5')) %>% 
    subset(subset = celltype_main == 'gdT', invert = T)
  return(seu)})
sce <- merge(list_seu[[1]], y = list_seu[[2]]) %>% as.SingleCellExperiment()
sce <- sce[rowSums(counts(sce) > 0) >= 20, ]
sce <- sce[, sce$response != 'NE']
sce$response_timepoint <- paste0(sce$response, '_', sce$time_point)
sce$interval <- as.numeric(sce$interval)
sce$prior <- 'No'
sce$prior[sce$dataset == 'BRCA_Bassez2'] <- 'Yes'
sce$prior[sce$dataset == 'BCC_Yost' & !sce$patient %in% c('BCC_Yost_su004')] <- 'Yes'
sce$prior[sce$dataset == 'SCC_Yost'] <- 'Yes'
sce$dataset <- as.character(sce$dataset)
sce$dataset[sce$dataset %in% c('BCC_Yost', 'SCC_Yost')] <- 'BCC/SCC_Yost'
sce$fraction_r2 <- as.numeric(sce$fraction_r2)
sce$modality <- ifelse(sce$treatment == 'aPD1+CTLA4', 'Dual', 'Mono')
sce$count_r2 <- as.numeric(sce$count_r2)
rm(seu, list_seu)

sample_id = "sample"
group_id = "response_timepoint"
celltype_id = "celltype_main"
# covariates = c("interval", "modality", "prior")
# covariates = "modality"
covariates =NA
batches = "dataset"
min_cells = 10
abundance_output = get_abundance_info(sce, sample_id, group_id, celltype_id, min_cells, covariates = covariates)
head(abundance_output$abundance_data)

contrasts_oi = c("'(OnE-PreE)-(OnNE-PreNE)','(OnNE-PreNE)-(OnE-PreE)'")
contrast_tbl = tibble(contrast =
                        c("(OnE-PreE)-(OnNE-PreNE)", "(OnNE-PreNE)-(OnE-PreE)"),
                      group = c("OnE","OnNE")) 

muscat_output = muscat_analysis(
  sce = sce,
  celltype_id = celltype_id,
  sample_id = sample_id,
  group_id = group_id,
  covariates = covariates,
  contrasts_oi = contrasts_oi,
  contrast_tbl = contrast_tbl)






















sce <- prepSCE(sce, 
               kid = "celltype_main", # subpopulation assignments
               gid = "response",  # group IDs (ctrl/stim)
               sid = "patient",   # sample IDs (ctrl/stim.1234)
               drop = TRUE); sce

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

t(table(sce$cluster_id, sce$sample_id))
sce <- runUMAP(sce, pca = 20)

.plot_dr <- function(sce, dr, col)
  plotReducedDim(sce, dimred = dr, colour_by = col) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_minimal() + theme(aspect.ratio = 1)
# downsample to max. 100 cells per cluster
cs_by_k <- split(colnames(sce), sce$cluster_id)
cs100 <- unlist(sapply(cs_by_k, function(u) 
  sample(u, min(length(u), 100))))

# plot t-SNE & UMAP colored by cluster & group ID
for (col in c("cluster_id", "group_id"))
  .plot_dr(sce[, cs100], 'UMAP', 'cluster_id')

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
# one sheet per subpopulation
assayNames(pb)

t(head(assay(pb)))
(pb_mds <- pbMDS(pb))

# run DS analysis
res <- pbDS(pb, verbose = FALSE)
# access results table for 1st comparison
tbl <- res$table[[1]]
# one data.frame per cluster
names(tbl)

# view results for 1st cluster
k1 <- tbl[[1]]
head(format(k1[, -ncol(k1)], digits = 2))

# construct design & contrast matrix
ei <- metadata(sce)$experiment_info
mm <- model.matrix(~ 0 + ei$group_id)
dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
contrast <- makeContrasts("RE-NR", levels = mm)
# run DS analysis
pbDS(pb, design = mm, contrast = contrast)

pb




