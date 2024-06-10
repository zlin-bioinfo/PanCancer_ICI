rm(list=ls())
pkgs <- c('qs','dplyr','ggplot2','SingleCellExperiment','Seurat','glmGamPoi')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

datasets <- list('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 'BCC_Yost', 'SCC_Yost', 'HNSC_IMCISION', 'HNSC_Luoma', 'NSCLC_Liu', 'CRC_Li', 'PCa_Hawley')
list_seu <- lapply(datasets, function(dataset){
  seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) 
  return(seu)
})

sce <- merge(x=list_seu[[1]], y=list_seu[2:length(datasets)], add.cell.ids = datasets) |> JoinLayers()|> as.SingleCellExperiment()
sce <- sce[rowSums(counts(sce)) > 30,]

library(scater)
qc <- perCellQCMetrics(sce)

sample_id = "sample"
group_id = "time_point"
celltype_id = "celltype_main"
covariates = "dataset"