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
rm(list_seu)
reduced_sce <- pseudobulk(sce, group_by = vars(patient, condition = time_point, celltype_main), n_cells = n())
rm(sce)

celltype_main <- unique(reduced_sce$celltype_main)
lapply(celltype_main, function(celltype){
  print(celltype)
  fit <- glm_gp(reduced_sce[, reduced_sce$celltype_main == celltype], 
                design = ~ condition + patient, size_factor = "ratio", verbose = TRUE)
  qsave(fit, paste0('/bigdata/zlin/Melanoma_meta/data/glmGamPoi_fit_', celltype, '.qs'))
  })

res_list <- lapply(celltype_main, function(celltype){
  print(celltype)
  fit <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/glmGamPoi_fit_', celltype, '.qs'))
  res <- test_de(fit, contrast = cond(condition = "Post") - cond(condition = "Pre"))
  return(res)
  })

qsave(res_list, '/bigdata/zlin/Melanoma_meta/data/glmGamPoi_list_DE.qs')

fit <- qread('/bigdata/zlin/Melanoma_meta/data/glmGamPoi_fit_Macro.qs')
res <- test_de(fit, contrast = cond( condition = "Post") - cond(condition = "Pre"))

ggplot(res, aes(x = lfc, y = - log10(pval))) + 
  geom_point(aes(color = adj_pval < 0.05 & (abs(lfc) > 1)), size = 0.5)

res |> filter(adj_pval<0.05 & (abs(lfc)>1)) |> 
  arrange(desc(lfc))




