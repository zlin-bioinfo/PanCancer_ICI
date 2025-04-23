pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','COSG','BPCells','SeuratExtend','MetBrewer','ggplot2','slingshot','CytoTRACE2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
setwd("/home/zlin/workspace/PanCancer_ICI")

# T NK
lapply(cohorts, function(cohort){
  print(cohort)
  seu <- qs_read(paste0('data/', cohort, '/seu_r2.qs2'))
  seu <- subset(seu, subset = celltype_main == "CD8+T")
  cytotrace2_result <- cytotrace2(seu, ncores = 16, is_seurat = T, slot_type = "counts", species = "human", seed = 123)
  print('done')
  write.csv(cytotrace2_result@meta.data, paste0('data/', cohort, '/cytotrace2_cd8.csv'))
})
cohorts <- c('SKCM_Becker', 'SKCM_Plozniak', 'BCC_Yost',
             'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 
             'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 'NSCLC_Yan', 
             'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo','RCC_Bi')
lapply(cohorts, function(cohort){
  print(cohort)
  seu <- qs_read(paste0('data/', cohort, '/seu_r2.qs2'))
  seu <- subset(seu, subset = celltype_main %in% c("B", "Plasma"))
  cytotrace2_result <- cytotrace2(seu, ncores = 16, is_seurat = T, slot_type = "counts", species = "human", seed = 123)
  print('done')
  write.csv(cytotrace2_result@meta.data, paste0('data/', cohort, '/cytotrace2_bplasma.csv'))
})



