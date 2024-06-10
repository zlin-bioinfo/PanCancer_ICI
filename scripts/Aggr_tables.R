rm(list=ls())
pkgs <- c('Seurat', 'qs', 'dplyr', 'tidyr',  'janitor')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1) 

# Load datasets 
datasets <- c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC_Yost', 'CRC_Li', 'PCa_Hawley', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_Franken', 'HNSC_IMCISION', 'HNSC_Luoma', 'SCC_Yost', 'NSCLC_Liu')
list_metadata <- lapply(datasets, function(dataset){ 
  seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs'))
  metadata <- seu@meta.data
  return(metadata)})
names(list_metadata) <- datasets
qsave(list_metadata, '/bigdata/zlin/Melanoma_meta/tables/meta_list.qs') # save as lists
list_metadata <- qread('/bigdata/zlin/Melanoma_meta/tables/meta_list.qs')
# Aggregate all metadata
# Frequency by lineage
t_nk <- c('CD4_Naive','CD4_Tm_CREM-','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM', 'CD4_Tm_CCL5', 
          'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
          'CD4_Treg_Early', 'CD4_Treg_ISG15', 'CD4_Treg_TNFRSF9', 
          'CD4_Th_ISG15', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif',
          'CD8_Prolif', 'CD8_Naive', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
          'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13',
          'CD8_Tex_ISG15', 'CD8_Temra_CX3CR1', 'CD8_NK-like', 
          'MAIT', 'gdT', 'NK_CD56loCD16hi', 'NK_CD56hiCD16lo')
# bplasma <- c('Naive B', 'Memory IgM+ B', 'Memory IgM- B', 'GC-like B', 'Plasmablast', 'Plasma')
bplasma <- c('B_Naive', 'B_ISG15', 'B_HSP', 'B_MT2A', 'ACB_EGR1', 'ACB_NR4A2', 'ACB_CCR7', 'B_Memory', 'B_AtM',
             'GCB_Pre', 'GCB_SUGCT', 'GCB_LMO2', 'GCB_Prolif', 'Plasmablast', 'Plasma_cell')
mye <- c('Mast','pDC','cDC1', 
         'cDC2_CD1C', 'cDC2_IL1B','cDC2_ISG15', 'cDC2_CXCL9', 'DC_LC-like', 'MigrDC', 'MoDC', 
         'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
         'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro_ISG15', 
         'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2')
nonimmune <- c('EC_lymphatic','EC_vascular','EndMT','CAF_inflammatory', 'CAF_adipogenic', 'CAF_PN', 'CAF_AP', 'Myofibroblast')
list_metadata <- list_metadata[c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC_Yost', 'CRC_Li', 'PCa_Hawley', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_Franken', 'HNSC_IMCISION', 'HNSC_Luoma', 'SCC_Yost', 'NSCLC_Liu')]
common_cols <- Reduce(intersect, lapply(list_metadata, colnames))
meta_list <- lapply(list_metadata, function(metadata){
  print(unique(metadata$dataset))
  metadata <- metadata[, common_cols]
  metadata$component <- 'celltype'
  metadata$component[metadata$celltype_r2 %in% t_nk] <- 'T_NK'
  metadata$component[metadata$celltype_r2 %in% bplasma] <- 'Bplasma'
  metadata$component[metadata$celltype_r2 %in% mye] <- 'Myeloids'
  metadata$component[metadata$celltype_r2 %in% nonimmune] <- 'Non-immune'
  metadata <- metadata |> 
    select(dataset, patient, sample, time_point, celltype_main, celltype_r2, interval, cancertype, response, res_metric, treatment, component, Phase) |> 
    group_by(sample) |> 
    dplyr::mutate(count_sample = n()) |> 
    group_by(celltype_main, .add = TRUE) |> 
    dplyr::mutate(count_main = n()) |> 
    dplyr::mutate(freq_main = count_main/count_sample) |> 
    group_by(celltype_r2, sample) |>   
    dplyr::mutate(count_r2 = n()) |>  
    dplyr::mutate(freq_r2 = count_r2/count_sample) |> 
    group_by(component, sample) |>  
    dplyr::mutate(count_component = n()) |> 
    dplyr::mutate(freq_r2_comp = count_r2/count_component) |> 
    ungroup()
})
meta_combi <- do.call(rbind, meta_list) |> sapply(as.character) |> as.data.frame() 
meta_combi$prior <- 'No'
meta_combi$prior[meta_combi$dataset == 'BRCA_Bassez2'] <- 'Yes'
meta_combi$prior[meta_combi$cancertype == 'BCC' & !meta_combi$patient %in% c('BCC_Yost_su004')] <- 'Yes'
meta_combi$prior[meta_combi$cancertype == 'SCC'] <- 'Yes'
meta_combi$dataset <- as.character(meta_combi$dataset)
meta_combi$freq_r2 <- as.numeric(meta_combi$freq_r2)
meta_combi$modality <- ifelse(meta_combi$treatment %in% c('aPD1+CTLA4', 'aPDL1+CTLA4'), 'Dual', 'Mono')
meta_combi$count_r2 <- as.numeric(meta_combi$count_r2)
meta_combi$interval <- as.numeric(meta_combi$interval)
rm_pt <- meta_combi |> tabyl(patient, time_point) |> filter(Pre < 200 | Post <200) |> pull(patient)
# "CRC_Li_P28" "HNSC_Franken_P16" "HNSC_Franken_P20" "SKCM_Becker_P10"  "SKCM_Becker_P6"  
meta_int <- meta_combi |> filter(!patient %in% rm_pt)
meta_int$freq_r2_comp <- as.numeric(meta_int$freq_r2_comp)
meta_int$interval <- as.numeric(meta_int$interval)
meta_int$int_cat <- ifelse(meta_int$interval < 21, '< 21d', '>= 21d')
write.csv(meta_int, '/bigdata/zlin/Melanoma_meta/tables/meta_int.csv', row.names = F)
# Set the relative frequency for major groups to 0 for some samples
# meta_filtered <- meta_int |> mutate(freq_r2_comp = case_when((count_component < 40 & count_r2 < 3) ~ 0, .default = freq_r2_comp))
# write.csv(meta_filtered, '/bigdata/zlin/Melanoma_meta/tables/meta_filtered.csv', row.names = F)

# Metadata by sample
meta_sample <- meta_int |> distinct(sample, .keep_all = T)
write.csv(meta_sample, '/bigdata/zlin/Melanoma_meta/tables/meta_sample.csv', row.names = F)

# Metadata by patient
meta_patient <- meta_int |> distinct(patient, .keep_all = T)
meta_patient$dataset[meta_patient$dataset %in% c('BCC_Yost','SCC_Yost')] <- 'BCC/SCC_Yost'
write.csv(meta_patient, '/bigdata/zlin/Melanoma_meta/tables/meta_patient.csv', row.names = F)
# Ro/e
roie <- function(meta_combi, by_component){
  roie_list <- lapply(unique(meta_combi$patient), function(p){
    meta_combi <- subset(meta_combi, patient == p)
    if (by_component == T){
      roie_list_c <- list()
      for (i in unique(meta_combi$component)){
        meta_sub <- meta_combi |> filter(component == i)
        observed_table <- table(meta_sub$time_point, meta_sub$celltype_r2) 
        expected_table <- chisq.test(observed_table)$expected
        roie <- data.frame(observed_table/expected_table)
        observed_df <- data.frame(observed_table)
        roie$count <- observed_df$Freq
        roie$component <- i
        roie_list_c[[i]] <- roie
      }
      roie <- do.call(rbind, roie_list_c)
    } else {
      observed_table <- table(meta_combi$time_point, meta_combi$celltype_r2) 
      expected_table <- chisq.test(observed_table)$expected
      roie <- data.frame(observed_table/expected_table)
      observed_df <- data.frame(observed_table)
      roie$count <- observed_df$Freq
    }
    names(roie)[1] <- 'time'
    roie$patient <- p
    roie$time <- factor(roie$time, levels = c('Pre','Post'))
    names(roie)[2:3] <- c('celltype','ratio')
    return(roie)
  })
  roie_mat <- do.call(rbind, roie_list)
  meta <- meta_combi |> 
    distinct(patient, .keep_all = T) |> 
    select(interval, cancertype, response, res_metric, treatment, patient, time_point, prior, dataset, modality)
  roie_mat <- left_join(roie_mat, meta, by = 'patient')
  roie_mat <- select(roie_mat, !time_point)
  return(roie_mat)
}
mat_roie <- roie(meta_int, by_component = T)
write.csv(mat_roie, '/bigdata/zlin/Melanoma_meta/tables/roie.csv', row.names = F)




