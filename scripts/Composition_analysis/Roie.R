rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','gtools','ComplexHeatmap','dittoSeq','RColorBrewer','ggpubr','tibble','epitools','dior','cowplot','ggnewscale','export','forcats','qs','purrr','effsize','ggrepel','lme4','lmerTest','rstatix','viridis','janitor')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Patient distribution
# add lineage component
t_nk <- c('CD4_Naive','CD4_Tm_TNF','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM',
          'CD4_Tm_CCL5', 'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
          'CD4_Treg_Early', 'CD4_Treg_ISG', 'CD4_Treg_TNFRSF9', 'CD4_Th_ISG', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif',
          'CD8_Prolif', 'CD8_Naive', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
          'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13', 'CD8_Tex_OXPHOS-', 
          'CD8_ISG', 'CD8_Temra_CX3CR1', 'CD8_NK-like', 'MAIT', 'gdT', 'NK_CD56loCD16hi', 'NK_CD56hiCD16lo')
bplasma <- c('Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Exhausted B cells', 'Plasma')
mye <- c('Mast','pDC','cDC1', 'cDC2_CD1C', 'cDC2_IL1B', 'cDC2_ISG15', 'cDC2_CXCL9', 'DC_LC-like', 'cCD2_MoDC', 'DC_Migr',
         'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
         'Macro_NLRP3','Macro_ISG15', 'Macro_INHBA', 'Macro_FN1', 'Macro_SPP1',
         'Macro_LYVE1','Macro_IL1B','Macro_C1QC','Macro_TREM2')
nonimmune <- c('EC_lymphatic','EC_vascular','EndMT','CAF_inflammatory', 'CAF_adipogenic', 'CAF_PN', 'CAF_AP', 'Myofibroblast')
list_metadata <- qread('/bigdata/zlin/Melanoma_meta/tables/meta_list.qs')
list_metadata <- list_metadata[c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC_Yost', 'CRC_Li', 'PCa_Hawley', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_IMCISION', 'HNSC_Luoma', 'SCC_Yost', 'NSCLC_Liu')]
meta_list <- lapply(list_metadata, function(metadata){
  print(unique(metadata$dataset))
  metadata$component <- 'celltype'
  metadata$component[metadata$celltype_r2 %in% t_nk] <- 'T_NK'
  metadata$component[metadata$celltype_r2 %in% bplasma] <- 'Bplasma'
  metadata$component[metadata$celltype_r2 %in% mye] <- 'Myeloids'
  metadata$component[metadata$celltype_r2 %in% nonimmune] <- 'Non-immune'
  metadata <- metadata |> 
    select(dataset, patient, sample, time_point, celltype_main, celltype_r2, interval, cancertype, response, res_metric, treatment, component, Phase) |> 
    group_by(sample) |> 
    mutate(count_sample = n()) |> 
    group_by(celltype_main, .add = TRUE) |> 
    mutate(count_main = n()) |> 
    mutate(freq_main = count_main/count_sample) |> 
    group_by(celltype_r2, sample) |>  
    mutate(count_r2 = n()) |>  
    mutate(freq_r2 = count_r2/count_sample) |> 
    group_by(component, sample) |>  
    mutate(count_component = n()) |> 
    mutate(freq_r2_comp = count_r2/count_component) 
})
meta_combi <- do.call(rbind, meta_list) |> sapply(as.character) |> as.data.frame() 
meta_combi$prior <- 'No'
meta_combi$prior[meta_combi$dataset == 'BRCA_Bassez2'] <- 'Yes'
meta_combi$prior[meta_combi$cancertype == 'BCC' & !meta_combi$patient %in% c('BCC/SCC_Yost_su004')] <- 'Yes'
meta_combi$prior[meta_combi$cancertype == 'SCC'] <- 'Yes'
meta_combi$dataset <- as.character(meta_combi$dataset)
meta_combi$freq_r2 <- as.numeric(meta_combi$freq_r2)
meta_combi$modality <- ifelse(meta_combi$treatment == 'aPD1+CTLA4', 'Dual', 'Mono')
meta_combi$count_r2 <- as.numeric(meta_combi$count_r2)
meta_combi$interval <- as.numeric(meta_combi$interval)
rm_pt <- meta_combi |> tabyl(patient, time_point) |> filter(Pre < 200 | Post <200) |> pull(patient)
meta_combi <- meta_combi |> filter(!patient %in% rm_pt)
write.csv(meta_combi, '/bigdata/zlin/Melanoma_meta/tables/meta_int.csv')

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
roie_mat <- roie(meta_combi, by_component = T)
write.csv(roie_mat, '/bigdata/zlin/Melanoma_meta/tables/roie_mat.csv')

# Load data
meta_combi <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_int.csv')
df_anno <- distinct(meta_combi, celltype_r2, celltype_main, .keep_all = T) |> select(celltype_r2, celltype_main)
roie_mat <- read.csv('/bigdata/zlin/Melanoma_meta/tables/roie_mat.csv') |> select(!X)

# Heatmap of freq
meta <- meta_combi |> 
  distinct(patient, .keep_all = T) 
timepoint <- 'Pre'
mat_pre <- meta_combi |> 
  filter(time_point == timepoint) |>
  distinct(patient, celltype_r2, .keep_all = T) |> 
  select(celltype_r2, patient, freq_r2_comp) |> 
  pivot_wider(names_from = celltype_r2, values_from = freq_r2_comp) |> 
  column_to_rownames(var = 'patient') |> 
  t() 
rowname_pre <- rownames(mat_pre)
mat_pre <- apply(mat_pre, 2, as.numeric)
rownames(mat_pre) <- rowname_pre
timepoint <- 'Post'
mat_post <- meta_combi |> 
  filter(time_point == timepoint) |>
  distinct(patient, celltype_r2, .keep_all = T) |> 
  select(celltype_r2, patient, freq_r2_comp) |> 
  pivot_wider(names_from = celltype_r2, values_from = freq_r2_comp) |> 
  column_to_rownames(var = 'patient') |> 
  t() 
# rowname_post <- rownames(mat_post)
# mat_post <- apply(mat_post, 2, as.numeric)
# rownames(mat_post) <- rowname_post
# column annotation
col_ha = HeatmapAnnotation(
  dataset = meta$dataset,
  cancertype = meta$cancertype,
  treatment = meta$treatment,
  prior = meta$prior,
  response = meta$response,
  res_metric = meta$res_metric,
  col = list(dataset = structure(names = unique(meta$dataset), dittoColors()[1:length(unique(meta$dataset))]),
             cancertype = structure(names = unique(meta$cancertype), pal_jco()(length(unique(meta$cancertype)))),
             treatment = structure(names = unique(meta$treatment), pal_nejm()(length(unique(meta$treatment)))),
             prior = structure(names = unique(meta$prior), pal_rickandmorty()(length(unique(meta$prior)))),
             response = structure(names = unique(meta$response), pal_startrek()(length(unique(meta$response)))),
             res_metric = structure(names = unique(meta$res_metric), pal_d3()(length(unique(meta$res_metric))))
  ),
  annotation_legend_param = list(
    dataset = list(title = "Dataset"),
    cancertype = list(title = "Cancer type"),
    treatment = list(title = "Treatment"),
    prior = list(title = "Prior-Tx"),
    response = list(title = "Response"),
    res_metric = list(title = "Response Metric")),
  show_annotation_name = FALSE
)
# row annotation
color_main <- c('#E31A1C', '#1F78B4', '#A6CEE3', '#CAB2D6', '#6A3D9A', '#FF7F00', '#B15928', '#FFFF99', '#B2DF8A', '#33A02C', '#FB9A99', '#FDBF6F')
row_ha = rowAnnotation(
  celltype_main = df_anno$celltype_main[match(rownames(mat_pre), df_anno$celltype_r2)],
  col = list(celltype_main = structure(names = c('CD4+T', 'CD8+T', 'NK', 'B', 'Plasma', 'pDC', 'Mast', 'cDC', 'Mono', 'Macro', 'Endo', 'CAF'), color_main)),
  annotation_legend_param = list(celltype_main = list(title = "Main")),
  show_annotation_name = FALSE
)
# celltype_order_main <- c('CD4+T', 'CD8+T', 'NK', 'B', 'Plasma', 'pDC','Mast','cDC', 'Mono', 'Macro','Endo','CAF')
celltype_order_r2 <- c('CD4_Naive','CD4_Tn_ADSL','CD4_Tm_TNF','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM',
                       'CD4_Tm_CCL5', 'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
                       'CD4_Treg_TNFRSF9-', 'CD4_Treg_S1PR1', 'CD4_Treg_TNFRSF9', 'CD4_Treg_ISG', 'CD4_Th_ISG', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_STR',
                       'CD8_Naive', 'CD8_MAIT_SLC4A10', 'CD8_Tm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
                       'CD8_Tex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13', 'CD8_Tex_OXPHOS-', 'CD8_STR', 'CD8_ISG', 'CD8_Temra_CX3CR1', 'CD8_NK-like',
                       'NK_CD56lowCD16hi','NK_CD56hiCD16low', 'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Exhausted B cells', 'Plasma', 
                       'Mast','pDC','cDC1', 'cDC2', 'cDC3', 'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
                       'Macro_NLRP3','Macro_IL1B','Macro_LYVE1','Macro_C1QC','Macro_TREM2', 'Macro_ISG15', 'Macro_INHBA', 'Macro_FN1', 'Macro_SPP1',
                       'EC_vascular', 'EC_lymphatic', 'EC_EndMT', 'Myofibroblast', 'CAF_inflammatory', 'CAF_adipogenic', 'CAF_EndMT', 'CAF_AP', 'CAF_PN')
pdf('/bigdata/zlin/Melanoma_meta/figures/Abundance/freq_pre.pdf', height = 8, width = 10)
Heatmap(mat_pre, na_col = 'white', name = 'Pre', column_title = 'Pre', 
        column_title_gp = gpar(fontsize = 20, fontface = 'bold'),
        show_column_names = F, show_row_names = F,
        row_order = match(celltype_order_r2, rownames(mat_pre)),
        cluster_rows = F, cluster_columns = F,
        col = rev(mako(n = 10)),
        column_names_gp = gpar(fontsize = 8),
        top_annotation = col_ha,
        left_annotation = row_ha,
        row_names_gp = gpar(fontsize = 6),
        width = ncol(mat_pre)*unit(0.1, "cm"), 
        height = nrow(mat_pre)*unit(2.3, "mm"),
        heatmap_legend_param = list(
          title = "Freq"
        )) 
annotation_titles = c(dataset = "Dataset",
                      cancertype = "Cancer type",
                      treatment = "Treatment",
                      prior = "Prior-Tx",
                      response = "Response",
                      res_metric = "Response Metric")
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
dev.off()

pdf('/bigdata/zlin/Melanoma_meta/figures/Abundance/freq_post.pdf', height = 8, width = 10)
Heatmap(mat_post, na_col = 'white', name = 'Post', column_title = 'Post', 
        column_title_gp = gpar(fontsize = 20, fontface = 'bold'),
        show_column_names = F, 
        row_order = match(celltype_order_r2, rownames(mat_post)),
        cluster_rows = F, cluster_columns = F,
        col = rev(mako(n = 10)),
        column_names_gp = gpar(fontsize = 8),
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 6),
        width = ncol(mat_post)*unit(0.1, "cm"), 
        height = nrow(mat_post)*unit(2.3, "mm"),
        heatmap_legend_param = list(
          title = "Freq"))
annotation_titles = c(dataset = "Dataset",
                      cancertype = "Cancer type",
                      treatment = "Treatment",
                      prior = "Prior-Tx",
                      response = "Response",
                      res_metric = "Response Metric")
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    # grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
dev.off()

# R(o/e)
meta <- roie_mat |> 
  distinct(patient, .keep_all = T) |> 
  select(!c(time:ratio))
roie_mat$ratio[roie_mat$ratio>2] = 2
timepoint <- 'Pre'
mat_pre <- roie_mat |> 
  filter(time == timepoint) |>
  select(celltype, patient, ratio) |> 
  pivot_wider(names_from = celltype, values_from = ratio) |> 
  column_to_rownames(var = 'patient') |> 
  t() 
timepoint <- 'Post'
mat_post <- roie_mat |> 
  filter(time == timepoint) |>
  select(celltype, patient, ratio) |> 
  pivot_wider(names_from = celltype, values_from = ratio) |> 
  column_to_rownames(var = 'patient') |> 
  t() 
# column annotation
col_ha = HeatmapAnnotation(
  dataset = meta$dataset,
  cancertype = meta$cancertype,
  treatment = meta$treatment,
  prior = meta$prior,
  response = meta$response,
  res_metric = meta$res_metric,
  col = list(dataset = structure(names = unique(meta$dataset), dittoColors()[1:length(unique(meta$dataset))]),
             cancertype = structure(names = unique(meta$cancertype), pal_jco()(length(unique(meta$cancertype)))),
             treatment = structure(names = unique(meta$treatment), pal_nejm()(length(unique(meta$treatment)))),
             prior = structure(names = unique(meta$prior), pal_rickandmorty()(length(unique(meta$prior)))),
             response = structure(names = unique(meta$response), pal_startrek()(length(unique(meta$response)))),
             res_metric = structure(names = unique(meta$res_metric), pal_d3()(length(unique(meta$res_metric))))
             ),
  annotation_legend_param = list(
    dataset = list(title = "Dataset"),
    cancertype = list(title = "Cancer type"),
    treatment = list(title = "Treatment"),
    prior = list(title = "Prior-Tx"),
    response = list(title = "Response"),
    res_metric = list(title = "Response Metric")),
  show_annotation_name = FALSE
)
# row annotation
color_main <- c('#E31A1C', '#1F78B4', '#A6CEE3', '#CAB2D6', '#6A3D9A', '#FF7F00', '#B15928', '#FFFF99', '#B2DF8A', '#33A02C', '#FB9A99', '#FDBF6F')
row_ha = rowAnnotation(
  celltype_main = df_anno$celltype_main[match(rownames(mat_pre), df_anno$celltype_r2)],
  col = list(celltype_main = structure(names = c('CD4+T', 'CD8+T', 'NK', 'B', 'Plasma', 'pDC', 'Mast', 'cDC', 'Mono', 'Macro', 'Endo', 'CAF'), color_main)),
  annotation_legend_param = list(celltype_main = list(title = "Main")),
  show_annotation_name = FALSE
)
pdf('/bigdata/zlin/Melanoma_meta/figures/Abundance/roie_pre.pdf', height = 8, width = 10)
Heatmap(mat_pre, na_col = 'white', name = 'Pre', column_title = 'Pre', 
        column_title_gp = gpar(fontsize = 20, fontface = 'bold'),
        show_column_names = F, show_row_names = F,
        row_order = match(celltype_order_r2, rownames(mat_pre)),
        cluster_rows = F, cluster_columns = F,
        col = rev(brewer.pal(10, 'RdBu')),
        column_names_gp = gpar(fontsize = 8),
        top_annotation = col_ha,
        left_annotation = row_ha,
        row_names_gp = gpar(fontsize = 6),
        width = ncol(mat_pre)*unit(0.1, "cm"), 
        height = nrow(mat_pre)*unit(2.3, "mm"),
        heatmap_legend_param = list(
          title = "R(o/e)", at = c(0, 1, 2), 
          labels = c("0", "1", ">2")
)) 
annotation_titles = c(dataset = "Dataset",
                      cancertype = "Cancer type",
                      treatment = "Treatment",
                      prior = "Prior-Tx",
                      response = "Response",
                      res_metric = "Response Metric")
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
dev.off()

pdf('/bigdata/zlin/Melanoma_meta/figures/Abundance/roie_post.pdf', height = 8, width = 10)
Heatmap(mat_post, na_col = 'white', name = 'Post', column_title = 'Post', 
        column_title_gp = gpar(fontsize = 20, fontface = 'bold'),
        show_column_names = F, 
        row_order = match(celltype_order_r2, rownames(mat_post)),
        cluster_rows = F, cluster_columns = F,
        col = rev(brewer.pal(10, 'RdBu')),
        column_names_gp = gpar(fontsize = 8),
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 6),
        width = ncol(mat_post)*unit(0.1, "cm"), 
        height = nrow(mat_post)*unit(2.3, "mm"),
        heatmap_legend_param = list(
          title = "R(o/e)", at = c(0, 1, 2), 
          labels = c("0", "1", "2")))
annotation_titles = c(dataset = "Dataset",
                      cancertype = "Cancer type",
                      treatment = "Treatment",
                      prior = "Prior-Tx",
                      response = "Response",
                      res_metric = "Response Metric")
for(an in names(annotation_titles)) {
  decorate_annotation(an, {
    # grid.text(annotation_titles[an], unit(-2, "mm"), just = "right")
    grid.rect(gp = gpar(fill = NA, col = "black"))
  })
}
dev.off()

# distribution across dataset
t_nk <- c('CD4_Naive','CD4_Tm_TNF','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM',
          'CD4_Tm_CCL5', 'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
          'CD4_Treg_Early', 'CD4_Treg_ISG', 'CD4_Treg_TNFRSF9', 'CD4_Th_ISG', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif',
          'CD8_Naive', 'CD8_Prolif', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
          'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13', 'CD8_Tex_OXPHOS-', 
          'CD8_ISG', 'CD8_Temra_CX3CR1', 'CD8_NK-like', 'MAIT', 'gdT', 'NK_CD56loCD16hi', 'NK_CD56hiCD16lo')
bplasma <- c('Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Exhausted B cells', 'Plasma')
mye <- c('Mast','pDC','cDC1', 'cDC2_CD1C', 'cDC2_IL1B', 'cDC2_ISG15', 'cDC2_CXCL9', 'DC_LC-like', 'cCD2_MoDC', 'DC_Migr',
         'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
         'Macro_NLRP3','Macro_ISG15', 'Macro_INHBA', 'Macro_FN1', 'Macro_SPP1',
         'Macro_LYVE1','Macro_IL1B','Macro_C1QC','Macro_TREM2')
nonimmune <- c('EC_lymphatic','EC_vascular','EndMT','CAF_inflammatory', 'CAF_adipogenic', 'CAF_PN', 'CAF_AP', 'Myofibroblast')
order_row <- c(t_nk, bplasma, mye, nonimmune)
# Pre
dist_median <- roie_mat |> 
  filter(time == 'Pre') |> 
  group_by(dataset) |> 
  mutate(n_patients = n_distinct(patient)) |> 
  filter(n_patients >= 5) |> 
  group_by(dataset, celltype) |> 
  summarize(median_roie = median(ratio)) |> 
  pivot_wider(values_from = median_roie, names_from = dataset) |> column_to_rownames(var = 'celltype')
overall <- roie_mat |> 
  filter(time == 'Pre') |> 
  group_by(dataset) |> 
  mutate(n_patients = n_distinct(patient)) |> 
  filter(n_patients >= 5) |> 
  group_by(celltype) |> 
  summarize(median_roie = median(ratio))
dist_median$`All Datasets` <- overall$median_roie[match(rownames(dist_median), overall$celltype)]
dist_median[dist_median == 0] <- NA

pdf('/bigdata/zlin/Melanoma_meta/figures/Abundance/roie_pre_median.pdf', height = 10, width = 6)
pheatmap(dist_median[order_row, c('All Datasets','SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'CRC_Li', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_IMCISION', 'HNSC_Luoma', 'NSCLC_Liu')], 
         cluster_rows = F, cluster_cols = F,
         breaks = c(seq(0, 2, by=0.5)), 
         border_color = NA, color = rev(brewer.pal(8,"PiYG")),
         fontsize_row = 7, fontsize_col = 8,
         name = 'R(o/e)', main = 'Pre (median)', display_numbers = T,
         gaps_row = c(19, 32, 34, 36, 41, 51, 63))
dev.off()

# Post
dist_median <- roie_mat |> 
  filter(time == 'Post') |> 
  group_by(dataset) |> 
  mutate(n_patients = n_distinct(patient)) |> 
  filter(n_patients >= 5) |> 
  group_by(dataset, celltype) |> 
  summarize(median_roie = median(ratio)) |> 
  pivot_wider(values_from = median_roie, names_from = dataset) |> column_to_rownames(var = 'celltype')
overall <- roie_mat |> 
  filter(time == 'Post') |> 
  group_by(dataset) |> 
  mutate(n_patients = n_distinct(patient)) |> 
  filter(n_patients >= 5) |> 
  group_by(celltype) |> 
  summarize(median_roie = median(ratio))
dist_median$`All Datasets` <- overall$median_roie[match(rownames(dist_median), overall$celltype)]
dist_median[dist_median == 0] <- NA

pdf('/bigdata/zlin/Melanoma_meta/figures/Abundance/roie_post_median.pdf', height = 10, width = 6)
pheatmap(dist_median[order_row, c('All Datasets','SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'CRC_Li', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_IMCISION', 'HNSC_Luoma', 'NSCLC_Liu')], 
         cluster_rows = F, cluster_cols = F,
         breaks = c(seq(0, 2, by=0.5)), 
         border_color = NA, color = rev(brewer.pal(8,"PiYG")),
         fontsize_row = 7, fontsize_col = 8,
         name = 'R(o/e)', main = 'Post (median)', display_numbers = T,
         gaps_row = c(19, 32, 34, 36, 41, 51, 63))
dev.off()

# # Univariate fixed-effect model
# uni_fe <- function(roie_mat , meta_combi, condition, values, n.sample = 15, covariates = c('ratio', 'interval', 'cancertype', 'modality', 'prior')){
#   df_check <- meta_combi |>
#     filter(patient %in% unique(roie_mat$patient)) |> 
#     select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
#     distinct(celltype_r2, patient, time_point, .keep_all = T) |>
#     pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
#     filter(Pre >= 5 | Post >= 5)
#   get_subtype <- function(value){
#     df_check |>
#       filter(eval(parse(text = sprintf("%s == '%s'", condition, value)))) |> 
#       group_by(celltype_r2) |>
#       summarise(count = n()) |> 
#       filter(count >= n.sample) |>
#       pull(celltype_r2)
#   }
#   if (is.na(condition) == FALSE){
#     test_list <- list()
#     subtypes <- values |> map(get_subtype) |> reduce(intersect)
#     for (value in values){
#       print(value)
#       roie <- filter(roie_mat, eval(parse(text = condition)) == value)
#       uni_models <- lapply(subtypes, function(subtype){
#         print(subtype)
#         pt <- filter(df_check, celltype_r2 == subtype, eval(parse(text = condition)) == value) |> pull(patient)
#         df <- filter(roie, celltype == subtype, patient %in% pt) 
#         df$timepoint <- ifelse(df$time == 'Pre', 0, 1)
#         valid_covariates <- sapply(covariates, function(x) length(unique(df[[x]])) > 1)
#         filtered_covariates <- covariates[valid_covariates]
#         formula <- as.formula(paste0("timepoint ~", paste(c(filtered_covariates, "(1 | patient)"), collapse = " + ")))
#         model <- lmer(formula, data = df, REML = FALSE)
#         model_summary <- summary(model)
#         coef_table <- coef(summary(model))
#         confint_table <- confint(model, level = 0.95)
#         subtype_results <- coef_table['ratio', , drop = FALSE]
#         subtype_confint <- confint_table['ratio', , drop = FALSE]
#         combined_results <- data.frame(
#           Celltypes = subtype,
#           Estimate = subtype_results[1],
#           StdError = subtype_results[2],
#           tValue = subtype_results[4],
#           pValue = coef_table['ratio', 'Pr(>|t|)'],
#           CI_lower = subtype_confint[1],
#           CI_upper = subtype_confint[2],
#           group = value
#         )
#         return(combined_results)
#       })
#       results <- do.call(rbind, uni_models) |> data.frame()
#       results$fdr <- p.adjust(results$pValue, method = 'fdr', n = nrow(results))
#       test_list[[which(values == value)]] <- results 
#     }
#     df <- do.call(rbind, test_list)
#   } else {
#     subtypes <- 
#       df_check |> 
#       group_by(celltype_r2) |>
#       summarise(count = n()) |> 
#       filter(count >= 15) |>
#       pull(celltype_r2)
#     uni_models <- lapply(subtypes, function(subtype){
#       print(subtype)
#       roie_mat$timepoint <- ifelse(roie_mat$time == 'Pre', 0, 1)
#       formula <- as.formula(paste("timepoint ~ ratio + interval + cancertype + modality + prior + (1 | patient)"))
#       model <- lmer(formula, data = roie_mat[roie_mat$celltype == subtype,], REML = FALSE)
#       model_summary <- summary(model)
#       coef_table <- coef(summary(model))
#       confint_table <- confint(model, level = 0.95)
#       subtype_results <- coef_table['ratio', , drop = FALSE]
#       subtype_confint <- confint_table['ratio', , drop = FALSE]
#       combined_results <- data.frame(
#         Celltypes = subtype,
#         Estimate = subtype_results[1],
#         StdError = subtype_results[2],
#         tValue = subtype_results[4],
#         pValue = coef_table['ratio', 'Pr(>|t|)'],
#         CI_lower = subtype_confint[1],
#         CI_upper = subtype_confint[2]
#       )
#       return(combined_results)
#     })
#     df <- do.call(rbind, uni_models) |> data.frame()
#     df$fdr <- p.adjust(df$pValue, method = 'fdr', n = nrow(df))
#   }
#   celltype_keep <- filter(df, fdr < 0.05) |> pull(Celltypes) |> unique()
#   df <- filter(df, Celltypes %in% celltype_keep)
#   return(df)
# }
# # Overall
# condition <- NA
# res_fe <- uni_fe(roie_mat, meta_combi, condition, values); res_fe
# order_row <- c(cd4t, cd8t, nk, bplasma, mye, endo, caf)[c(cd4t, cd8t, nk, bplasma, mye, endo, caf) %in% unique(res_fe_mono$Celltypes)]
# pdf('/bigdata/zlin/Melanoma_meta/figures/uni_fe_all.pdf', height = 5, width = 5)
# p <- ggplot(res_fe, aes(x= factor(Celltypes, levels = rev(order_row)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(fdr))) +
#   geom_linerange(size=1, color = '#A4DBFF', position=position_dodge(width = 0.5)) +
#   geom_hline(yintercept=0, lty=2) +
#   geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
#   scale_size(breaks = c(-log10(0.001), -log10(0.01), -log10(0.05)), labels = c('<0.001', '<0.01', '<0.05'), range = c(2,5)) +
#   scale_y_continuous(name= "Effect Size", limits = c(-0.4, 0.5)) + xlab("") +
#   coord_flip() + theme_minimal() +
#   labs(size = "FDR") +
#   guides(size = guide_legend(override.aes = list(fill = "black"))) +
#   theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
#         axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
#         axis.text.y = element_text(size = 9, colour = "black"))
# grid.draw(p)
# grid.segments(
#   x = unit(0.4, "npc"), 
#   x1 = unit(0.3, "npc"), 
#   y = unit(0.22, "npc"), 
#   y1 = unit(0.22, "npc"), 
#   arrow = arrow(type = "open", length = unit(0.05, "inches"))
# )
# grid.segments(
#   x = unit(0.7, "npc"), 
#   x1 = unit(0.8, "npc"), 
#   y = unit(0.22, "npc"), 
#   y1 = unit(0.22, "npc"), 
#   arrow = arrow(type = "open", length = unit(0.05, "inches"))
# )
# grid.text("Pre", x = unit(0.35, "npc"), y = unit(0.19, "npc"), gp = gpar(fontsize = 8))
# grid.text("Post", x = unit(0.75, "npc"), y = unit(0.19, "npc"), gp = gpar(fontsize = 8))
# dev.off()
# 
# # Response
# condition <- 'response'
# values <- unique(eval(parse(text = paste0('roie_mat$', condition))))[1:2]; values
# res_fe <- uni_fe(roie_mat, meta_combi, condition, values); res_fe
# order_row <- c(cd4t, cd8t, nk, bplasma, mye, endo, caf)[c(cd4t, cd8t, nk, bplasma, mye, endo, caf) %in% unique(res_fe_mono$Celltypes)]
# pdf('/bigdata/zlin/Melanoma_meta/figures/uni_fe_res.pdf', height = 6, width = 5)
# #define colours for dots and bars y
# dotCOLS = c("#a6d8f0","#DB5F75")
# barCOLS = c("#0555A8","#B10041")
# p <- ggplot(res_fe, aes(x= factor(Celltypes, levels = rev(order_row)), y=Estimate, ymin=CI_lower, ymax=CI_upper,col=group,fill=group, size = -log10(fdr))) +
#   #specify position here
#   geom_linerange(size=1,position=position_dodge(width = 0.5)) +
#   geom_hline(yintercept=0, lty=2) +
#   #specify position here too
#   geom_point(shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.5)) +
#   scale_size(breaks = c(-log10(0.001), -log10(0.01), -log10(0.05)), labels = c('<0.001', '<0.01', '<0.05'), range = c(0,4)) +
#   scale_fill_manual(values=barCOLS) +
#   scale_color_manual(values=dotCOLS) +
#   scale_y_continuous(name= "Effect Size", limits = c(-0.7, 0.9)) + xlab("") +
#   coord_flip() +
#   theme_minimal() +
#   labs(fill = "Response", color = "Response", size = "FDR") +
#   guides(size = guide_legend(override.aes = list(fill = "black"))) +
#   theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
#         axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
#         axis.text.y = element_text(size = 9, colour = "black"))
# grid.draw(p)
# grid.segments(
#   x = unit(0.3, "npc"), 
#   x1 = unit(0.2, "npc"), 
#   y = unit(0.13, "npc"), 
#   y1 = unit(0.13, "npc"), 
#   arrow = arrow(type = "open", length = unit(0.05, "inches"))
# )
# # Add right arrow
# grid.segments(
#   x = unit(0.7, "npc"), 
#   x1 = unit(0.8, "npc"), 
#   y = unit(0.13, "npc"), 
#   y1 = unit(0.13, "npc"), 
#   arrow = arrow(type = "open", length = unit(0.05, "inches"))
# )
# grid.text("Pre", x = unit(0.25, "npc"), y = unit(0.11, "npc"), gp = gpar(fontsize = 8))
# grid.text("Post", x = unit(0.75, "npc"), y = unit(0.11, "npc"), gp = gpar(fontsize = 8))
# dev.off()
# 
# # Subgroup (Monotherapy)
# df_mono <- roie_mat |> filter(modality == 'Mono')
# df_meta_mono <- meta_combi |> filter(modality == 'Mono')
# condition <- 'response'
# values <- c("NR", "RE")
# res_fe_mono <- uni_fe(df_mono, df_meta_mono, condition, values, n.sample = 10)
# order_row <- c(cd4t, cd8t, nk, bplasma, mye, endo, caf)[c(cd4t, cd8t, nk, bplasma, mye, endo, caf) %in% unique(res_fe_mono$Celltypes)]
# pdf('/bigdata/zlin/Melanoma_meta/figures/uni_fe_mono.pdf', height = 6, width = 5)
# dotCOLS = c("#CAB2D6","#FDBF6F")
# barCOLS = c("#6A3D9A","#FF7F00")
# p <- ggplot(res_fe_mono, aes(x= factor(Celltypes, levels = rev(order_row)), y=Estimate, ymin=CI_lower, ymax=CI_upper,col=group, fill=group, size = -log10(fdr))) +
#   #specify position here
#   geom_linerange(size=1,position=position_dodge(width = 0.5)) +
#   geom_hline(yintercept=0, lty=2) +
#   #specify position here too
#   geom_point(shape=21, colour="white", stroke = 0.5, position=position_dodge(width = 0.5)) +
#   scale_size(breaks = c(-log10(0.001), -log10(0.01), -log10(0.05)), labels = c('<0.001', '<0.01', '<0.05'), range = c(0,4)) +
#   scale_fill_manual(values=barCOLS) +
#   scale_color_manual(values=dotCOLS) +
#   scale_y_continuous(name= "Effect Size", limits = c(-1, 1)) + xlab("") +
#   coord_flip() + 
#   theme_minimal() +
#   labs(fill = "Response", color = "Response", size = "FDR") +
#   guides(size = guide_legend(override.aes = list(fill = "black"))) +
#   theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
#         axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
#         axis.text.y = element_text(size = 9, colour = "black"))
# grid.draw(p)
# grid.segments(
#   x = unit(0.4, "npc"), 
#   x1 = unit(0.3, "npc"), 
#   y = unit(0.13, "npc"), 
#   y1 = unit(0.13, "npc"), 
#   arrow = arrow(type = "open", length = unit(0.05, "inches"))
# )
# # Add right arrow
# grid.segments(
#   x = unit(0.7, "npc"), 
#   x1 = unit(0.8, "npc"), 
#   y = unit(0.13, "npc"), 
#   y1 = unit(0.13, "npc"), 
#   arrow = arrow(type = "open", length = unit(0.05, "inches"))
# )
# grid.text("Pre", x = unit(0.35, "npc"), y = unit(0.11, "npc"), gp = gpar(fontsize = 8))
# grid.text("Post", x = unit(0.75, "npc"), y = unit(0.11, "npc"), gp = gpar(fontsize = 8))
# dev.off()

# Univariate Logistic
meta_combi$int_cat <- ifelse(meta_combi$interval <= 21, 'Short','Long')
roie_mat$int_cat <- ifelse(roie_mat$interval <= 21, 'Short','Long')
uni_logi <- function(roie_mat , meta_combi, n = 5){
  df_check <- meta_combi |>
    filter(patient %in% unique(roie_mat$patient)) |>
    select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
    distinct(celltype_r2, patient, time_point, .keep_all = T) |>
    pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) 
    # filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5))
  subtypes <- df_check |>
    group_by(celltype_r2, response) |>
    mutate(n_res = n()) |>
    filter((response == 'RE' & n_res >= n) | (response == 'NR' & n_res >= n)) |>
    distinct(celltype_r2, .keep_all = T) |>
    select(celltype_r2, response) |>
    group_by(celltype_r2) |>
    summarize(n = n()) |>
    filter(n == 2) |>
    pull(celltype_r2)
  uni_models <- lapply(subtypes, function(subtype){
    tryCatch({
      print(subtype)
      pt <- filter(df_check, celltype_r2 == subtype) |> pull(patient)
      df <- filter(roie_mat, celltype == subtype, patient %in% pt, time == 'Pre')
      df$response <- ifelse(df$response == 'RE', 1, 0)
      variables_to_use <- c('ratio', 'modality'[length(unique(df$modality)) > 1], 'int_cat'[length(unique(df$int_cat)) > 1], 'cancertype'[length(unique(df$cancertype)) > 1])
      formula <- as.formula(paste("response ~", paste((variables_to_use), collapse = " + ")))
      model <- glm(formula, data = df, family = binomial)
      coef_table <- coef(summary(model))
      confint_table <- confint(model, level = 0.95)
      subtype_results <- coef_table['ratio', , drop = FALSE]
      subtype_confint <- confint_table['ratio', , drop = FALSE]
      combined_results <- data.frame(
        Celltypes = subtype,
        Estimate = subtype_results[1],
        StdError = subtype_results[2],
        pValue = coef_table['ratio', 'Pr(>|z|)'],
        CI_lower = subtype_confint[1],
        CI_upper = subtype_confint[2]
      )
      return(combined_results)
      },  error = function(e) {
      # Simply return NULL if an error occurs
      NULL
      })
    })
  # Filter out NULL elements from the list
  uni_models <- Filter(Negate(is.null), uni_models)
  results <- do.call(rbind, uni_models) |> data.frame()
  results$fdr <- p.adjust(results$pValue, method = 'fdr', n = nrow(results))
  return(results)
}
res_logi <- uni_logi(roie_mat, meta_combi); res_logi
res_filt <- filter(res_logi, pValue < 0.5) |> arrange(desc(Estimate))
order_row <- res_filt$Celltypes

pdf('/bigdata/zlin/Melanoma_meta/figures/uni_logi_pre.pdf', height = 6, width = 6)
p <- ggplot(res_filt, aes(x= factor(Celltypes, levels = rev(order_row)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  #specify position here
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  #specify position here too
  geom_point(shape=21, colour="white", fill = '#005D89', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.001), -log10(0.01), -log10(0.05)), labels = c('<0.001', '<0.01', '<0.05'), range = c(0,6)) +
  scale_y_continuous(name= "Log(Odd Ratio)", limits = c(-5, 4), breaks = c(-5, -2.5, 0, 2.5), labels = c("<-5", "-2.5", "0", "5")) + 
  xlab("") + ggtitle("R o/e (Pre-Tx)") +
  coord_flip() + 
  theme_minimal() + 
  labs(size = "P value") +
  guides(size = guide_legend(override.aes = list(fill = "#005D89"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))

grid.draw(p)
grid.segments(
  x = unit(0.45, "npc"), 
  x1 = unit(0.35, "npc"), 
  y = unit(0.13, "npc"), 
  y1 = unit(0.13, "npc"), 
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
# Add right arrow
grid.segments(
  x = unit(0.7, "npc"), 
  x1 = unit(0.8, "npc"), 
  y = unit(0.13, "npc"), 
  y1 = unit(0.13, "npc"), 
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Non-responder", x = unit(0.4, "npc"), y = unit(0.11, "npc"), gp = gpar(fontsize = 8))
grid.text("Responder", x = unit(0.75, "npc"), y = unit(0.11, "npc"), gp = gpar(fontsize = 8))
dev.off()

p1 <- filter(roie_mat, celltype == 'Plasma') |> 
  select(time, ratio, patient) |> 
  pivot_wider(names_from = time, values_from = ratio) |> 
  ggplot(aes(x = Pre, y = Post)) +
  geom_point() + geom_smooth(method = lm) +
  ggtitle('Plasma')

p2 <- filter(roie_mat, celltype == 'CAF_inflammatory') |> 
  select(time, ratio, patient) |> 
  pivot_wider(names_from = time, values_from = ratio) |> 
  ggplot(aes(x = Pre, y = Post)) +
  geom_point() + geom_smooth(method = lm) +
  ggtitle('CAF_inflammatory')

p3 <- filter(roie_mat, celltype == 'Macro_C1QC') |> 
  select(time, ratio, patient) |> 
  pivot_wider(names_from = time, values_from = ratio) |> 
  ggplot(aes(x = Pre, y = Post)) +
  geom_point() + geom_smooth(method = lm) +
  ggtitle('Macro_C1QC')
p1+p2+p3

# Pairwise comparison (frequency)
meta_combi$time_point[meta_combi$time_point == 'Post' & meta_combi$interval < 21] <- 'On'
meta_combi$int_cat <- ifelse(meta_combi$interval < 21, 'Short','Long')
meta_combi$time_point <- factor(meta_combi$time_point, levels = c('Pre', 'On', 'Post'))

# RE
# T-test early
meta_re_early <- meta_combi |>  
  filter(celltype_r2 == subtype, response == 'RE', count_r2 >= 5, int_cat == 'Early') |> 
  distinct(sample, .keep_all = T) |> 
  group_by(patient) |> 
  mutate(pt_count = n()) |> 
  ungroup() |> 
  filter(pt_count == 2) |> 
  mutate(freq_rela = freq_r2/freq_main) |> 
  select(freq_rela, time_point, patient, int_cat) |> 
  pivot_wider(values_from = freq_rela, names_from = time_point)
stat_res_re_early <- t.test(meta_re_early$Pre, meta_re_early$On, paired = T)
# T-test post
meta_re_post <- meta_combi |>  
  filter(celltype_r2 == subtype, response == 'RE', count_r2 >= 5, int_cat == 'Post') |> 
  distinct(sample, .keep_all = T) |> 
  group_by(patient) |> 
  mutate(pt_count = n()) |> 
  ungroup() |> 
  filter(pt_count == 2) |> 
  mutate(freq_rela = freq_r2/freq_main) |> 
  select(freq_rela, time_point, patient, int_cat) |> 
  pivot_wider(values_from = freq_rela, names_from = time_point)
stat_res_re_post <- t.test(meta_re_post$Pre, meta_re_post$Post, paired = T)
# boxplot RE
meta_combi |> 
  filter(celltype_r2 == subtype, response != 'NE', response == 'RE', count_r2 >= 5) |> 
  distinct(sample, .keep_all = T) |> 
  group_by(patient) |> 
  mutate(pt_count = n()) |> 
  ungroup() |> 
  filter(pt_count == 2) |> 
  mutate(freq_rela = freq_r2/freq_main) |> 
  ggplot(aes(x = time_point, y = freq_rela)) +
  geom_boxplot(color = 'black', outlier.shape = NA) +
  geom_line(aes(group = patient), color = "gray",linetype = "dashed") +
  geom_point(aes(color = cancertype, shape = modality), size=3, alpha = 0.6) +
  scale_color_manual(values = brewer.pal(length(unique(meta_combi$dataset)),"Set1")) +
  theme_classic2() + xlab("") + ylab("Frequency") + ggtitle('RE') +
  labs(color = "Cancer Type", shape = "Modality") +
  geom_text(aes(label = paste0("p.adj = ", formatC(p.adjust(stat_res_re_early$p.value, method = 'fdr', 2), digits = 3)), 
                x = 'Post', y = 0.5), vjust = -0.5) +
  geom_text(aes(label = paste0("p.adj = ", formatC(p.adjust(stat_res_re_post$p.value, method = 'fdr', 2), digits = 3)), 
                x = 'On', y = 0.5), vjust = -0.5)
# NR
# T-test early
meta_re_early <- meta_combi |>  
  filter(celltype_r2 == subtype, response == 'NR', count_r2 >= 5, int_cat == 'Early') |> 
  distinct(sample, .keep_all = T) |> 
  group_by(patient) |> 
  mutate(pt_count = n()) |> 
  ungroup() |> 
  filter(pt_count == 2) |> 
  mutate(freq_rela = freq_r2/freq_main) |> 
  select(freq_rela, time_point, patient, int_cat) |> 
  pivot_wider(values_from = freq_rela, names_from = time_point)
stat_res_re_early <- t.test(meta_re_early$Pre, meta_re_early$On, paired = T)
# T-test post
meta_re_post <- meta_combi |>  
  filter(celltype_r2 == subtype, response == 'NR', count_r2 >= 5, int_cat == 'Post') |> 
  distinct(sample, .keep_all = T) |> 
  group_by(patient) |> 
  mutate(pt_count = n()) |> 
  ungroup() |> 
  filter(pt_count == 2) |> 
  mutate(freq_rela = freq_r2/freq_main) |> 
  select(freq_rela, time_point, patient, int_cat) |> 
  pivot_wider(values_from = freq_rela, names_from = time_point)
stat_res_re_post <- t.test(meta_re_post$Pre, meta_re_post$Post, paired = T)
# boxplot NR
meta_combi |> 
  filter(celltype_r2 == subtype, response != 'NE', response == 'NR', count_r2 >= 5) |> 
  distinct(sample, .keep_all = T) |> 
  group_by(patient) |> 
  mutate(pt_count = n()) |> 
  ungroup() |> 
  filter(pt_count == 2) |> 
  mutate(freq_rela = freq_r2/freq_main) |> 
  ggplot(aes(x = time_point, y = freq_rela)) +
  geom_boxplot(color = 'black', outlier.shape = NA) +
  geom_line(aes(group = patient), color = "gray",linetype = "dashed") +
  geom_point(aes(color = cancertype, shape = modality), size=3, alpha = 0.6) +
  scale_color_manual(values = brewer.pal(length(unique(meta_combi$dataset)),"Set1")) +
  theme_classic2() + xlab("") + ylab("Frequency") + ggtitle('RE') +
  labs(color = "Cancer Type", shape = "Modality") +
  geom_text(aes(label = paste0("p.adj = ", formatC(p.adjust(stat_res_nr_early$p.value, method = 'fdr', 2), digits = 3)), 
                x = 'Post', y = 0.5), vjust = -0.5) +
  geom_text(aes(label = paste0("p.adj = ", formatC(p.adjust(stat_res_nr_post$p.value, method = 'fdr', 2), digits = 3)), 
                x = 'On', y = 0.5), vjust = -0.5)

# Filter subtypes
n.cell <- 5
df_check <- meta_combi |>
  select(celltype_r2, patient, sample, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
  distinct(celltype_r2, sample, time_point, .keep_all = T) |>
  pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
  filter(Pre >= n.cell | On >= n.cell| Post >= n.cell)
process_subtype <- function(subtype, df_check){
  df <- filter(df_check, celltype_r2 == subtype) |> 
    group_by(patient) |> 
    mutate(pt_count = n()) |> 
    ungroup() |> 
    filter(pt_count == 2) |> 
    distinct(patient, .keep_all = T) |> 
    filter((sum((response == "RE" & interval < 21) >= 3) | (sum(response == "RE" & interval >= 21) >= 3)) &
             ((sum(response == "NR" & interval < 21) >= 3) | (sum(response == "NR" & interval >= 21) >= 3))) 
}
results <- lapply(unique(df_check$celltype_r2), process_subtype, df_check)
combined_results <- do.call(rbind, results)
subtypes <- unique(combined_results$celltype_r2)

filter_data <- function(subtype, res, cat, n.cell = 5) {
  return(
    meta_combi |>
      filter(celltype_r2 == subtype, response == res, int_cat == cat) |>
      distinct(sample, .keep_all = TRUE) |>
      group_by(patient, sufficient_count_r2 = count_r2 >= n.cell) |> 
      mutate(pt_count = n()) |> 
      ungroup() |> 
      filter(pt_count == 2, sufficient_count_r2)
  )
}
responses <- c('RE', 'NR')
int_cats <- c('Short', 'Long')

pvalues <- data.frame(matrix(NA, nrow = length(subtypes), ncol = 4, 
                             dimnames = list(subtypes, paste0(rep(responses, each = 2), '_', int_cats))))
effectsizes <- data.frame(matrix(NA, nrow = length(subtypes), ncol = 4, 
                                 dimnames = list(subtypes, paste0(rep(responses, each = 2), '_', int_cats))))

pw_comp <- function(meta, int_cat){
  if (nrow(meta) < 6) {
    return(list(pvalue = NA, effectsize = NA))
  } else {
    meta_wide <- meta |>
      select(freq_r2, time_point, patient) |>
      pivot_wider(values_from = freq_r2, names_from = time_point)
  if (int_cat == 'Early'){
    ttest <- t.test(meta_wide$Pre, meta_wide$On, paired = TRUE)
    meta$time_point <- factor(meta$time_point, levels = c('Pre',  'On'))
    effectsize <- meta |> 
      filter(int_cat == 'Early') |> 
      cohens_d(freq_r2 ~ time_point, paired = TRUE)
    } else {
    ttest <- t.test(meta_wide$Pre, meta_wide$Post, paired = TRUE)
    meta$time_point <- factor(meta$time_point, levels = c('Pre',  'Post'))
    effectsize <- meta |> 
      filter(int_cat == 'Post') |> arrange(time_point, patient) |> 
      cohens_d(freq_r2 ~ time_point, paired = T)
    }
    pvalue <- ttest$p.value
    effectsize <- effectsize$effsize
    # if (ttest$estimate < 0) {effectsize <- effectsize * (-1)}
    return(list(pvalue = pvalue, effectsize = effectsize))
  }
}

for (subtype in subtypes) {
  for (int_cat in int_cats) {
    for (response in responses) {
      meta <- filter_data(subtype, res = response, cat = int_cat)
      result <- pw_comp(meta, int_cat)
      col_name <- paste(int_cat, response, sep = "_")
      pvalues[subtype, col_name] <- result$pvalue
      effectsizes[subtype, col_name] <- result$effectsize
    }
  }
}

pvalues <- filter(pvalues, RE_Early < 0.05 | RE_Post < 0.05 | NR_Early < 0.05 | NR_Post < 0.05)
pheatmap(effectsizes[rownames(pvalues),], 
         cluster_cols = F, cluster_rows = F,
         breaks = c(seq(-1, 1, by=0.5)))



# t-test
paired_comp <- function(roie_mat , meta_combi, condition, values){
  test_list <- list()
  df_check <- meta_combi |>
    select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
    distinct(celltype_r2, patient, time_point, .keep_all = T) |>
    pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
    filter(Pre >= 5 | Post >= 5)
  get_subtype <- function(value){
    df_check |>
      filter(eval(parse(text = sprintf("%s == '%s'", condition, value)))) |> 
      group_by(celltype_r2) |>
      summarise(count = n()) |> 
      filter(count >= 3) |>
      pull(celltype_r2)
  }
  subtypes <- values |> map(get_subtype) |> reduce(intersect)
  # paired t test with multiple hypothosis testing
  for (value in values){
    # calculate for each value in the condition
    roie <- filter(roie_mat, eval(parse(text = condition)) == value)
    dep_comp <- lapply(subtypes, function(subtype){
      pt <- filter(df_check, celltype_r2 == subtype, eval(parse(text = condition)) == value) |> pull(patient)
      df <- filter(roie, celltype == subtype, patient %in% pt) |>
        pivot_wider(names_from = time, values_from = ratio)
      df <- filter(roie, celltype == subtype, patient %in% pt)
      glm()
      t <- t.test(df$Post, df$Pre, paired = T, alternative = 'two.sided')
      cohen_d <- cohen.d(df$Post, df$Pre, paired = TRUE, hedges.correction = TRUE)
      comp <- c(subtype, t$statistic, t$p.value, cohen_d$estimate)
      return(comp)
    })
    results <- do.call(rbind, dep_comp) |> data.frame()
    rownames(results) <- subtypes
    colnames(results) <- c('celltype','t_score', 'pvalue','effect_size')
    results$fdr <- p.adjust(results$pvalue, method = 'fdr', n = nrow(results))
    test_list[[which(values == value)]] <- results 
  }
  names(test_list) <- values
  df <- reduce(test_list, full_join, by = "celltype", suffix = paste0('_', names(test_list)))
  return(df)
}
res_comp <- paired_comp(roie_mat, meta_combi, condition, values) |> .[.$fdr_RE < 0.05 | .$fdr_NR < 0.05,];res_comp

mat_es <- res_comp[,str_detect(colnames(res_comp), 'effect_size')] |> sapply(as.numeric)
rownames(mat_es) <- res_comp$celltype
colnames(mat_es) <- values
mat_sig <- res_comp[,str_detect(colnames(res_comp), 'fdr')] |> sapply(as.numeric)
mat_es[mat_es>0.8] = 0.8
mat_es[mat_es< -0.8] = -0.8

pdf('/bigdata/zlin/Melanoma_meta/figures/res_roie_0.05.pdf', height = 6, width = 6)
Heatmap(mat_es, name = "mat", col = rev(brewer.pal(n=8, name = 'RdBu')),
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "right", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Effect Size", at = c(-0.8, 0, 0.8), labels = c("-0.8","0","0.8")),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_es)*unit(7, "mm"), 
        height = nrow(mat_es)*unit(7, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
grid.text("P adjusted: \n *** < 0.001 \n** < 0.01 \n* < 0.05", x = 0.65, y = 0.33, gp = gpar(fontsize = 10))
dev.off()

# 
df_check <- meta_combi |>
  select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
  distinct(celltype_r2, patient, time_point, .keep_all = T) |>
  pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
  filter(Pre >= 5 | Post >= 5)
get_subtype <- function(value){
  df_check |>
    filter(eval(parse(text = sprintf("%s == '%s'", condition, value)))) |> 
    group_by(celltype_r2) |> 
    summarise(count = n()) |> 
    filter(count >= 3) |>
    pull(celltype_r2)
}

subtypes <- values |> map(get_subtype) |> reduce(intersect)
variables <- c('interval', 'cancertype', 'prior', 'modality')

# T test
test_list <- list()
for (value in values){
  results <- data.frame(
    Subtype = character(),
    coef = character(),
    pval = character(),
    stringsAsFactors = FALSE
  )
  # calculate for each value in the condition
  univ_logi <- lapply(subtypes, function(subtype){
    pt <- filter(df_check, celltype_r2 == subtype, eval(parse(text = condition)) == value) |> pull(patient)
    df <- filter(roie_mat, eval(parse(text = condition)) == value, celltype == subtype, patient %in% pt) 
    # Check if there are any variables in df with only one unique value and remove them
    variables_to_use <- c('ratio', variables[sapply(df[variables], function(x) length(unique(x))) > 1])

    # If no variables are left after the check, print a message and skip to the next iteration
    if (length(variables_to_use) == 0) {
      warning("No variables left to use in glm for subtype: ", subtype)
      next
    }
    
    # Create the formula for glm
    formula_glm <- as.formula(paste("time ~", paste((variables_to_use), collapse = " + ")))
    
    model <- glm(formula_glm, data = df, family = binomial)
    summary_model <- summary(model)
    
    # Extract the coefficients
    coefs <- coef(summary_model)

    # Save the results in the data frame
    results <- rbind(results, data.frame(
      Subtype = subtype,
      coef = coefs["ratio", "Estimate"],
      pval = coefs["ratio", "Pr(>|z|)"]
    ))
    return(results)
  })
  results <- do.call(rbind, univ_logi) |> data.frame()
  results$fdr <- p.adjust(results$pval, method = 'fdr', n = nrow(results))
  test_list[[which(values == value)]] <- results 
}

df <- reduce(test_list, full_join, by = "Subtype", suffix = paste0('_', values)) |> 
  .[.$fdr_RE < 0.05 | .$fdr_NR < 0.05,] |> 
  mutate(Category = case_when(
    fdr_RE < 0.05 & fdr_NR < 0.05 ~ "Both < 0.05",
    fdr_RE < 0.05 ~ "RE < 0.05",
    fdr_NR < 0.05  ~ "NR < 0.05"
  ))
df |> ggplot(aes(x = coef_RE, y = coef_NR, fill = Category)) +
  geom_point(aes(shape = Category), size = 3) +
  geom_text_repel(aes(label = Subtype), box.padding = 0.5, point.padding = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  theme_bw() 

mat_coef <- df[,str_detect(colnames(df), 'coef')] |> sapply(as.numeric)
rownames(mat_coef) <- df$Subtype
colnames(mat_coef) <- values
mat_sig <- df[,str_detect(colnames(df), 'fdr')] |> sapply(as.numeric)
mat_es[mat_es>5] = 5
mat_es[mat_es< -5] = -5
pdf('/bigdata/zlin/Melanoma_meta/figures/res_time_0.05.pdf', height = 10, width = 6)
Heatmap(mat_coef, name = "mat", col = rev(brewer.pal(n=8, name = 'RdBu')),
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "right", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Coef", at = c(-5, 0, 5), labels = c("-5","0","5")),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_coef)*unit(7, "mm"), 
        height = nrow(mat_coef)*unit(7, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
grid.text("P adjusted: \n *** < 0.001 \n ** < 0.01 \n* < 0.05", x = 0.65, y = 0.33, gp = gpar(fontsize = 10))
dev.off()



  
