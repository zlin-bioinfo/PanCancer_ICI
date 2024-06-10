rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','gtools','dittoSeq','RColorBrewer','ggpubr','tibble','epitools','dior','cowplot','ggnewscale','export','forcats','qs','purrr','effsize','corrplot','lmerTest','rstatix','pheatmap','janitor','ggmosaic','gridExtra','grid','ggalluvial','MetBrewer','ggmosaic','ComplexHeatmap','ggrepel')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# loading
meta_int <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_int.csv') 
pt_df <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_patient.csv') 
pt_df <- pt_df |> filter(dataset != 'NSCLC_Liu')
pt_df$dataset[pt_df$dataset %in% c('BCC_Yost','SCC_Yost')] <- 'BCC/SCC_Yost'
# Immune
immune_cell <- c('CD4_Naive','CD4_Tm_CREM-','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM', 'CD4_Tm_CCL5', 
                 'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
                 'CD4_Treg_Early', 'CD4_Treg_ISG15', 'CD4_Treg_TNFRSF9', 
                 'CD4_Th_ISG15', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif',
                 'CD8_Prolif', 'CD8_Naive', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
                 'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13',
                 'CD8_Tex_ISG15', 'CD8_Temra_CX3CR1', 'CD8_NK-like', 
                 'MAIT', 'gdT', 'NK_CD56loCD16hi', 'NK_CD56hiCD16lo',
                 'B_Naive', 'B_ISG15', 'B_HSP', 'B_MT2A', 'ACB_EGR1', 'ACB_NR4A2', 'ACB_CCR7', 'B_Memory', 'B_AtM',
                 'GCB_Pre', 'GCB_SUGCT', 'GCB_LMO2', 'GCB_Prolif', 'Plasmablast', 'Plasma_cell',
                 'Mast','pDC','cDC1', 
                 'cDC2_CD1C', 'cDC2_IL1B','cDC2_ISG15', 'cDC2_CXCL9', 'DC_LC-like', 'MigrDC', 'MoDC', 
                 'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
                 'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro_ISG15', 
                 'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2')
nonimmune <- c('EC_lymphatic','EC_vascular','EndMT','CAF_inflammatory', 'CAF_adipogenic', 'CAF_PN', 'CAF_AP', 'Myofibroblast')
# Heatmap
cor_mat <- function(meta_int, pt_df, timepoint, filtering = F){
  freq_filt <- meta_int |> 
    select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
    distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
    pivot_wider(values_from = count_r2, names_from = time_point, values_fill = 0)
  if (filtering == T){
    freq_filt <- freq_filt |> 
      filter(abs(Pre-Post) >= 3)
      # filter(abs(Pre-Post) >= 3, (Pre >= 3 | Post >= 3))
  }
  freq_filt$pt_r2 <- paste0(freq_filt$patient, '_', freq_filt$celltype_r2)
  freq_wide <- meta_int |> 
    select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
    distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
    pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
    mutate(change = log2((Post + 0.01)/(Pre + 0.01)), diff = (Post - Pre))
  freq_wide$pt_r2 <- paste0(freq_wide$patient, '_', freq_wide$celltype_r2)
  freq_wide <- filter(freq_wide, pt_r2 %in% freq_filt$pt_r2)
  mat_freq <- freq_wide |> 
    filter(dataset != 'NSCLC_Liu', 
           celltype_r2 %in% immune_cell) |> 
    select(patient, timepoint, celltype_r2) |> 
    pivot_wider(values_from = timepoint, names_from = patient, values_fill = 0) |> 
    column_to_rownames(var = 'celltype_r2')
  M <- cor(mat_freq)
  return(M)
}
pal <- colorRampPalette(brewer.pal(10, "RdBu"))
ht_cor_pt <- function(cor_matrix, n.cluster = 3, clust_method = 'complete'){
  pt_df$cluster <- cutree(hclust(dist(cor_matrix), method = clust_method), k=n.cluster)
  col_ha = HeatmapAnnotation(
    Dataset = pt_df$dataset,
    `Cancer Type` = pt_df$cancertype,
    Treatment = pt_df$treatment,
    Modality = pt_df$modality,
    `Response Metrics` = pt_df$res_metric,
    Response = pt_df$response,
    Cluster = pt_df$cluster,
    col = list(Dataset = structure(names = unique(pt_df$dataset), dittoColors()[1:length(unique(pt_df$dataset))]),
               `Cancer Type` = c("PCa" = "#a82203", 
                                 "CRC" = "#208cc0",     
                                 "NSCLC" = "#f1af3a",  
                                 "HNSC" = "#cf5e4e",
                                 "SCC" = "#3B7546",     
                                 "BCC" = "#0092A5", 
                                 "SKCM" = "#000000", 
                                 "TNBC" = "#FF0196", 
                                 "ER+BC" = "#F073EA",    
                                 "HER2+BC" = "#FFB0FF"),
               Treatment = structure(names = unique(pt_df$treatment), pal_npg()(length(unique(pt_df$treatment)))),
               Modality = c('Mono' = "#FF7F00", 'Dual' = "#6A3D9A"),
               Response = structure(names = c('RE','NR'), pal_startrek()(length(unique(pt_df$response)))),
               `Response Metrics` = structure(names = sort(unique(pt_df$res_metric)), met.brewer("Juarez", length(unique(pt_df$res_metric)))),
               Cluster = structure(names = unique(pt_df$cluster), pal_d3()(length(unique(pt_df$cluster))))
    )
    # annotation_legend_param = list(
    #   Dataset = list(title = "Dataset"),
    #   `Cancer Type` = list(title = "Cancer type"),
    #   treatment = list(title = "Treatment"),
    #   prior = list(title = "Prior-Tx"),
    #   res_metric = list(title = "Response Metric"),
    #   response = list(title = "Response"),
    #   cluster = list(title = "Cluster")),
    # show_annotation_name = T
  )
  Heatmap(cor_matrix, col = rev(pal(10)), name = "Correlation\n(Pearson's ρ)",
          show_row_names = F, show_column_names = F, 
          show_row_dend = F, show_column_dend = T,
          clustering_method_rows = clust_method,
          clustering_method_columns = clust_method,
          heatmap_legend_param = list(
            legend_direction = "horizontal", 
            legend_width = unit(2, "cm"), at = c(-1, 0, 1),
            legend_side = 'bottom',
            title_position = "topcenter"),
          width = ncol(cor_matrix)*unit(1, "mm"), 
          height = nrow(cor_matrix)*unit(1, "mm"),
          top_annotation = col_ha)
}

cor_matrix <- cor_mat(meta_int, pt_df, timepoint = 'diff', filtering = F) # c('Pre', 'diff', 'Post')
n <- 3
c_method <- 'ward.D2' # 'complete' (default) 'ward.D2'
ht_cor_pt(cor_matrix, n.cluster=n, clust_method = c_method)
pt_df$cluster_diff <- cutree(hclust(dist(cor_matrix), method = c_method), k=n)
pt_df$cluster_diff <- factor(pt_df$cluster_diff)
pt_df |> janitor::tabyl(response, cluster_diff)
test_result <- pt_df |> janitor::tabyl(response, cluster_diff) |> column_to_rownames(var = 'response') |> select(-2) |> chisq.test()
pt_df$response <- factor(pt_df$response, levels = c('RE','NR'))
pt_df |> 
  ggplot() +
  geom_mosaic(aes(x = product(response, cluster_diff), fill = response)) +
  labs(x = "", y = "") +
  scale_fill_brewer("Response", palette = "Set1") + theme_mosaic() +
  theme(axis.text.x = element_text(size = 14, colour = 'black'),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate("text", x = 1, y = 1.2, label = paste("Chi-squared Test (1 vs 3)\np-value:", round(test_result$p.value, 4)), # "Chi-square Test\np-value:"
           hjust = 1.05, vjust = 1.5, size = 4, colour = "black") +
  geom_mosaic_text(aes(x = product(response, cluster_diff), label = after_stat(.wt)), size = 3.5) 
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/mosaic_change.pdf', height = 5, width = 5)

pdf('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/ht_change.pdf', height = 10, width = 20)
ht_cor_pt(cor_matrix, n.cluster=n, clust_method = c_method)
dev.off()

# 
freq_wide <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(change = log2((Post + 0.01)/(Pre + 0.01)), diff = (Post - Pre))
freq_wide$cluster <- pt_df$cluster_diff[match(freq_wide$patient, pt_df$patient)]
freq_wide <- filter(freq_wide, !cluster  %in% c(2, NA))
freq_wide$group <- ifelse(freq_wide$cluster == 1, 'better', 'worse')
list_res <- lapply(unique(freq_wide$celltype_r2), function(subtype){
  # subtype <- 'Macro_IL1B'
  print(subtype)
  pvalue <- freq_wide |> 
    filter(celltype_r2 == subtype) |>
    distinct(patient, .keep_all = T) |> 
    t_test(diff ~ group) |> pull(p)
  effectsize <- freq_wide |> 
    filter(celltype_r2 == subtype) |> 
    distinct(patient, .keep_all = T) |> 
    cohens_d(diff~group) |> pull(effsize)
  combined_results <- data.frame(Celltypes = subtype,
                                 effsize = effectsize*1,
                                 pValue = pvalue)
  return(combined_results)
})
pdf('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/vol_subtypes.pdf', height = 5, width = 6)
do.call(rbind, list_res) |> 
  # mutate(padj = p.adjust(pValue, method = 'fdr', n = unique(freq_wide$celltype_r2)),
  #        significant = ifelse(pValue < 0.05 & effsize > 0.5, 'C1',
  #                             ifelse(pValue < 0.05 & effsize < -0.5, 'C3', 'Not')),
  #        celltype_label = case_when(abs(effsize) > 0.5 ~ Celltypes,
  #                                   abs(effsize) <= 0.5 ~ '')) |> 
  mutate(
         significant = ifelse(pValue < 0.05 & effsize > 0.5, 'C1',
                              ifelse(pValue < 0.05 & effsize < -0.5, 'C3', 'Not')),
         celltype_label = case_when(abs(effsize) > 0.5 ~ Celltypes,
                                    abs(effsize) <= 0.5 ~ '')) |>
  mutate(significant = factor(significant, levels = c('C3','Not','C1'))) |> 
  arrange(desc(effsize)) |> 
  ggplot(aes(x = effsize, y = -log10(pValue))) +
  geom_point(aes(size = abs(effsize), color = significant)) + 
  scale_size_continuous(range = c(0.5, 5), breaks = c(0.2, 0.5, 0.8, 1.2), labels = c("0.2", "0.5", "0.8", "1.2"), name = 'Effect Size') +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", size = 0.2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.2) +
  geom_text_repel(aes(label = celltype_label), box.padding = 0.6, max.overlaps = Inf, size = 2.6) +
  scale_color_manual(values = c('#2CA02CFF','gray','#1F77B4FF'), guide = FALSE) +
  theme_classic() + 
  # theme(panel.grid.major = element_line(color = "gray", size = 0.1),
  #       panel.grid.minor = element_line(color = "lightgray", size = 0.1)) +
  theme(panel.grid.major = element_line(color = "lightgray", size = 0.1)) +
  xlab("Effect Size \n(Cohen's d)")
y_text = 0.03
x1 = 0.15
x2 = 0.8
grid.segments(
  x = unit(x1 + 0.05, "npc"),
  x1 = unit(x1 - 0.05, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x2 - 0.05, "npc"),
  x1 = unit(x2 + 0.05, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("C3\n(NR-enriched)", x = unit(x1, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
grid.text("C1\n(RE-enriched)", x = unit(x2, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
dev.off()

subtypes <- c('Plasma_cell','CD4_Naive')
freq_wide |> 
  filter(celltype_r2 %in% subtypes) |>
  distinct(patient, .keep_all = T) |> 
  mutate(cluster = ifelse(cluster == 1, 'C1 ', 'C3'), 
         celltype_r2 = factor(celltype_r2, levels = subtypes)) |> 
  ggplot(aes(x = cluster, y = diff)) +
  geom_violin(aes(fill = cluster), alpha = 0.8) +
  geom_boxplot(width = 0.3, alpha = 0.2) +
  scale_fill_manual(values = c('#2CA02CFF', '#1F77B4FF'), guide = FALSE) +
  geom_point(size = 0.5) +
  facet_wrap(.~ celltype_r2, scales = 'free_y') +
  xlab("") + ylab("∆ Relative Frequency") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 10, colour = 'black'),
        strip.text = element_text(size = 12, margin = margin(2,1,2,1)),
        panel.grid.major = element_line(color = "gray", size = 0.1)) +
  stat_compare_means(method = 't.test')
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/box_subtypes.pdf', width = 5, height = 3)
  
order_column <- do.call(rbind, list_res) |> 
  mutate(padj = p.adjust(pValue, method = 'fdr', n = length(immune_cell)),
         significant = ifelse(pValue < 0.05 & effsize > 0.2, 'C2',
                              ifelse(pValue < 0.05 & effsize < -0.2, 'C3', 'Not')),
         celltype_label = case_when(abs(effsize) > 0.5 ~ Celltypes,
                                    abs(effsize) <= 0.5 ~ '')) |> 
  arrange(desc(effsize)) |> 
  filter(celltype_label != '') |> 
  pull(celltype_label)

# Supervised 
# Univariate logistic models
freq_wide <- meta_int |>
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |>
  distinct(patient, time_point, celltype_r2, .keep_all = T) |>
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  mutate(change = log2((Post + 0.01)/(Pre + 0.01)), diff = (Post - Pre))
# freq_wide$cluster <- pt_df$cluster_diff[match(freq_wide$patient, pt_df$patient)]
# freq_wide <- filter(freq_wide, !cluster  %in% c(1, NA))
# freq_wide$cluster <- ifelse(freq_wide$cluster == 2, 1, 0)
# freq_wide$cluster <- as.factor(freq_wide$cluster)
freq_wide$res <- ifelse(freq_wide$response == 'RE', 1, 0)
freq_wide$res <- as.factor(freq_wide$res)
list_res <- lapply(immune_cell[!immune_cell %in% 'B_MT2A'], function(subtype){
  print(subtype)
  df <- freq_wide |>
    filter(celltype_r2 == subtype) |>
    distinct(patient, .keep_all = T) |> 
    mutate(diff_scale = scale(diff))
  df$int_cat <- ifelse(df$interval == '<21', 0, 1)
  model <- glm(res ~ diff_scale + int_cat + modality, data = df,family = "binomial")
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  subtype_results <- coef_table['diff_scale', , drop = FALSE]
  subtype_confint <- confint_table['diff_scale', , drop = FALSE]
  combined_results <- data.frame(
    Celltypes = subtype,
    Estimate = subtype_results[1],
    StdError = subtype_results[2],
    tValue = subtype_results[4],
    pValue = coef_table['diff_scale', 'Pr(>|z|)'],
    CI_lower = subtype_confint[1],
    CI_upper = subtype_confint[2])
  return(combined_results)
})

do.call(rbind, list_res) |>
  arrange(desc(Estimate)) |>
  ggplot(aes(x = Estimate, y = -log10(pValue))) +
  geom_point() + 
  geom_label_repel(aes(label = Celltypes))

res_logi <- do.call(rbind, list_res) |>
  arrange(desc(Estimate)) |> 
  filter(pValue < 0.1)
ggplot(data = res_logi, aes(x= factor(Celltypes, levels = rev(res_logi$Celltypes)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  coord_flip() + theme_minimal()




test_result <- pt_df |> janitor::tabyl(modality, cluster_diff)  |> chisq.test()
pt_df$modality <- factor(pt_df$modality, levels = c('Mono', 'Dual')) 
pt_df |> 
  ggplot() +
  geom_mosaic(aes(x = product(modality, cluster_diff), fill = modality)) +
  labs(x = "", y = "") +
  scale_fill_manual(values = c( "#FF7F00", "#6A3D9A"), name = 'Modality') + 
  theme_mosaic() +
  theme(axis.text.x = element_text(size = 12, colour = 'black'),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate("text", x = 1, y = 1.2, label = paste("Chi-squared Test\np-value:", round(test_result$p.value, 4)), 
           hjust = 1.05, vjust = 1.5, size = 4, colour = "black") +
  geom_mosaic_text(aes(x = product(modality, cluster_diff), label = after_stat(.wt)),size = 3.5)
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/mosaic_modality.pdf', height = 5, width = 5)

test_result <- pt_df |> janitor::tabyl(treatment, cluster_diff)  |> chisq.test()
pt_df |> 
  ggplot() +
  geom_mosaic(aes(x = product(treatment, cluster_diff), fill = treatment)) +
  labs(x = "", y = "") +
  scale_fill_manual(values = structure(names = unique(pt_df$treatment), pal_npg()(length(unique(pt_df$treatment)))), name = 'Treatment') +
  theme_mosaic() +
  theme(axis.text.x = element_text(size = 12, colour = 'black'),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate("text", x = 1, y = 1.25, label = paste("Chi-squared Test\np-value:", round(test_result$p.value, 4)), 
           hjust = 1.05, vjust = 1.5, size = 4, colour = "black") +
  geom_mosaic_text(aes(x = product(treatment, cluster_diff), label = after_stat(.wt)),size = 3.5)
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/mosaic_tx.pdf', height = 5, width = 5)

pt_df$primarysite <- ifelse(pt_df$cancertype %in% c('SKCM', 'BCC', 'SCC', 'CRC', 'PCa'), 'Skin&Others',
                            ifelse(pt_df$cancertype %in% c('TNBC', 'HER2+BC', 'ER+BC'), 'Breast', 'Head and Neck'))
test_result <- pt_df |> janitor::tabyl(primarysite, cluster_diff)  |> chisq.test()
pt_df |> 
  ggplot() +
  geom_mosaic(aes(x = product(primarysite, cluster_diff), fill = primarysite)) +
  labs(x = "", y = "") +
  scale_fill_manual(values = structure(names = unique(pt_df$primarysite), pal_lancet()(length(unique(pt_df$primarysite)))), name = 'Primary Site') +
  theme_mosaic() +
  theme(axis.text.x = element_text(size = 12, colour = 'black'),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate("text", x = 1, y = 1.25, label = paste("Chi-squared Test\np-value:", round(test_result$p.value, 4)), 
           hjust = 1.05, vjust = 1.5, size = 4, colour = "black") +
  geom_mosaic_text(aes(x = product(primarysite, cluster_diff), label = after_stat(.wt)),size = 3.5)
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/mosaic_primsite.pdf', height = 5, width = 5)





# filtering <- F
# freq_filt <- meta_int |> 
#   select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
#   distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
#   pivot_wider(values_from = count_r2, names_from = time_point, values_fill = 0)
# if (filtering == T){
#   freq_filt <- freq_filt |> filter(abs(Pre-Post) >= 3, (Pre >= 3 | Post >= 3))
# }
# freq_filt$pt_r2 <- paste0(freq_filt$patient, '_', freq_filt$celltype_r2)
freq_wide <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component, int_cat) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(change = log2((Post + 0.01)/(Pre + 0.01)), diff = (Post - Pre))
freq_wide$pt_r2 <- paste0(freq_wide$patient, '_', freq_wide$celltype_r2)
# freq_wide <- filter(freq_wide, pt_r2 %in% freq_filt$pt_r2)
freq_wide$cluster <- pt_df$cluster_diff[match(freq_wide$patient, pt_df$patient)]
df_change <- freq_wide |> 
  drop_na(cluster) |> 
  filter(celltype_r2 %in% immune_cell) |> 
  group_by(cluster, celltype_r2) |> 
  summarize(median_diff = median(diff), .groups = 'drop') |> 
  pivot_wider(names_from = cluster, values_from = median_diff) |> 
  column_to_rownames(var = 'celltype_r2') 

re_celltype <- rownames(df_change)[apply(df_change, 1, function(row) which.max(row) == 1)]
re_up <- intersect(rownames(df_change)[df_change$`1` > 0], re_celltype)
re_down <- intersect(rownames(df_change)[df_change$`1` < 0], re_celltype)
nr_celltype <- rownames(df_change)[apply(df_change, 1, function(row) which.max(row) == 3)]
nr_up <- intersect(rownames(df_change)[df_change$`3` > 0], nr_celltype)
nr_down <- intersect(rownames(df_change)[df_change$`3` < 0], nr_celltype)
col_row <- ifelse(rownames(df_change) %in% c(re_up, nr_up), '#A2001F', 
                  ifelse(rownames(df_change) %in% c(re_down, nr_down), '#005D89', 'black'))
ht_change <- Heatmap(t(scale(t(df_change))), 
                     # column_title = 'Change',
                     cluster_columns = F, 
                     clustering_method_rows = "ward.D2",
                     clustering_distance_rows = "manhattan",
                     show_row_dend = F, 
                     column_names_rot = 0,
                     column_names_centered = T,
                     name = 'Median\n(z-score)', col = rev(pal(100)), 
                     rect_gp = gpar(col = "white", lwd = 1),
                     row_names_gp = gpar(fontsize = 9, col = col_row), 
                     heatmap_legend_param = list(
                       legend_direction = "horizontal", 
                       legend_width = unit(2, "cm"), at = c(-1, 0, 1),
                       title_position = "topcenter"),
                     width = ncol(df_change)*unit(5, "mm"), 
                     height = nrow(df_change)*unit(3.5, "mm"))
pdf('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/ht_celltype_change.pdf', height = 11, width = 3)
draw(ht_change, heatmap_legend_side = 'bottom')
dev.off()

df_change <- freq_wide |> 
  filter(celltype_r2 %in% immune_cell) |> 
  group_by(response, celltype_r2) |> 
  summarize(median_diff = median(diff), .groups = 'drop') |> 
  pivot_wider(names_from = response, values_from = median_diff) |> 
  column_to_rownames(var = 'celltype_r2') 
re_celltype <- rownames(df_change)[apply(df_change, 1, function(row) which.max(row) == 2)]
re_up <- intersect(rownames(df_change)[df_change$RE > 0], re_celltype)
re_down <- intersect(rownames(df_change)[df_change$RE < 0], re_celltype)
nr_celltype <- rownames(df_change)[apply(df_change, 1, function(row) which.max(row) == 1)]
nr_up <- intersect(rownames(df_change)[df_change$NR > 0], nr_celltype)
nr_down <- intersect(rownames(df_change)[df_change$NR < 0], nr_celltype)
col_row <- ifelse(rownames(df_change) %in% c(re_up, nr_up), '#A2001F', 
                  ifelse(rownames(df_change) %in% c(re_down, nr_down), '#005D89', 'black'))
ht_change <- Heatmap(t(scale(t(df_change))), 
                     # column_title = 'Change',
                     cluster_columns = F, 
                     clustering_method_rows = "ward.D2",
                     clustering_distance_rows = "manhattan",
                     show_row_dend = F, 
                     column_names_rot = 0,
                     column_names_centered = T,
                     name = 'Median\n(z-score)', col = rev(pal(100)), 
                     rect_gp = gpar(col = "white", lwd = 1),
                     row_names_gp = gpar(fontsize = 9, col = col_row), 
                     heatmap_legend_param = list(
                       legend_direction = "horizontal", 
                       legend_width = unit(2, "cm"), at = c(-1, 0, 1),
                       title_position = "topcenter"),
                     width = ncol(df_change)*unit(5, "mm"), 
                     height = nrow(df_change)*unit(3.5, "mm"))




celltype_re <- rownames(df_change)[apply(df_change, 1, function(row) which.max(row) == 2)]
celltype_nr <- rownames(df_change)[apply(df_change, 1, function(row) which.max(row) == 3)]
write.csv(celltype_re, '/bigdata/zlin/Melanoma_meta/tables/celltype_re.csv', row.names = F)
write.csv(celltype_nr, '/bigdata/zlin/Melanoma_meta/tables/celltype_nr.csv', row.names = F)
# list_res <- lapply(celltype_re, function(subtype){
#   df <- freq_wide |> 
#     filter(celltype_r2 == subtype) |> 
#     group_by(dataset) |> 
#     mutate(diff_scale = scale(diff)) |> 
#     ungroup()
#   df$res <- ifelse(df$response == 'NR', 0, 1)
#   df$int_cat <- ifelse(df$int_cat == '< 21d', 0, 1)
#   formula <- as.formula("res ~ diff + dataset + int_cat + modality")
#   model <- glm(formula, data = df, family = 'binomial')
#   model_summary <- summary(model)
#   coef_table <- coef(summary(model))
#   confint_table <- confint(model, level = 0.95)
#   subtype_results <- coef_table['diff', , drop = FALSE]
#   subtype_confint <- confint_table['diff', , drop = FALSE]
#   combined_results <- data.frame(
#     Celltypes = subtype,
#     Estimate = subtype_results[1],
#     StdError = subtype_results[2],
#     tValue = subtype_results[4],
#     pValue = coef_table['diff', 'Pr(>|z|)'],
#     CI_lower = subtype_confint[1],
#     CI_upper = subtype_confint[2]
#   )
#   return(combined_results)
# })
# 
# res <- do.call(rbind, list_res) |> as.data.frame() 
# res |> filter(pValue<0.05)








# df_pre <- freq_wide |> 
#   drop_na(cluster) |> 
#   filter(celltype_r2 %in% immune_cell) |> 
#   group_by(cluster, celltype_r2) |> 
#   summarize(median_pre = median(Pre), .groups = 'drop') |> 
#   pivot_wider(names_from = cluster, values_from = median_pre) |> 
#   column_to_rownames(var = 'celltype_r2')
# ht_pre <- Heatmap(t(scale(t(df_pre[row_order(ht_change),]))),
#                   column_title = 'Pre-Tx',
#                   cluster_columns = F, cluster_rows = F,
#                   column_names_rot = 0,
#                   show_row_names = F,
#                   name = 'Z-score (Median)', col = rev(pal(100)), 
#                   rect_gp = gpar(col = "white", lwd = 1),
#                   show_heatmap_legend = F,
#                   width = ncol(df_change)*unit(4, "mm"), 
#                   height = nrow(df_change)*unit(4, "mm"))
# pdf('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/ht_pre.pdf', height = 12, width = 3)
# ht_pre
# dev.off()

celltype <- rownames(df_change)[apply(df_change, 1, function(row) which.max(row) == 1)]
freq_wide |> 
  drop_na(cluster) |> 
  filter(celltype_r2 %in% celltype) |> 
  select(Pre, diff, cluster, celltype_r2, patient) |> 
  group_by(celltype_r2) |> 
  mutate(across(c(Pre, diff), ~scale(.) %>% as.vector)) |> 
  group_by(cluster, patient) |> 
  summarise(Pre = sum(Pre, na.rm = TRUE),
            Change = sum(diff, na.rm = TRUE)) |> 
  pivot_longer(cols = -c('cluster','patient'), names_to = 'time', values_to = 'Enrichment') |> 
  mutate(time = fct_relevel(time, "Pre", "Change"),
         cluster = factor(cluster)) |> 
  ggplot(aes(x = cluster, y = Enrichment)) +
  geom_violin(aes(fill = cluster)) + 
  geom_boxplot(width = 0.2, alpha = 0.2) +
  scale_fill_d3(name = "Cluster") +
  facet_wrap(.~ time) + 
  stat_compare_means() +
  xlab("Cluster \n(dynamic changes)") + ylab("Enrichment (z-score)")+
  theme_classic2() + 
  theme(axis.text.x = element_text(size = 12, colour = 'black'),
                           legend.position = 'none')
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/enrich_pre_change_C1.pdf', height = 3.5, width = 6)

celltype <- rownames(df_change)[apply(df_change, 1, function(row) which.max(row) == 3)]
freq_wide |> 
  drop_na(cluster) |> 
  filter(celltype_r2 %in% celltype) |> 
  select(Pre, diff, cluster, celltype_r2, patient) |> 
  group_by(celltype_r2) |> 
  mutate(across(c(Pre, diff), ~scale(.) %>% as.vector)) |> 
  group_by(cluster, patient) |> 
  summarise(Pre = sum(Pre, na.rm = TRUE),
            Change = sum(diff, na.rm = TRUE)) |> 
  pivot_longer(cols = -c('cluster','patient'), names_to = 'time', values_to = 'Enrichment') |> 
  mutate(time = fct_relevel(time, "Pre", "Change"),
         cluster = factor(cluster)) |> 
  ggplot(aes(x = cluster, y = Enrichment)) +
  geom_violin(aes(fill = cluster)) + 
  geom_boxplot(width = 0.2, alpha = 0.2) +
  scale_fill_d3(name = "Cluster") +
  facet_wrap(.~ time) + 
  stat_compare_means() +
  xlab("Cluster \n(dynamic changes)") + ylab("Enrichment (z-score)")+
  theme_classic2() + 
  theme(axis.text.x = element_text(size = 12, colour = 'black'),
        legend.position = 'none')
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/enrich_pre_change_C3.pdf', height = 3.5, width = 6)

cor_matrix <- cor_mat(meta_int, pt_df, timepoint = 'Pre', filtering = F) # c('Pre', 'diff', 'Post')
n <- 3
c_method <- 'ward.D2' # 'complete' (default) 'ward.D2'
ht_cor_pt(cor_matrix, n.cluster=n, clust_method = c_method)
pt_df$cluster_pre <- cutree(hclust(dist(cor_matrix), method = c_method), k=n)
pt_df$cluster_pre <- factor(pt_df$cluster_pre)
pt_df |> janitor::tabyl(response, cluster_pre)
# test_result <- pt_df |> janitor::tabyl(response, cluster_diff)  |> chisq.test()
pt_df$response <- factor(pt_df$response, levels = c('RE','NR'))
pt_df |> 
  ggplot() +
  geom_mosaic(aes(x = product(response, cluster_pre), fill = response)) +
  labs(x = "", y = "") +
  scale_fill_brewer("Response", palette = "Set1") + theme_mosaic() +
  theme(axis.text.x = element_text(size = 14, colour = 'black'),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  # annotate("text", x = 1, y = 1.2, label = paste("Chi-squared Test\np-value:", round(test_result$p.value, 4)), # "Chi-square Test\np-value:"
  #          hjust = 1.05, vjust = 1.5, size = 4, colour = "black") +
  geom_mosaic_text(aes(x = product(response, cluster_pre), label = after_stat(.wt)), size = 3.5) 
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/mosaic_pre.pdf', height = 5, width = 5)

pdf('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/ht_pre.pdf', height = 10, width = 20)
ht_cor_pt(cor_matrix, n.cluster=n, clust_method = c_method)
dev.off()

cor_matrix <- cor_mat(meta_int, pt_df, timepoint = 'Post', filtering = F) # c('Pre', 'diff', 'Post')
n <- 3
c_method <- 'ward.D2' # 'complete' (default) 'ward.D2'
ht_cor_pt(cor_matrix, n.cluster=n, clust_method = c_method)
pt_df$cluster_post <- cutree(hclust(dist(cor_matrix), method = c_method), k=n)
pt_df$cluster_post <- factor(pt_df$cluster_post)
pt_df |> janitor::tabyl(response, cluster_post)
# test_result <- pt_df |> janitor::tabyl(response, cluster_diff)  |> chisq.test()
pt_df$response <- factor(pt_df$response, levels = c('RE','NR'))
pt_df |> 
  ggplot() +
  geom_mosaic(aes(x = product(response, cluster_post), fill = response)) +
  labs(x = "", y = "") +
  scale_fill_brewer("Response", palette = "Set1") + theme_mosaic() +
  theme(axis.text.x = element_text(size = 14, colour = 'black'),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()) +
  # annotate("text", x = 1, y = 1.2, label = paste("Chi-squared Test\np-value:", round(test_result$p.value, 4)), # "Chi-square Test\np-value:"
  #          hjust = 1.05, vjust = 1.5, size = 4, colour = "black") +
  geom_mosaic_text(aes(x = product(response, cluster_post), label = after_stat(.wt)), size = 3.5) 
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/mosaic_post.pdf', height = 5, width = 5)

pdf('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/ht_post.pdf', height = 10, width = 20)
ht_cor_pt(cor_matrix, n.cluster=n, clust_method = c_method)
dev.off()



pt_df$response <- factor(pt_df$response, levels = c('RE','NR')) 
pt_df$dataset <- factor(pt_df$dataset, levels = c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao', 'TNBC_Zhang', 'BCC/SCC_Yost', 'HNSC_IMCISION', 'HNSC_Luoma', 'CRC_Li', 'PCa_Hawley'))
pt_df$cancertype <- factor(pt_df$cancertype, levels = c("SKCM", "BCC", "SCC", "TNBC", "ER+BC", "HER2+BC", "NSCLC", "HNSC", "PCa", "CRC"))
pt_df$modality <- factor(pt_df$modality, levels = c('Mono', 'Dual'))
pt_df$int_cat <- ifelse(pt_df$interval < 21, '< 21d', '>= 21d')
pt_df|> 
  mutate(row_order = rev(row_number())) |>
  ggplot(aes(y = row_order, axis1 = cluster_pre, axis2 = cluster_diff, axis3 = cluster_post)) +
  geom_alluvium(aes(fill = response), width = 1/12) +
  geom_stratum(width = 1/12, alpha = .25) +
  geom_label(aes(label = after_stat(stratum)), stat = 'stratum', alpha = 0.3) +
  scale_x_discrete(limits = c("Pre-Tx", "Change", "Post_Tx"), expand = c(.05, .05)) +
  scale_fill_brewer(type = 'qual', palette = 'Set1') + 
  theme_minimal() + ggtitle("") + ylab("") + xlab("") + labs(fill = 'Response') +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"))
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/alluvium_pre_change_post.pdf', width = 6, height = 4)

pt_df$response <- factor(pt_df$response, levels = c('RE','NR')) 
pt_df$dataset <- factor(pt_df$dataset, levels = c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao', 'TNBC_Zhang', 'BCC/SCC_Yost', 'HNSC_Franken', 'HNSC_IMCISION', 'HNSC_Luoma', 'CRC_Li', 'PCa_Hawley'))
pt_df$cancertype <- factor(pt_df$cancertype, levels = c("SKCM", "BCC", "SCC", "TNBC", "ER+BC", "HER2+BC", "NSCLC", "HNSC", "PCa", "CRC"))
pt_df$modality <- factor(pt_df$modality, levels = c('Mono', 'Dual'))
pt_df$int_cat <- ifelse(pt_df$interval < 21, '< 21d', '>= 21d')
pt_df |> 
  mutate(row_order = rev(row_number())) |>
  ggplot(aes(y = row_order, axis1 = dataset, axis2 = cancertype, axis3 = modality, axis4 = cluster_pre, axis5 = int_cat, axis6 = cluster_diff, axis7 = cluster_post)) +
  geom_alluvium(aes(fill = response), width = 1/12) +
  geom_stratum(width = 1/12, alpha = .25) +
  geom_label(aes(label = after_stat(stratum)), stat = 'stratum', alpha = 0.3) +
  scale_x_discrete(limits = c("Cohort", "Cancer Type", "Modality", "Cluster_pre", "Interval","Cluster_Change", "Cluster_post"), expand = c(.05, .05)) +
  scale_fill_brewer(type = 'qual', palette = 'Set1') + 
  theme_minimal() + ggtitle("") + ylab("") + xlab("") + labs(fill = 'Response') +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black")) 
ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/alluvium_comprehensive.pdf', width = 16, height = 8)

# df <- pt_df |> rename(C_pre = cluster_pre, C_change = cluster_change)|> make_long(C_pre, C_change)
# ggplot(df, aes(x = x,                        
#                next_x = next_x,                                     
#                node = node,
#                next_node = next_node,        
#                fill = factor(node))) +
#   geom_sankey(flow.alpha = 0.3,               #This Creates the transparency of your node 
#               node.color = "black",   # This is your node color        
#               show.legend = TRUE)  +   # This determines if you want your legend to show
#   scale_fill_d3() + xlab("") + 
#   theme_sankey(base_size = 12) + guides(fill = guide_legend(title = "Cluster"))
# ggsave('/bigdata/zlin/Melanoma_meta/figures/Cluster_response/sankey_pre_change.pdf')


