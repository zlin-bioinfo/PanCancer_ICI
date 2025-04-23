rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','tibble','pheatmap','qs','MetBrewer','viridis','ComplexHeatmap','colorRamp2','corrr','ggnewscale','NMF','ggpubr','corrplot','rstatix','ggforce','ggnetwork','igraph')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
setwd("/home/zlin/workspace/PanCancer_ICI")
source('scripts/Celltype_classification.R')
metadata <- read.csv('tables/meta_all.csv') 
patient_info <- metadata |> distinct(patient, .keep_all = T)
metadata$freq_r2_comp[metadata$celltype_r2 == 'Malignant(CNA+)'] <- metadata$freq_r2[metadata$celltype_r2 == 'Malignant(CNA+)']
pt_paired <- metadata |>
  distinct(sample, .keep_all = TRUE) |>
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |>
  pull(patient)
df_diff <- metadata |> 
  filter(!celltype_r2 %in% c('Melanocytes(CNA-)','Epithelial(CNA-)','Cycling T','Cycling NK', "GCB-cycling", "PC-cycling", 'Cycling non-immune', 'Cycling myeloids','Cycling T/NK'),
         # subset %in% c('All TME', 'CD45+sorted'),
         subset %in% c('All TME'),
         # cohort != 'RCC_Bi',
         patient %in% pt_paired) |> 
  select(patient,  sample, freq_r2_comp, celltype_r2, time_point) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  select(-sample) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(diff=On-Pre) |> 
  select(!Pre:On) |> 
  pivot_wider(values_from = diff, names_from = celltype_r2, values_fill = 0) |> 
  column_to_rownames(var = 'patient')

# cor.mat <- cor_mat(df_diff)
# cor.mat.p <- cor.mat |>  cor_get_pval()
# cor.mat <- cor.mat |> column_to_rownames(var = 'rowname')
# cor.mat.p <- cor.mat.p |> column_to_rownames(var = 'rowname')
# cor_plot(cor.mat, p.mat = cor.mat.p)

cor.res <- cor(df_diff) |> as.data.frame()
cor.test.res <- cor.mtest(df_diff, conf.level = 0.95) 
cor.test.res <- cor.test.res$p

pdf('figures/Co-regulating/ht_tme.pdf', height =12, width = 13)
Heatmap(cor.res, col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), 
        cluster_rows = T, cluster_columns = T, column_names_rot = 45, row_names_side = "right", column_names_side = 'bottom',
        show_column_dend = F, show_row_dend = F,
        heatmap_legend_param = list(title = "Pearson's \nCorrelation", 
                                    at = c(-0.8, 0, 0.8), 
                                    labels = c("0.8", "0", "0.8"),
                                    legend_direction = "horizontal", 
                                    legend_width = unit(1.5, "cm"), 
                                    legend_side = 'bottom',
                                    title_position = "topcenter"),
        rect_gp = gpar(col ="grey", lwd = 0.5),
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 7),
        width = ncol(cor.res)*unit(3, "mm"),
        height = nrow(cor.res)*unit(3, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(cor.test.res[i, j])) {
            if(cor.test.res[i, j] < 0.05) {
              gb = textGrob("*", gp = gpar(fontsize=14))
              gb_h = convertHeight(grobHeight(gb), "mm")
              gb_w = convertWidth(grobWidth(gb), "mm")
              grid.text('*', x, y - gb_h*0.5 + gb_w*0.4, gp=gpar(fontsize=14))
            }
          }
        })
dev.off()

# Create a tidy data frame of correlations
tidy_cors <- df_diff |> 
  correlate() |> 
  stretch()
df_p <- cor.test.res
diag(df_p) <- NA
df <- df_p |> as.data.frame() |> 
  rownames_to_column(var = 'x') |> 
  pivot_longer(values_to = 'p', names_to = 'y', cols = -x) |> 
  mutate(FDR = p.adjust(p, method = "fdr")) |> 
  merge(tidy_cors, by = c('x','y'))
graph_cors <- df |> 
  filter(FDR<0.05, r > 0.3) |> 
  graph_from_data_frame(directed = FALSE)
V(graph_cors)$`Cell type` <- ifelse(V(graph_cors)$name %in% t_nk, 'T/NK cells',
                                    ifelse(V(graph_cors)$name %in% mye, 'Myeloid cells',
                                           ifelse(V(graph_cors)$name %in% bplasma, 'B lymphocytes', 'Non-immune')))
V(graph_cors)$`Cell type` <- factor(V(graph_cors)$`Cell type`, levels = c("T/NK cells", "B lymphocytes", "Myeloid cells","Non-immune"))
communities <- cluster_leiden(graph_cors, resolution = 0.01)
communities$nb_clusters
V(graph_cors)$Cluster <- as.factor(communities$membership)
df_cluster <- data.frame(cluster = V(graph_cors)$Cluster, celltype = V(graph_cors)$name, cellgroup = V(graph_cors)$`Cell type`)
cluster_to_mark <- df_cluster |>
  distinct(cluster, cellgroup, .keep_all = T) |>
  group_by(cluster) |>
  dplyr::summarize(count = n()) |>
  filter(count > 1) |>
  pull(cluster)
p_pos <- graph_cors |> 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = r, size = -log10(FDR))) +
  scale_size(name = 'Edge Width\n -log(FDR)') +
  scale_colour_gradientn(
    name = "Pearson Correlation", limits = c(0, 1), breaks = c(0,1), 
    colors = colorRampPalette(c("white", "#B2182B"))(50),  # Replaced #CF0034
    guide = guide_colorbar(direction = "horizontal", 
                           title.position = "top", 
                           title.hjust = 0.5)) +
  new_scale_color() + 
  geom_nodes(aes(color = `Cell type`), size = 3) +
  scale_color_manual(values=c(
    'T/NK cells' = "#1D3557",   # Deep Navy Blue (Avoiding #154999)
    'Myeloid cells' = '#2A9D8F', # Teal Green
    'B lymphocytes' = '#D4A373', # Warm Earthy Brown (Distinct from Purple)
    'Non-immune' = '#F4A261'     # Warm Orange (More contrast from Yellow)
  ), name = 'Node Color') +
  # scale_color_manual(values=met.brewer("Juarez", 6)[-c(1,4)], name = 'Node Color') +
  # geom_mark_hull(
  #   aes(x, y, group = Cluster, filter = Cluster %in% cluster_to_mark),
  #   concavity = 3,
  #   expand = unit(2.5, "mm"),
  #   color = 'darkgray') +
  geom_nodetext_repel(aes(label = name), size = 3) +
  guides(alpha = "none") +
  theme_blank() + 
  ggtitle('TME \n(FDR<0.05 & r>0.25)'); p_pos
ggsave('figures/Co-regulating/graph_TME_pos.pdf', height = 6, width = 7)

graph_cors <- df |> 
  filter(FDR<0.05, r< -0.3) |> 
  graph_from_data_frame(directed = FALSE)
V(graph_cors)$`Cell type` <- ifelse(V(graph_cors)$name %in% t_nk, 'T/NK cells',
                                    ifelse(V(graph_cors)$name %in% mye, 'Myeloid cells',
                                           ifelse(V(graph_cors)$name %in% bplasma, 'B lymphocytes', 'Non-immune')))
V(graph_cors)$`Cell type` <- factor(V(graph_cors)$`Cell type`, levels = c("T/NK cells", "B lymphocytes", "Myeloid cells","Non-immune"))
p_neg <- graph_cors |> 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = r, size = -log10(FDR))) +
  scale_size(name = 'Edge Width\n -log(p)') +
  scale_colour_gradientn(
    name = "Pearson Correlation", limits = c(-1, 0), breaks = c(-1,0,1), 
    colors = colorRampPalette(c("#336699", "white"))(50),
    guide = guide_colorbar(direction = "horizontal", 
                           title.position = "top", 
                           title.hjust = 0.5)) +
  new_scale_color() + 
  geom_nodes(aes(color = `Cell type`), size = 5) +
  scale_color_manual(values=c(
    'T/NK cells' = "#1D3557",   # Deep Navy Blue (Avoiding #154999)
    'Myeloid cells' = '#2A9D8F', # Teal Green
    'B lymphocytes' = '#D4A373', # Warm Earthy Brown (Distinct from Purple)
    'Non-immune' = '#F4A261'     # Warm Orange (More contrast from Yellow)
  ), name = 'Node Color') +
  # geom_mark_hull(
  #   aes(x, y, group = Cluster, filter = Cluster %in% cluster_to_mark),
  #   concavity = 5,
  #   expand = unit(2.5, "mm"), 
  #   color = 'darkgray') +
  geom_nodetext_repel(aes(label = name), size = 3) +
  guides(alpha = "none") +
  theme_blank() + 
  ggtitle('TME \n(FDR<0.05, r< -0.3)'); p_neg
ggsave('figures/Co-regulating/graph_TME_neg.pdf', height = 6, width = 7)

df_diff_sub <- df_diff[,c('Macro_C1QC','CD4_Treg','CAF-desmo','CD4_T-naive','CD8_T-naive','B-naive','CAF_SFRP2','Endo-lymphatic','cDC2_CD1C')]
df_diff_sub <- df_diff_sub[order(df_diff_sub$Macro_C1QC-df_diff_sub$`CD4_T-naive`),]
# score_pos <- apply(df_diff_sub[,c('Macro_C1QC','CD4_Treg','CAF-desmo')], 1, sum)
# score_neg <- apply(df_diff_sub[,c('CD4_T-naive','CD8_T-naive','B-naive','CAF_SFRP2','Endo-lymphatic','cDC2_CD1C')], 1, sum)
# df_diff_sub <- df_diff_sub[order(score_pos-score_neg),]

patient_info_sub <- patient_info |> filter(patient %in% rownames(df_diff_sub)) 
patient_info_sub <- patient_info_sub[match(rownames(df_diff_sub), patient_info_sub$patient),]
col_ha = HeatmapAnnotation(
  `Cancer Type` = patient_info_sub$subtype,
  # Modality = patient_info_sub$modality,
  # `Response Metrics` = patient_info_sub$res_metric,
  Response = patient_info_sub$response,
  col = list(`Cancer Type` = structure(names = unique(patient_info_sub$subtype), met.brewer('Juarez',length(unique(patient_info_sub$subtype)))),
             # Treatment = structure(names = unique(patient_info_sub$treatment), pal_npg()(length(unique(patient_info_sub$treatment)))),
             # Modality = c('Mono' = "#FF7F00", 'Dual' = "#6A3D9A"),
             Response = c('R' = '#CC0C00FF','NR' = '#5C88DAFF','NE' = '#84BD00FF'),
             # `Response Metrics` = structure(names = sort(unique(patient_info_sub$res_metric)), met.brewer("Juarez", length(unique(patient_info_sub$res_metric)))),
             Cluster = structure(names = unique(patient_info_sub$cluster), pal_d3()(length(unique(patient_info_sub$cluster))))
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
pdf('figures/Co-regulating/ht_cd4t_c1qc.pdf', height = 3, width = 10)
df_diff_sub |> t() |> scale() |> 
  Heatmap(show_column_names = F, cluster_columns = F, 
          col = circlize::colorRamp2(c(-2, 0, 2), c("#154999", "white", "#CF0034")),
          heatmap_legend_param = list(title = "Dynamics \n(On-Pre)", 
                                      at = c(-2, 0, 2), 
                                      labels = c("Min", "", "Max")),
          top_annotation = col_ha)
dev.off()
library(tidyplots)
my_colors <- c('R' = '#CC0C00FF','NR' = '#5C88DAFF')
patient_info_sub$group <- 'Mid'
patient_info_sub <- filter(patient_info_sub, response != 'NE')
patient_info_sub$group[round(152*3/4):152] <- 'Macro_C1QC.high'
patient_info_sub$group[1:round(152/4)] <- 'CD4_T-naive.high'
patient_info_sub |> tabyl(response,group) |> data.frame() |> select(-Mid) |> 
  pivot_longer(cols = c(Macro_C1QC.high,CD4_T.naive.high),values_to = 'count', names_to = 'Group') |> 
  tidyplot(x = Group, y = count, color = response) |>
  add_barstack_relative(reverse = F) |> 
  adjust_colors(my_colors) |>
  theme_minimal_x() |>
  adjust_size(50, 50) |>
  remove_legend_title() |> 
  remove_y_axis_title()
ggsave('figures/Co-regulating/cd4tvsc1qc.pdf', height = 3, width = 3.5)

