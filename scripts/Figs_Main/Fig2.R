pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','BPCells','SeuratExtend','MetBrewer','ggplot2','CytoTRACE2','AnnotationDbi','org.Hs.eg.db','effsize','metafor','rstatix','ggpubr','rcompanion')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(max.print = 10000)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024) 

# CD4+T cells
seu <- qs_read('data/seu_CD4T.qs2')
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
seu$celltype_r2 <- factor(seu$celltype_r2, levels =  c('CD4_T-naive','CD4_Tcm','CD4_Tstr','CD4_Tctl','CD4_Tfh','CD4_Th17','CD4_Treg','CD4_T-ISG','Cycling'))
p <- DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', 
             cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=1, override.aes = list(size = 3)))
p + theme(legend.position = "none") 
ggsave('figures/UMAP/UMAP_CD4T.png', height = 4, width = 4, dpi = 300)
get_legend(p) |> as_ggplot()
ggsave('figures/UMAP/legend_CD4T.pdf')

# CD8+T/NK cells
seu <- qs_read('data/seu_CD8TNK.qs2')
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
seu$celltype_r2 <- factor(seu$celltype_r2, levels = c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
                                                      "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
                                                      "CD8_NK-like",'MAIT','gdT','NK_CD56loCD16hi','NK_CD56hiCD16lo','Cycling'))
p <- DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', 
             cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
p + theme(legend.position = "none") 
ggsave('figures/UMAP/UMAP_CD8T.png', height = 4, width = 4, dpi = 300)
get_legend(p) |> as_ggplot()
ggsave('figures/UMAP/legend_CD8T.pdf')

# Myeloid cells
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
seu$celltype_r2[seu$celltype_r2 %in% c('cDC2_CD1C','cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'MoDC')] <- 'cDC2'
myeloids <- c('Mast','pDC','cDC1', 'cDC2', 'mregDC',
              'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
              'Macro_IL1B', 'Macro_INHBA', 'Macro_FN1', 'Macro_SPP1', 'Macro-ISG', 
              'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2', 'Cycling')
seu$celltype_r2 <- factor(seu$celltype_r2, levels = myeloids)
p <- DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
             cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
p + theme(legend.position = "none") 
ggsave('figures/UMAP/UMAP_myeloids.png', height = 4, width = 4, dpi = 300)
get_legend(p) |> as_ggplot()
ggsave('figures/UMAP/legend_myeloids.pdf')

# B/plasma cells
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
seu$celltype_r2[seu$celltype_r2 == 'PC-early_LTB'] <- "PC_IGHG"
seu$celltype_r2 <- factor(seu$celltype_r2, levels =  c("B-naive", "B-ISG", "B-HSP", "B_MT2A", 
                                                       "ACB_EGR1", "ACB_NR4A2", "ACB_CCR7", "B-memory", "B-AtM", 
                                                       "GCB-pre", "GCB-DZ_SUGCT", "GCB-LZ_LMO2",
                                                       "GCB-cycling", "PC-cycling","PC-trans",
                                                       "PC-early_RGS13", "PC_IGHG", "PC_IGHA"))
p <- DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
             cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
p + theme(legend.position = "none") 
ggsave('figures/UMAP/UMAP_Bplasma.png', height = 4, width = 4, dpi = 300)
get_legend(p) |> as_ggplot()
ggsave('figures/UMAP/legend_Bplasma.pdf')

# Non-immune cells
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
non_immune <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein","Pericytes","SMC", "Myofibroblasts","CAF-desmo",  "CAF_SFRP2", 
                "CAF-prog", "iCAF_MMP1", "iCAF_IL6", "CAF-ap","Cycling")
seu$celltype_r2 <- factor(seu$celltype_r2, levels = non_immune)
p <- DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
             cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
p + theme(legend.position = "none") 
ggsave('figures/UMAP/UMAP_nonimmune.png', height = 4, width = 4, dpi = 300)

# Compositional Analysis 
source('scripts/Celltype_classification.R')
metadata <- read.csv('tables/meta_all.csv')
unwanted_celltypes <- c('Melanocytes(CNA-)', 'Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Epithelial(CNA-)', 'Malignant(CNA+)')
celltypes <- metadata$celltype_r2 |> unique() |> setdiff(unwanted_celltypes)
metadata <- metadata[-which(metadata$cohort %in% c('NSCLC_Liu')),]
metadata$cohort[metadata$cohort == 'BCC&SCC_Yost'] <- 'BCC_Yost'

df <- metadata |> 
  filter(!celltype_r2 %in% unwanted_celltypes, 
         count_immune >= 50, non_malignant_count >= 100, 
         subset %in% c("CD45+sorted", 'TME'),
         tx_status %in% c("Baseline", "Treated")) |> 
  select(patient, sample, celltype_r2, freq_r2_comp, subtype, modality) |>
  distinct(sample, celltype_r2, freq_r2_comp, subtype, modality, .keep_all = T) 

df$proportion <- df$freq_r2_comp
df$response <- metadata$response[match(df$sample, metadata$sample)]
df$tx_status <- metadata$tx_status[match(df$sample, metadata$sample)]
df$cohort <- metadata$cohort[match(df$sample, metadata$sample)]
df$treatment <- metadata$treatment[match(df$sample, metadata$sample)]

# Per-cohort PAIRED analysis
# patient with paired samples
pt <- df |>
  distinct(sample, .keep_all = TRUE) |>
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |>
  pull(patient) 
# patients in Becker without matched tumors
unmatched_pt <- c("SKCM_this study_Patient2", "SKCM_this study_Patient3")
set.seed(1234)
wilcox_per_cohort_results_list <- list()
for (cohort_name in unique(df$cohort)) {
  print(cohort_name)
  for (celltype in setdiff(unique(df$celltype_r2), c('pDC','Mast'))) {
    sub_df <- df |> filter(cohort == cohort_name, celltype_r2 == celltype, patient %in% setdiff(pt, unmatched_pt))
    
    # Check for paired data structure
    if(length(unique(sub_df$tx_status)) < 2 || nrow(sub_df) < 2) next
    
    tryCatch({
      # Pivot to wide format for paired analysis WITH zero-filling
      wide_df <- sub_df |>
        select(patient, tx_status, freq_r2_comp) |>
        pivot_wider(names_from = tx_status, values_from = freq_r2_comp, values_fill = 0) 
      
      # Ensure both time points are present (create if missing)
      if(!"Baseline" %in% colnames(wide_df)) wide_df$Baseline <- 0
      if(!"Treated" %in% colnames(wide_df)) wide_df$Treated <- 0
      
      # Check if we have enough subjects with data
      if(nrow(wide_df) < 2) next
      
      # Perform PAIRED Wilcoxon test
      wilcox_res <- wilcox.test(wide_df$Treated, wide_df$Baseline, 
                                paired = TRUE, exact = FALSE, conf.int = TRUE)
      
      # Calculate effect size for paired data - Matched-pairs rank biserial
      n_pairs <- nrow(wide_df)
      differences <- wide_df$Treated - wide_df$Baseline
      
      # Matched pairs rank biserial correlation
      effect_size <- wilcox_res$statistic / (n_pairs * (n_pairs + 1) / 2) - 0.5
      effect_size <- effect_size * 2  # Scale to range [-1, 1]
      
      # For SE calculation - use the CI method for better accuracy
      # Since we have paired data with potential zeros, bootstrap might be better
      # But for simplicity, we'll use the correlation approximation
      effect_se <- sqrt((1 - effect_size^2) / (n_pairs - 2))
      
      wilcox_per_cohort_results_list[[length(wilcox_per_cohort_results_list) + 1]] <- tibble(
        cohort = cohort_name,
        celltype_r2 = celltype,
        p_value = wilcox_res$p.value,
        effect_size = effect_size,
        effect_se = effect_se,
        n_pairs = n_pairs,
        median_baseline = median(wide_df$Baseline),
        median_treated = median(wide_df$Treated),
        mean_diff = mean(differences),
        sd_diff = sd(differences),
        # Add min/max for data quality assessment
        min_baseline = min(wide_df$Baseline),
        max_baseline = max(wide_df$Baseline),
        min_treated = min(wide_df$Treated),
        max_treated = max(wide_df$Treated)
      )
    }, error = function(e) {
      message(paste("Error processing", celltype, "in cohort", cohort_name, ":", e$message))
    })
  }
}
final_cohort_results <- bind_rows(wilcox_per_cohort_results_list)

# Meta-analysis 
wilcox_meta_results_list <- list()
for (celltype in unique(final_cohort_results$celltype_r2)) {
  print(celltype)
  cell_data <- final_cohort_results |> filter(celltype_r2 == celltype)
  
  if (nrow(cell_data) >= 2) {
    meta_res <- tryCatch({
      rma(
        yi = effect_size,
        sei = effect_se,
        data = cell_data,
        method = "REML"
      )
    }, error = function(e) {
      message(paste("Meta-analysis failed for", celltype, ":", e$message))
      NULL
    })
    
    if (!is.null(meta_res)) {
      wilcox_meta_results_list[[celltype]] <- tibble(
        celltype_r2 = celltype,
        pooled_effect = meta_res$b[1],
        pooled_effect_se = meta_res$se,
        ci_lower = meta_res$ci.lb,
        ci_upper = meta_res$ci.ub,
        p_value = meta_res$pval,
        i2 = meta_res$I2,
        tau2 = meta_res$tau2,
        n_cohorts = nrow(cell_data),
        total_pairs = sum(cell_data$n_pairs)  
      )
    }
  }
}

final_results <- bind_rows(wilcox_meta_results_list) |> arrange(desc(pooled_effect))
write.csv(final_results |> arrange(desc(pooled_effect)) |> data.frame(), 'tables/wicox_paired_r2_per_cohort_tx_status_meta.csv', row.names = F)
# Significant celltype
celltypes_filtered <- final_results |> arrange(desc(pooled_effect)) |> filter(p_value < 0.1) |> pull(celltype_r2) |> unique()
library(grid)
pdf('Dynamics/fp_tx_status_paired.pdf', height = 6, width = 4.5)
p <- final_results |> arrange(desc(pooled_effect)) |> filter(p_value < 0.1) |> 
  mutate(effect_direction = ifelse(pooled_effect > 0, "Treated", "Baseline")) |> 
  ggplot(aes(x = pooled_effect, y = reorder(celltype_r2, pooled_effect), color = effect_direction)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, color = 'black') +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Baseline" = "#04a3bd", "Treated" = "#f0be3d")) +
  # Customize the plot labels and title
  labs(
    x = "Effect Size \n(pooled)",
    y = "Cell Type",
    color = "Effect Direction" # Rename the legend title
  ) +
  # Use a clean theme and adjust text size for better readability
  theme_minimal() +
  theme(plot.margin = unit(c(0.5, 2.5, 1.5, 0.5), "lines"),
        axis.title.x = element_text(size = 18, color = 'black', vjust = -5),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.text.x = element_text(size = 10, color = 'black'),
        axis.title.y = element_blank(), 
        legend.position = "none" 
  ) +
  coord_cartesian(xlim = c(-0.7, 0.7), clip = "off") 
grid.draw(p)
y_text = 0.07
x1 = 0.46
x2 = 0.745
grid.segments(
  x = unit(x1 + 0.07, "npc"),
  x1 = unit(x1 - 0.1, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x2 - 0.07, "npc"),
  x1 = unit(x2 + 0.1, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Baseline", x = unit(x1, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 11))
grid.text("Treated", x = unit(x2, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 11))
dev.off()
# dotplot
annotation_data <- data.frame(
  cohort = unique(final_cohort_results$cohort),
  treatment = metadata$treatment[match(unique(final_cohort_results$cohort), metadata$cohort)]
)

col_order <- c("SKCM_this study", "SKCM_Plozniak", 
               "BCC_Yost", "BRCA_Bassez1", "BRCA_Bassez2", 
               "TNBC_Bassez1", "TNBC_Bassez2", "TNBC_Zhang", "TNBC_Shiao", 
               "HNSC_Franken", "HNSC_vanderLeun", "HNSC_Luoma(Combo)", "HNSC_Luoma", 
               "CRC_Li", "CRC_Chen", 
               "NSCLC_Yan", "PCa_Hawley", "RCC_Bi", "HCC_Guo", "HCC_Ma")
annotation_data$cohort <- factor(annotation_data$cohort, levels = col_order)

annotation_plot <- ggplot(annotation_data, aes(x = cohort, y = 1)) +
  geom_tile(aes(fill = treatment), color = "white", size = 0.5) +
  scale_fill_manual(values = met.brewer('Austria', n = length(unique(annotation_data$treatment)))) +
  theme_void() + # Remove all plot elements
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = unit(c(0.5, 0, 0, 0), "lines"),
    legend.key.size = unit(0.6, "cm")) +
  labs(fill = "Treatment")
dotplot <- final_cohort_results |> 
  filter(celltype_r2 %in% celltypes_filtered) |>
  ggplot(aes(x = factor(cohort, levels = col_order), y = factor(celltype_r2, levels = rev(celltypes_filtered)))) +
  geom_point(aes(color = pmax(-0.7, pmin(0.7, effect_size)), 
                 size = pmax(1, pmin(6, 1 / ifelse(is.na(effect_se), 1e+10, effect_se))))) + 
  geom_point(data = final_cohort_results |> filter(p_value < 0.05, celltype_r2 %in% celltypes_filtered),
             aes(size = pmax(1, pmin(6, 1 / ifelse(is.na(effect_se), 1e+10, effect_se)))), 
             shape = 1,
             color = "black",
             fill = NA,
             stroke = 1) +
  scale_color_gradient2(low = "#04a3bd", mid = "white", high = "#f0be3d", midpoint = 0, limits = c(-0.7, 0.7), na.value = "lightgray") +
  scale_size_continuous(range = c(0, 5), 
                        limits = c(0, 6),
                        breaks = c(2, 4, 6),  # Changed to 2, 4, 6
                        labels = c("2", "4", "6")) +
  labs(
    x = "",
    y = "",
    color = "Effect size",
    size = "Confidence"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 10),
        legend.position = "right",legend.key.size = unit(0.5, "cm"), # Adjust size of the keys
        legend.title = element_text(size = 10), # Adjust title text size
        legend.text = element_text(size = 10),
        plot.margin = unit(c(0.5, 0, 1.5, 0), 'lines'))
annotation_plot / dotplot +
  plot_layout(heights = c(0.05, 1), guides = 'collect')
ggsave('Dynamics/dot_tx_status_paired.pdf', height = 6, width = 6)


