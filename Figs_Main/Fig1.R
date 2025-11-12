pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','tibble','qs2','MetBrewer','janitor','effsize','metafor','rstatix','ggpubr','rcompanion')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
source('scripts/Celltype_classification.R')

# UMAP Main
datasets <- c('SKCM_Plozniak', 'BCC_Yost',
              'SCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'NSCLC_Yan', 
              'NSCLC_Liu',
              'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo','HCC_Ma','RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2'))
  seu <- subset(seu, subset = celltype_main != 'Neutrophils')
  seu <- seu[,!is.na(seu$celltype_r2)]
  seu$cell.id <- colnames(seu)
  seu <- seu |>
    subset(subset = sample %in% filter_sample) |>
    NormalizeData() |>
    FindVariableFeatures()
  return(seu)
})
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)]) |>
  SketchData(ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
seu <- seu |> 
  FindVariableFeatures() |> 
  ScaleData() |>
  RunPCA(verbose=T) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony', verbose = T) |> 
  FindNeighbors(reduction = "harmony", dims = 1:30) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:30, reduction = 'harmony')
seu <- ProjectIntegration(seu, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
seu <- ProjectData(seu, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                   full.reduction = "harmony.full", dims = 1:30, refdata = list(seurat_clusters_full = 'seurat_clusters'))
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",
               reduction.key = "UMAP_full_")
seu$celltype_main[seu$celltype_main %in% c('Melanocytes(CNA+)','Melanocytes(CNA-)')] <- 'Melanocytes'
seu$celltype_main[seu$celltype_main %in% c('Epithelial(CNA+)','Epithelial(CNA-)')] <- 'Epithelial'
seu$celltype_main[seu$celltype_main == 'Malignant(CNA+)'] <- 'Epithelial'
qs_save(seu, 'data/seu_full_integrated.qs2')
seu <- qs_read('data/seu_full_integrated.qs2')

seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% t_nk] <- 'Cycling T/NK'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c("B-naive", "B-ISG", "B-HSP", "B_MT2A", 
                                                                          "ACB_EGR1", "ACB_NR4A2", "ACB_CCR7", "B-memory", "B-AtM", 
                                                                          "GCB-pre", "GCB-DZ_SUGCT", "GCB-LZ_LMO2",
                                                                          "GCB-cycling")] <- 'B'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c("PC-cycling","PC-trans","PC-early_RGS13", "PC_IGHG", "PC_IGHA")] <- 'Plasma'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein")] <- 'Endo'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c("Pericytes","SMC")] <- 'Mural'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c("Myofibroblasts", "CAF_SFRP2","CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap")] <- 'CAF'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c('Mast','pDC','cDC1',
                                                                          'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'mregDC', 'MoDC',
                                                                          'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
                                                                          'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro-ISG',
                                                                          'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2')] <- 'Cycling myeloids'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c('Epithelial','Melanocytes','Malignant(CNA+)')] <- 'Epithelial/melanocytes'
seu$celltype_main[seu$celltype_main %in% c('Epithelial','Melanocytes','Malignant(CNA+)')] <- 'Epithelial/melanocytes'
seu$celltype_main <- factor(seu$celltype_main, levels = c("CD4+T", "CD8+T", "NK", "Cycling T/NK", "B", "Plasma",
                                                          "Mast", "pDC", "cDC", "Mono/macro", 'Cycling myeloids',
                                                          "Endo", "Mural", "CAF", "Epithelial/melanocytes"))
seu <- subset(seu, subset = cohort %in% c('SCC_Yost', 'NSCLC_Liu'), invert = T)
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
# seu # 1445482 cells
library(ggrastr)
p <- DimPlot(seu, group.by = 'celltype_main', reduction = 'umap.full', alpha = 0.1, 
             cols = rev(color_pro(length(unique(seu$celltype_main)), 1)), 
             label = T, label.size = 4.5, repel = F, raster = T) + 
  ggtitle('') + theme_void() + 
  theme(legend.text = element_text(size=14), 
        legen= unit(0.5, 'cm')) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1]) , y = min(seu@reductions$umap.full@cell.embeddings[,2]),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])) + 3, yend = min(seu@reductions$umap.full@cell.embeddings[,2])),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1]) , y = min(seu@reductions$umap.full@cell.embeddings[,2]),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])) , yend =min(seu@reductions$umap.full@cell.embeddings[,2]) + 3),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) + 1.5, y = min(seu@reductions$umap.full@cell.embeddings[,2]) -1, label = "UMAP_1",
           color="black", size = 3) + 
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) -1, y = min(seu@reductions$umap.full@cell.embeddings[,2]) + 1.5, label = "UMAP_2",
           color="black",size = 3, angle=90)
# rasterize(p, layers='Point', dpi=300)
ggsave('figures/UMAP/UMAP_main.png', height = 6, width = 7.5, dpi=300)
# seu$cohort[seu$cohort %in% c("BCC_Yost", "SCC_Yost")] <- "BCC&SCC_Yost"
DimPlot(seu, group.by = 'cohort', reduction = 'umap.full', alpha=0.9, 
        cols = rev(color_pro(length(unique(seu$cohort)), 1))) +
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5),
                                     legend.text=element_text(size=16)) +
  guides(color = guide_legend(ncol=1, override.aes = list(size = 4)))
ggsave('figures/UMAP/UMAP_cohort.png', height = 5, width = 7, dpi = 300)genes_to_check = list(c('CD3D','CD4','CD8A'), # T cells 'CD8B'
                                                                                              c('KLRD1','FCGR3A'), 
                                                                                              c('MKI67', 'TOP2A'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                                                                                              c('CD79A','CD19', 'MS4A1'),  # B cells 
                                                                                              c('MZB1','JCHAIN'),
                                                                                              c('KIT','TPSAB1'),
                                                                                              c('LILRA4','PLD4'),
                                                                                              c('CLEC9A','CD1C','LAMP3'), 
                                                                                              c('CD68', 'LYZ', 'CD14'),  
                                                                                              c('PECAM1','VWF', 'ENG'),
                                                                                              c("RGS5",'ACTA2'),
                                                                                              c('COL1A1','FAP'),
                                                                                              c('KRT19', 'EPCAM'),
                                                                                              c('MLANA','TYR'),
                                                                                              c('RPL11','RPL10A')
)
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis()

# Compositional Analysis (main level)
metadata <- read.csv('tables/meta_all.csv')
metadata <- metadata[-which(metadata$cohort %in% c('NSCLC_Liu')),]
metadata$cohort[metadata$cohort == 'BCC&SCC_Yost'] <- 'BCC_Yost'

unwanted_celltypes <- c('Melanocytes(CNA-)', 'Epithelial(CNA-)')
celltypes <- metadata$celltype_main |> unique() |> setdiff(unwanted_celltypes)
df <- metadata |> 
  filter(count_immune >= 50, non_malignant_count >= 100, 
         # response %in% c('R','NR'), 
         tx_status %in% c("Baseline", "Treated"),
         subset %in% c('TME')) |> 
  mutate(celltype_main = case_when(celltype_main == 'Cycling' ~ celltype_r2, .default = celltype_main)) |> 
  mutate(celltype_main = case_when(celltype_main %in% c('B', 'B-AtM', 'B-ISG', 'GCB-cycling') ~ 'B cells',
                                   # celltype_main == 'Mono/macro' ~ 'Mono/macro',
                                   celltype_main %in% c('Cycling T', "Cycling NK") ~ 'Cycling T/NK',
                                   celltype_main == 'CAF' ~ 'CAFs',
                                   celltype_main == 'CD8+T' ~ 'CD8+T cells',
                                   celltype_main %in% c('PC-cycling', 'Plasma', "PC_IGHG", "PC_IGHA", "PC-trans", "PC-early_RGS13") ~ 'Plasma cells',
                                   # celltype_main == 'pDC' ~ 'pDCs',
                                   celltype_main == 'CD4+T' ~ 'CD4+T cells',
                                   # celltype_main == 'cDC' ~ 'cDC',
                                   celltype_main == 'Endo' ~ 'Endothelial cells',
                                   celltype_main == 'Cycling non-immune' ~ 'Cycling non-immune',
                                   celltype_main == 'NK' ~ 'NK cells',
                                   celltype_main == 'Mast' ~ 'Mast cells',
                                   celltype_main == 'Mural' ~ 'Mural cells',
                                   celltype_main %in% c('Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Malignant(CNA+)') ~ 'Malignant cells',
                                   .default = celltype_main)) |> 
  filter(!celltype_main %in% unwanted_celltypes) |> 
  group_by(sample, celltype_main) |> 
  mutate(count_main = n()) |> 
  ungroup() |> 
  mutate(freq_main = count_main/count_sample) |>
  distinct(sample, celltype_main, freq_main) 

df$proportion <- df$freq_main
df$response <- metadata$response[match(df$sample, metadata$sample)]
df$tx_status <- metadata$tx_status[match(df$sample, metadata$sample)]
df$cohort <- metadata$cohort[match(df$sample, metadata$sample)]
df$treatment <- metadata$treatment[match(df$sample, metadata$sample)]
df <- df |> mutate(
  tx_status = factor(tx_status, levels = c("Baseline", "Treated"))
)
celltypes <- unique(df$celltype_main) 
# Per-cohort
wilcox_per_cohort_results_list <- list()
for (cohort_name in unique(df$cohort)) {
  print(cohort_name)
  for (celltype in unique(df$celltype_main)) {
    sub_df <- df |> filter(cohort == cohort_name, celltype_main == celltype)
    if(length(unique(sub_df$tx_status)) < 2 || nrow(sub_df) < 2) next
    tryCatch({
      sub_df$tx_status = factor(sub_df$tx_status, levels = c("Treated", "Baseline"))
      # Perform Wilcoxon test to get p-value
      wilcox_res <- wilcox.test(freq_main ~ tx_status, data = sub_df, exact = FALSE, conf.int = TRUE)
      # Calculate rank-biserial correlation with confidence intervals
      rbc <- rcompanion::wilcoxonRG(x = sub_df$freq_main, 
                                    g = sub_df$tx_status, 
                                    ci = TRUE,  
                                    conf = 0.95,
                                    R = 1000)   
      wilcox_per_cohort_results_list[[length(wilcox_per_cohort_results_list) + 1]] <- tibble(
        cohort = cohort_name,
        celltype_main = celltype,
        p_value = wilcox_res$p.value,
        effect_size = rbc$r,
        effect_se = (rbc$upper - rbc$lower) / (2 * 1.96),
        n_baseline = sum(sub_df$tx_status == 'Baseline'),
        n_traeted = sum(sub_df$tx_status == 'Treated')
      )
    }, error = function(e) {
      message(paste("Error processing", celltype, "in cohort", cohort_name, ":", e$message))
    })
  }
}
final_cohort_results <- bind_rows(wilcox_per_cohort_results_list)
write.csv(final_cohort_results |> data.frame(), 'tables/wicox_main_per_cohort_tx_status.csv', row.names = F)
final_cohort_results <- read.csv('tables/wicox_main_per_cohort_tx_status.csv')

# Meta-analysis
set.seed(1234)
wilcox_meta_results_list <- list()
for (celltype in unique(final_cohort_results$celltype_main)) {
  print(celltype)
  # Filter data for this cell type across all cohorts
  cell_data <- final_cohort_results |> filter(celltype_main == celltype)
  # Only proceed if we have results from at least 2 cohorts
  if (nrow(cell_data) >= 2) {
    # Perform proper random-effects meta-analysis
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
        celltype_main = celltype,
        pooled_effect = meta_res$b[1],
        pooled_effect_se = meta_res$se,
        ci_lower = meta_res$ci.lb,
        ci_upper = meta_res$ci.ub,
        p_value = meta_res$pval,
        i2 = meta_res$I2,
        tau2 = meta_res$tau2,
        n_cohorts = nrow(cell_data)
      )
    }
  }
}
final_results <- bind_rows(wilcox_meta_results_list) |> arrange(desc(pooled_effect))
roworder <- final_results |> arrange(desc(pooled_effect)) |> filter(n_cohorts > 10) |> pull(celltype_main)
write.csv(final_results |> arrange(desc(pooled_effect)) |> data.frame(), 'tables/wicox_main_per_cohort_tx_status_meta.csv', row.names = F)
final_results <- read.csv('tables/wicox_main_per_cohort_tx_status_meta.csv')
library(grid)
pdf('Dynamics/fp_main_tx_status.pdf', height = 5, width = 4.5)
xlim_min <- -0.3
xlim_max <- 0.3
plot_data <- final_results |>
  arrange(desc(pooled_effect)) |>
  filter(n_cohorts > 10) |>
  mutate(
    # Create new columns for the clipped confidence intervals
    ci_lower_clipped = pmax(ci_lower, xlim_min),
    ci_upper_clipped = pmin(ci_upper, xlim_max),
    
    # Create logical flags to determine if arrows are needed
    needs_left_arrow = ci_lower < xlim_min,
    needs_right_arrow = ci_upper > xlim_max
  )
p <- ggplot(plot_data, aes(x = pooled_effect, y = reorder(celltype_main, pooled_effect))) +
  
  # Draw the clipped error bars using the new columns
  geom_errorbarh(aes(xmin = ci_lower_clipped, xmax = ci_upper_clipped), height = 0.2, color = 'black') +
  
  # Add the left-pointing arrows for CIs that exceed the lower limit
  geom_segment(
    data = subset(plot_data, needs_left_arrow),
    aes(x = xlim_min, xend = xlim_min - 0.01, y = reorder(celltype_main, pooled_effect), yend = reorder(celltype_main, pooled_effect)),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "black"
  ) +
  
  # Add the right-pointing arrows for CIs that exceed the upper limit
  geom_segment(
    data = subset(plot_data, needs_right_arrow),
    aes(x = xlim_max, xend = xlim_max + 0.01, y = reorder(celltype_main, pooled_effect), yend = reorder(celltype_main, pooled_effect)),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"), color = "black"
  ) +
  
  # Plot the pooled effect points
  geom_point(aes(color = ifelse(pooled_effect > 0, "Treated", "Baseline")), size = 3) +
  
  # Add the vertical line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  
  # Set the colors and axis limits
  scale_color_manual(values = c("Baseline" = "#04a3bd", "Treated" = "#f0be3d")) +
  
  # Customize the plot labels and title
  labs(
    x = "Effect Size \n(pooled)",
    y = "Cell Type",
    color = "Effect Direction"
  ) +
  
  # Use a clean theme and adjust text size for better readability
  theme_minimal() +
  theme(plot.margin = unit(c(0.5, 2.5, 1.5, 0.5), "lines"),
        axis.text.x = element_text(size = 10, color = 'black'),
        axis.title.x = element_text(size = 14, color = 'black', vjust = -5),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title.y = element_blank(),
        legend.position = "none"
  ) +
  # Set the x-axis limits and ensure clipping is off for arrows
  coord_cartesian(xlim = c(xlim_min, xlim_max), clip = "off")

grid.draw(p)
y_text = 0.06
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
grid.text("Baseline", x = unit(x1, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 10))
grid.text("Treated", x = unit(x2, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 10))
dev.off()
