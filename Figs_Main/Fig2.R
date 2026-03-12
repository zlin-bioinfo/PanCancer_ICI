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
# Per-cohort PAIRED analysis
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggnewscale)

# 1. Helper Functions
# Compute common row order across subgroups
# Merges cell types from both groups and sorts by pooled effect (1-cycle preferred)
compute_common_roworder <- function(celltypes_subset,
                                     results_1cycle, results_2plus,
                                     min_cohorts_1cycle, min_cohorts_2plus) {
  ct_1cycle <- results_1cycle |>
    filter(celltype_r2 %in% celltypes_subset, n_cohorts >= min_cohorts_1cycle) |>
    pull(celltype_r2)
  ct_2plus <- results_2plus |>
    filter(celltype_r2 %in% celltypes_subset, n_cohorts >= min_cohorts_2plus) |>
    pull(celltype_r2)
  common_ct <- union(ct_1cycle, ct_2plus)
  if (length(common_ct) == 0) return(character(0))

  eff_1 <- results_1cycle |>
    filter(celltype_r2 %in% common_ct) |>
    select(celltype_r2, pooled_effect)
  eff_2 <- results_2plus |>
    filter(celltype_r2 %in% common_ct) |>
    select(celltype_r2, pooled_effect)
  combined <- eff_1 |>
    full_join(eff_2, by = "celltype_r2", suffix = c("_1cycle", "_2plus")) |>
    mutate(sort_effect = ifelse(!is.na(pooled_effect_1cycle), pooled_effect_1cycle, pooled_effect_2plus)) |>
    arrange(sort_effect)

  return(unique(combined$celltype_r2))
}

# Calculate plot dimensions from cell type and cohort counts
calc_dimensions <- function(n_celltypes, n_cohorts) {
  height <- 1.9 + n_celltypes * 0.25
  width <- 3 + n_cohorts * 0.33
  return(list(width = width, height = height))
}

# Compute shared xlim across two subgroups for a set of cell types
compute_shared_xlim <- function(results_1, results_2, common_roworder,
                                 min_cohorts_1, min_cohorts_2) {
  ci_1 <- results_1 |> filter(celltype_r2 %in% common_roworder, n_cohorts >= min_cohorts_1)
  ci_2 <- results_2 |> filter(celltype_r2 %in% common_roworder, n_cohorts >= min_cohorts_2)
  ci_min <- min(c(ci_1$ci_lower, ci_2$ci_lower), na.rm = TRUE)
  ci_max <- max(c(ci_1$ci_upper, ci_2$ci_upper), na.rm = TRUE)
  ci_range <- ci_max - ci_min
  c(min(ci_min - 0.05 * ci_range, -0.05),
    max(ci_max + 0.05 * ci_range, 0.05))
}

# Compute xlim from a data frame with ci_lower and ci_upper columns
compute_xlim <- function(ci_lower, ci_upper, pad = 0.05, min_abs = 0.1) {
  ci_min <- min(ci_lower, na.rm = TRUE)
  ci_max <- max(ci_upper, na.rm = TRUE)
  ci_rng <- ci_max - ci_min
  c(min(ci_min - pad * ci_rng, -min_abs),
    max(ci_max + pad * ci_rng,  min_abs))
}

# Generate striped background data (ggforestplot-style alternating rows)
make_stripe_df <- function(n_rows) {
  data.frame(
    ymin = seq(0.5, n_rows + 0.5, by = 2),
    ymax = pmin(seq(1.5, n_rows + 1.5, by = 2), n_rows + 0.5)
  )
}

# 2. Overlay forest plot (Pooled + 1c + 2+c on same axes)
# Build overlay forest plot with striped background
build_overlay_forest <- function(merged_data,
                                  roworder,
                                  xlim,
                                  color_bl = "#04a3bd",
                                  color_tx = "#f0be3d",
                                  shape_vals = c("1 cycle" = 16, "2+ cycles" = 17),
                                  point_size = 2.5,
                                  dodge_width = 0.55,
                                  stripe_fill = "#ededed",
                                  base_size = 12) {

  n_ct <- length(roworder)
  stripe_df <- make_stripe_df(n_ct)

  # Clip CIs to xlim
  plot_data <- merged_data |>
    mutate(
      ci_lower_clip = pmax(ci_lower, xlim[1]),
      ci_upper_clip = pmin(ci_upper, xlim[2]),
      direction = ifelse(pooled_effect > 0, "Treated > BL", "BL > Treated"),
      pt_alpha  = ifelse(p_value < 0.05, 1, 0.35)
    )

  x_breaks <- pretty(xlim, n = 5)

  ggplot(plot_data, aes(x = pooled_effect, y = celltype_r2, group = cycle)) +
    geom_rect(data = stripe_df,
              aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
              fill = stripe_fill, color = NA, inherit.aes = FALSE) +
    geom_hline(yintercept = c(0.5, n_ct + 0.5), color = "black", linewidth = 0.3) +
    geom_vline(xintercept = x_breaks, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_errorbarh(
      aes(xmin = ci_lower_clip, xmax = ci_upper_clip),
      height = 0.08, linewidth = 0.35, color = "black",
      position = position_dodge(width = dodge_width, reverse = TRUE)
    ) +
    geom_point(aes(shape = cycle, color = direction, alpha = pt_alpha),
               size = point_size,
               position = position_dodge(width = dodge_width, reverse = TRUE)) +
    scale_color_manual(
      values = c("BL > Treated" = color_bl, "Treated > BL" = color_tx),
      guide = "none"
    ) +
    scale_shape_manual(
      values = shape_vals, name = "Treatment cycle",
      guide = guide_legend(override.aes = list(color = "black", alpha = 1, size = point_size))
    ) +
    scale_alpha_identity() +
    scale_x_continuous(
      breaks = x_breaks,
      labels = function(x) sprintf("%.1f", x)
    ) +
    scale_y_discrete(drop = FALSE) +
    coord_cartesian(xlim = xlim, clip = "on") +
    labs(x = "Pooled effect size\n(rank-biserial r)", y = NULL) +
    theme_bw(base_size = base_size) +
    theme(
      panel.background   = element_rect(fill = "white", color = NA),
      panel.border       = element_blank(),
      axis.text.y        = element_text(size = 11, color = "black"),
      axis.text.x        = element_text(size = 10, color = "black"),
      axis.title.x       = element_text(size = 11, color = "black", margin = margin(t = 8)),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin        = margin(t = 5, r = 5, b = 25, l = 5)
    )
}

# 3. Combined dotplot (all cohorts, color = effect, size = confidence)
# Build combined dotplot for per-cohort results
build_overlay_dotplot <- function(dot_data,
                                   roworder,
                                   cohort_order,
                                   color_bl = "#04a3bd",
                                   color_tx = "#f0be3d",
                                   effect_clip = 0.7,
                                   base_size = 11) {

  dd <- dot_data |>
    mutate(
      celltype_r2    = factor(celltype_r2, levels = roworder),
      cohort         = factor(cohort, levels = cohort_order),
      effect_clipped = pmax(-effect_clip, pmin(effect_clip, effect_size)),
      weight         = pmax(1, pmin(6, 1 / ifelse(is.na(effect_se), 1e10, effect_se))),
      significant    = p_value < 0.05
    )

  ggplot(dd, aes(x = cohort, y = celltype_r2)) +
    geom_point(aes(color = effect_clipped, size = weight), na.rm = TRUE) +
    geom_point(data = subset(dd, significant),
               aes(size = weight), shape = 21, color = "black", fill = NA, stroke = 0.5) +
    scale_color_gradient2(
      low = color_bl, mid = "white", high = color_tx,
      midpoint = 0, limits = c(-effect_clip, effect_clip),
      na.value = "lightgray", name = "Effect size"
    ) +
    scale_size_continuous(range = c(0, 5), name = "Confidence",
                          limits = c(0, 6), breaks = c(2, 4, 6)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9, color = "black"),
      axis.text.y    = element_blank(),
      axis.ticks.y   = element_line(color = "black"),
      panel.grid.minor = element_blank(),
      plot.margin    = margin(t = 5, r = 0, b = 25, l = 0)
    )
}

# 4. Combined annotation (cancer type + treatment + cycle bars above dotplot)
build_overlay_annotation <- function(ann_data,
                                      cancertype_colors,
                                      treatment_colors,
                                      cycle_colors = c("1 cycle" = "#a6cee3", "2+ cycles" = "#fb9a99"),
                                      all_subtypes = NULL,
                                      all_treatments = NULL,
                                      label_margin_left = 55,
                                      base_size = 11) {

  if (is.null(all_subtypes))  all_subtypes  <- sort(unique(ann_data$subtype))
  if (is.null(all_treatments)) all_treatments <- sort(unique(ann_data$treatment))

  ggplot(ann_data, aes(x = cohort)) +
    geom_tile(aes(y = 1, fill = subtype), color = "white", linewidth = 0.5) +
    scale_fill_manual(values = cancertype_colors, name = "Cancer type",
                      limits = all_subtypes) +
    new_scale_fill() +
    geom_tile(aes(y = 2, fill = treatment), color = "white", linewidth = 0.5) +
    scale_fill_manual(values = treatment_colors, name = "Treatment",
                      limits = all_treatments) +
    new_scale_fill() +
    geom_tile(aes(y = 3, fill = cycle_group), color = "white", linewidth = 0.5) +
    scale_fill_manual(values = cycle_colors, name = "Treatment cycle") +
    annotate("text", x = 0.35, y = 1, label = "Cancer",    hjust = 1, size = 5) +
    annotate("text", x = 0.35, y = 2, label = "Treatment", hjust = 1, size = 5) +
    annotate("text", x = 0.35, y = 3, label = "Treatment cycle",     hjust = 1, size = 5) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = base_size) +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(), panel.border = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = label_margin_left)
    )
}

# 5. Patchwork composition (forest + dotplot + annotation)
# Compose overlay summary: annotation on top-right, forest left, dotplot right
compose_overlay <- function(p_forest, p_dotplot, p_annot,
                             n_celltypes, n_cohorts,
                             w_forest = 7,
                             h_annot = 0.6,
                             w_per_cohort = 0.33,
                             h_per_ct = 0.25,
                             h_base = 1.1,
                             w_legend = 3) {

  h_content <- n_celltypes * h_per_ct + h_base
  w_dot     <- n_cohorts * w_per_cohort

  design <- "
#A
BC
"

  p_combined <- p_annot +
    p_forest + p_dotplot +
    plot_layout(
      design = design,
      heights = c(h_annot, h_content),
      widths  = c(w_forest, w_dot),
      guides  = "collect"
    ) &
    theme(legend.position = "right")

  total_w <- w_forest + w_dot + w_legend
  total_h <- h_annot + h_content

  return(list(plot = p_combined, width = total_w, height = total_h))
}


# 6. Single-subgroup forest+dotplot+annotation (for side-by-side layout)
# Build forest+dotplot+annotation for one subgroup (1-cycle or 2+-cycle)
build_summary_plot <- function(final_results_sub, final_cohort_results_sub, df_sub,
                                common_roworder, title_label,
                                min_cohorts_sub, shared_xlim = NULL,
                                effect_clip = 0.7,
                                color_baseline = "#04a3bd",
                                color_treated = "#f0be3d",
                                master_cohort_order = NULL,
                                cancertype_colors = NULL,
                                treatment_colors = NULL) {

  plot_data_raw <- final_results_sub |>
    filter(celltype_r2 %in% common_roworder, n_cohorts >= min_cohorts_sub) |>
    arrange(pooled_effect)

  if (nrow(plot_data_raw) == 0) return(NULL)

  # xlim
  if (!is.null(shared_xlim)) {
    xlim_min <- shared_xlim[1]
    xlim_max <- shared_xlim[2]
  } else {
    ci_min <- min(plot_data_raw$ci_lower, na.rm = TRUE)
    ci_max <- max(plot_data_raw$ci_upper, na.rm = TRUE)
    ci_range <- ci_max - ci_min
    xlim_min <- min(ci_min - 0.05 * ci_range, -0.05)
    xlim_max <- max(ci_max + 0.05 * ci_range,  0.05)
  }

  x_breaks <- pretty(c(xlim_min, xlim_max), n = 3)
  if (!0 %in% x_breaks) x_breaks <- sort(unique(c(0, x_breaks)))
  x_breaks <- x_breaks[x_breaks >= xlim_min & x_breaks <= xlim_max]

  plot_data <- plot_data_raw |>
    mutate(
      sig_label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      ci_lower_clipped  = pmax(ci_lower, xlim_min),
      ci_upper_clipped  = pmin(ci_upper, xlim_max),
      needs_left_arrow  = ci_lower < xlim_min,
      needs_right_arrow = ci_upper > xlim_max,
      celltype_r2 = factor(celltype_r2, levels = common_roworder),
      significant = p_value < 0.05
    )

  forest_plot <- ggplot(plot_data, aes(x = pooled_effect, y = celltype_r2)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = ci_lower_clipped, xmax = ci_upper_clipped),
                   height = 0.2, color = "black", linewidth = 0.4) +
    geom_segment(data = subset(plot_data, needs_left_arrow),
                 aes(x = xlim_min, xend = xlim_min - 0.02, y = celltype_r2, yend = celltype_r2),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black") +
    geom_segment(data = subset(plot_data, needs_right_arrow),
                 aes(x = xlim_max, xend = xlim_max + 0.02, y = celltype_r2, yend = celltype_r2),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black") +
    geom_point(data = subset(plot_data, significant),
               aes(color = ifelse(pooled_effect > 0, "Treated", "Baseline")),
               size = 3, shape = 16) +
    geom_point(data = subset(plot_data, !significant),
               aes(color = ifelse(pooled_effect > 0, "Treated", "Baseline")),
               size = 3, shape = 1, stroke = 1) +
    geom_text(aes(x = xlim_max + 0.06, label = sig_label),
              hjust = 0, size = 6, color = "black") +
    scale_color_manual(values = c("Baseline" = color_baseline, "Treated" = color_treated),
                       guide = "none") +
    scale_x_continuous(breaks = x_breaks, labels = function(x) sprintf("%.1f", x)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Pooled Effect Size \n(rank-biserial r)", y = NULL, title = "") +
    coord_cartesian(xlim = c(xlim_min, xlim_max), clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      plot.margin      = margin(t = 5, r = 30, b = 25, l = 5),
      axis.text.x      = element_text(size = 10, color = "black"),
      axis.title.x     = element_text(size = 11, color = "black", margin = margin(t = 8)),
      axis.text.y      = element_text(size = 11, color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )

  # Dotplot
  cohorts_in_data <- final_cohort_results_sub |>
    filter(celltype_r2 %in% common_roworder) |> pull(cohort) |> unique()
  if (!is.null(master_cohort_order)) {
    cohort_order_filtered <- master_cohort_order[master_cohort_order %in% cohorts_in_data]
  } else {
    cohort_order_filtered <- cohorts_in_data
  }

  cohort_treatment <- data.frame(
    cohort    = cohort_order_filtered,
    treatment = df_sub$treatment[match(cohort_order_filtered, df_sub$cohort)]
  )

  dotplot_data <- final_cohort_results_sub |>
    filter(celltype_r2 %in% common_roworder) |>
    left_join(cohort_treatment, by = "cohort") |>
    mutate(
      celltype_r2    = factor(celltype_r2, levels = common_roworder),
      cohort         = factor(cohort, levels = cohort_order_filtered),
      effect_clipped = pmax(-effect_clip, pmin(effect_clip, effect_size)),
      weight         = pmax(1, pmin(6, 1 / ifelse(is.na(effect_se), 1e+10, effect_se))),
      significant    = p_value < 0.05
    )

  dotplot <- ggplot(dotplot_data, aes(x = cohort, y = celltype_r2)) +
    geom_point(aes(color = effect_clipped, size = weight), na.rm = TRUE) +
    geom_point(data = subset(dotplot_data, significant),
               aes(size = weight), shape = 21, color = "black", fill = NA, stroke = 0.5) +
    scale_color_gradient2(low = color_baseline, mid = "white", high = color_treated,
                          midpoint = 0, limits = c(-effect_clip, effect_clip),
                          na.value = "lightgray", name = "Effect size") +
    scale_size_continuous(range = c(0, 5), name = "Confidence",
                          limits = c(0, 6), breaks = c(2, 4, 6)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9, color = "black"),
      axis.text.y    = element_blank(),
      axis.ticks.y   = element_line(color = "black"),
      panel.grid.minor = element_blank(),
      plot.margin    = margin(t = 5, r = 0, b = 25, l = 0)
    )

  # Annotation bars
  annotation_data <- data.frame(
    cohort    = factor(cohort_order_filtered, levels = cohort_order_filtered),
    treatment = df_sub$treatment[match(cohort_order_filtered, df_sub$cohort)],
    subtype   = df_sub$subtype[match(cohort_order_filtered, df_sub$cohort)]
  )

  combined_annotation <- ggplot(annotation_data, aes(x = cohort)) +
    geom_tile(aes(y = 1, fill = subtype), color = "white", linewidth = 0.5) +
    scale_fill_manual(values = cancertype_colors, name = "Cancer type") +
    new_scale_fill() +
    geom_tile(aes(y = 2, fill = treatment), color = "white", linewidth = 0.5) +
    scale_fill_manual(values = treatment_colors, name = "Treatment") +
    annotate("text", x = 0.35, y = 1, label = "Cancer type", hjust = 1, size = 2.8) +
    annotate("text", x = 0.35, y = 2, label = "Treatment",   hjust = 1, size = 2.8) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_line(color = "black"),
      panel.grid = element_blank(), panel.border = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 55)
    )

  return(list(forest = forest_plot, dotplot = dotplot, annotation = combined_annotation,
              title_label = title_label,
              n_celltypes = length(common_roworder), n_cohorts = length(cohort_order_filtered)))
}

# Compose side-by-side layout from two subgroup component lists
compose_combined <- function(result_1, result_2, title, w_forest, h_annot = 0.5) {
  if (is.null(result_1$forest) || is.null(result_2$forest)) return(NULL)

  n_cohorts_1 <- result_1$n_cohorts
  n_cohorts_2 <- result_2$n_cohorts

  title_1 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = result_1$title_label,
                                  size = 4.5, fontface = "bold") +
    theme_void() + coord_cartesian(clip = "off")
  title_2 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = result_2$title_label,
                                  size = 4.5, fontface = "bold") +
    theme_void() + coord_cartesian(clip = "off")

  design <- "AABB\nCDEF\nGHIJ"

  h_title_in   <- 0.3
  h_annot_in   <- h_annot
  n_celltypes  <- result_1$n_celltypes
  h_content_in <- n_celltypes * 0.25 + 1.1
  height <- h_title_in + h_annot_in + h_content_in
  width  <- 3 + (n_cohorts_1 + n_cohorts_2) * 0.33 + 4

  p <- title_1 + title_2 +
    plot_spacer() + result_1$annotation +
    plot_spacer() + result_2$annotation +
    result_1$forest + result_1$dotplot +
    result_2$forest + result_2$dotplot +
    plot_layout(design = design,
                heights = c(h_title_in, h_annot_in, h_content_in),
                widths = c(w_forest, n_cohorts_1, w_forest, n_cohorts_2),
                guides = "collect") +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))) &
    theme(legend.position = "right")

  return(list(plot = p, width = width, height = height))
}

# plotting
source('scripts/Celltype_classification.R')
set.seed(1234)
metadata <- read.csv('tables/meta_all.csv')

unwanted_celltypes <- c('Melanocytes(CNA-)', 'Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Epithelial(CNA-)', 'Malignant(CNA+)')
# patients in Becker without matched tumors
unmatched_pt <- c("SKCM_this study_Patient2", "SKCM_this study_Patient3")
df <- metadata |>
  filter(count_immune >= 50, non_malignant_count >= 100,
         !patient %in% unmatched_pt,
         tx_status %in% c("Baseline", "Treated"),
         subset %in% c('TME', 'CD45+sorted')) |>
  filter(!celltype_r2 %in% unwanted_celltypes) |>
  distinct(sample, celltype_r2, freq_r2_comp, .keep_all = TRUE)

# Cohort adjustments for subtype-specific analysis
df$cohort[df$cohort == 'BRCA_Bassez1' & df$subtype == 'TNBC'] <- 'TNBC_Bassez1'
df$cohort[df$cohort == 'BRCA_Bassez2' & df$subtype == 'TNBC'] <- 'TNBC_Bassez2'
df$cohort[df$cohort == 'CRC_Chen' & df$tx_cycle == "2 cycles"] <- 'CRC_Chen_2cycles'
df$cohort[df$cohort == 'CRC_Chen' & df$tx_cycle == "3+ cycles"] <- 'CRC_Chen_3cycles+'

# Pool SKCM cohorts and split by mono vs combo
df$cohort[df$cohort %in% c('SKCM_this study', 'SKCM_Plozniak') & df$treatment == 'aPD1'] <- 'SKCM_Mono'
df$cohort[df$cohort %in% c('SKCM_this study', 'SKCM_Plozniak') & df$treatment == 'aPD1+CTLA4'] <- 'SKCM_Combo'

# Split HNSCC_Franken by mono vs combo
df$cohort[df$cohort == 'HNSCC_Franken' & df$treatment == 'aPDL1'] <- 'HNSCC_Franken_Mono'
df$cohort[df$cohort == 'HNSCC_Franken' & df$treatment == 'aPDL1+CTLA4'] <- 'HNSCC_Franken_Combo'

# HNSCC_Luoma too small to split (3+3), label as aPD1±CTLA4
df$treatment[df$cohort == 'HNSCC_Luoma'] <- 'aPD1\u00b1CTLA4'

# Colors (consistent with r2_overview.R)
metadata_for_colors <- metadata |>
  mutate(treatment = case_when(
    cohort %in% c('HNSCC_Luoma', 'SKCM_Plozniak', 'SKCM_this study') ~ 'aPD1\u00b1CTLA4',
    TRUE ~ treatment
  ))
metadata_for_colors$treatment[metadata_for_colors$cohort == 'HNSCC_Franken'] <- 'aPDL1\u00b1CTLA4'
treatment_colors <- brewer.pal(n = length(unique(metadata_for_colors$treatment)), name = "Dark2")
names(treatment_colors) <- unique(metadata_for_colors$treatment)
cancerytpe_colors <- as.character(met.brewer("Austria", n = length(unique(metadata$subtype))))
names(cancerytpe_colors) <- unique(metadata$subtype)

df$proportion <- df$freq_r2_comp
df <- df |> mutate(
  tx_status = factor(tx_status, levels = c("Baseline", "Treated"))
)
celltypes <- unique(df$celltype_r2)

# Subgroup analysis by tx_cycle: 1 cycle vs 2+ cycles
# Exclude cohorts with "Not specified" tx_cycle
df <- df |> filter(tx_cycle != "Not specified")

# Create tx_cycle_group variable (tx_cycle is character)
df <- df |>
  mutate(tx_cycle_group = case_when(
    tx_cycle == "1 cycle" ~ "1 cycle",
    tx_cycle %in% c("2 cycles", "2-4 cycles", "3 cycles", "3+ cycles", "4 cycles") ~ "2+ cycles",
    TRUE ~ NA_character_
  ))

df_1cycle <- df |> filter(tx_cycle_group == "1 cycle" | tx_status == "Baseline") |>
  filter(cohort %in% unique(final_cohort_results_1cycle$cohort))
df_2plus <- df |> filter(tx_cycle_group == "2+ cycles" | tx_status == "Baseline") |>
  filter(cohort %in% unique(final_cohort_results_2plus$cohort))

# Analysis for 1 cycle subgroup
# Identify cohorts that have treated samples with 1 cycle
cohorts_with_1cycle <- df |>
  filter(tx_status == "Treated", tx_cycle_group == "1 cycle") |>
  pull(cohort) |>
  unique()

# Filter to only cohorts with 1 cycle data, then select Baseline + Treated (1 cycle)
df_1cycle <- df |>
  filter(cohort %in% cohorts_with_1cycle) |>
  filter(tx_status == "Baseline" | (tx_status == "Treated" & tx_cycle_group == "1 cycle"))

# Identify patients with matched Baseline and Treated (1 cycle) samples
matched_patients_1cycle <- df_1cycle |>
  dplyr::group_by(cohort, patient) |>
  dplyr::summarise(has_both = dplyr::n_distinct(tx_status) == 2, .groups = "drop") |>
  dplyr::filter(has_both) |>
  dplyr::pull(patient)

# Per-cohort PAIRED Wilcoxon tests - 1 cycle
wilcox_per_cohort_results_list_1cycle <- list()
for (cohort_name in unique(df_1cycle$cohort)) {
  print(cohort_name)
  for (celltype in setdiff(unique(df_1cycle$celltype_r2), c('pDC', 'Mast'))) {
    sub_df <- df_1cycle |>
      filter(cohort == cohort_name, celltype_r2 == celltype, patient %in% matched_patients_1cycle)

    # Check for paired data structure
    if(length(unique(sub_df$tx_status)) < 2 || nrow(sub_df) < 2) next

    tryCatch({
      # Pivot to wide format for paired analysis
      wide_df <- sub_df |>
        select(patient, tx_status, freq_r2_comp) |>
        pivot_wider(names_from = tx_status, values_from = freq_r2_comp, values_fill = 0)

      # Ensure both time points are present
      if(!"Baseline" %in% colnames(wide_df)) wide_df$Baseline <- 0
      if(!"Treated" %in% colnames(wide_df)) wide_df$Treated <- 0

      # Check if we have enough paired subjects
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

      # SE calculation using correlation approximation
      effect_se <- sqrt((1 - effect_size^2) / (n_pairs - 2))

      wilcox_per_cohort_results_list_1cycle[[length(wilcox_per_cohort_results_list_1cycle) + 1]] <- tibble(
        cohort = cohort_name,
        celltype_r2 = celltype,
        p_value = wilcox_res$p.value,
        effect_size = effect_size,
        effect_se = effect_se,
        n_pairs = n_pairs,
        median_baseline = median(wide_df$Baseline),
        median_treated = median(wide_df$Treated),
        mean_diff = mean(differences)
      )
    }, error = function(e) {
      message(paste("Error processing", celltype, "in cohort", cohort_name, ":", e$message))
    })
  }
}
final_cohort_results_1cycle <- bind_rows(wilcox_per_cohort_results_list_1cycle)

# Meta-analysis - 1 cycle
# Filter out problematic rows: n_pairs < 3 (SE undefined), SE = 0/Inf/NA
final_cohort_results_1cycle_filtered <- final_cohort_results_1cycle |>
  filter(n_pairs >= 3,
         is.finite(effect_se),
         effect_se > 0)

wilcox_meta_results_list_1cycle <- list()
for (celltype in unique(final_cohort_results_1cycle_filtered$celltype_r2)) {
  print(celltype)
  cell_data <- final_cohort_results_1cycle_filtered |> filter(celltype_r2 == celltype)
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
      wilcox_meta_results_list_1cycle[[celltype]] <- tibble(
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
final_results_1cycle <- bind_rows(wilcox_meta_results_list_1cycle) |> arrange(desc(pooled_effect))

# Analysis for 2+ cycles subgroup
# Identify cohorts that have treated samples with 2+ cycles
cohorts_with_2plus <- df |>
  dplyr::filter(tx_status == "Treated", tx_cycle_group == "2+ cycles") |>
  dplyr::pull(cohort) |>
  unique()

# Filter to only cohorts with 2+ cycles data, then select Baseline + Treated (2+ cycles)
df_2plus <- df |>
  dplyr::filter(cohort %in% cohorts_with_2plus) |>
  dplyr::filter(tx_status == "Baseline" | (tx_status == "Treated" & tx_cycle_group == "2+ cycles"))

# Identify patients with matched Baseline and Treated (2+ cycles) samples
matched_patients_2plus <- df_2plus |>
  dplyr::group_by(cohort, patient) |>
  dplyr::summarise(has_both = dplyr::n_distinct(tx_status) == 2, .groups = "drop") |>
  dplyr::filter(has_both) |>
  dplyr::pull(patient)

# Per-cohort PAIRED Wilcoxon tests - 2+ cycles
wilcox_per_cohort_results_list_2plus <- list()
for (cohort_name in unique(df_2plus$cohort)) {
  print(cohort_name)
  for (celltype in setdiff(unique(df_2plus$celltype_r2), c('pDC', 'Mast'))) {
    sub_df <- df_2plus |>
      dplyr::filter(cohort == cohort_name, celltype_r2 == celltype, patient %in% matched_patients_2plus)

    if(length(unique(sub_df$tx_status)) < 2 || nrow(sub_df) < 2) next

    tryCatch({
      wide_df <- sub_df |>
        dplyr::select(patient, tx_status, freq_r2_comp) |>
        pivot_wider(names_from = tx_status, values_from = freq_r2_comp, values_fill = 0)

      if(!"Baseline" %in% colnames(wide_df)) wide_df$Baseline <- 0
      if(!"Treated" %in% colnames(wide_df)) wide_df$Treated <- 0
      if(nrow(wide_df) < 2) next

      wilcox_res <- wilcox.test(wide_df$Treated, wide_df$Baseline,
                                paired = TRUE, exact = FALSE, conf.int = TRUE)

      n_pairs <- nrow(wide_df)
      differences <- wide_df$Treated - wide_df$Baseline

      # Matched pairs rank biserial correlation
      effect_size <- wilcox_res$statistic / (n_pairs * (n_pairs + 1) / 2) - 0.5
      effect_size <- effect_size * 2  # Scale to range [-1, 1]

      # SE calculation using correlation approximation
      effect_se <- sqrt((1 - effect_size^2) / (n_pairs - 2))

      wilcox_per_cohort_results_list_2plus[[length(wilcox_per_cohort_results_list_2plus) + 1]] <- tibble(
        cohort = cohort_name,
        celltype_r2 = celltype,
        p_value = wilcox_res$p.value,
        effect_size = effect_size,
        effect_se = effect_se,
        n_pairs = n_pairs,
        median_baseline = median(wide_df$Baseline),
        median_treated = median(wide_df$Treated),
        mean_diff = mean(differences)
      )
    }, error = function(e) {
      message(paste("Error processing", celltype, "in cohort", cohort_name, ":", e$message))
    })
  }
}
final_cohort_results_2plus <- bind_rows(wilcox_per_cohort_results_list_2plus)

# Meta-analysis - 2+ cycles
# Filter out problematic rows: n_pairs < 3 (SE undefined), SE = 0/Inf/NA
final_cohort_results_2plus_filtered <- final_cohort_results_2plus |>
  filter(n_pairs >= 3,
         is.finite(effect_se),
         effect_se > 0)

print(paste("Cohorts with 2+ cycle data (after filtering):", paste(unique(final_cohort_results_2plus_filtered$cohort), collapse = ", ")))
print(paste("Number of cohorts:", length(unique(final_cohort_results_2plus_filtered$cohort))))

wilcox_meta_results_list_2plus <- list()
for (celltype in unique(final_cohort_results_2plus_filtered$celltype_r2)) {
  print(celltype)
  cell_data <- final_cohort_results_2plus_filtered |> dplyr::filter(celltype_r2 == celltype)
  print(paste("  -> n_cohorts:", nrow(cell_data)))
  if (nrow(cell_data) >= 2) {
    meta_res <- tryCatch({
      rma(yi = effect_size, sei = effect_se, data = cell_data, method = "REML")
    }, error = function(e) {
      message(paste("Meta-analysis failed for", celltype, ":", e$message))
      NULL
    })
    if (!is.null(meta_res)) {
      wilcox_meta_results_list_2plus[[celltype]] <- tibble(
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
final_results_2plus <- bind_rows(wilcox_meta_results_list_2plus) |> arrange(desc(pooled_effect))
print(paste("Cell types with meta-analysis results:", nrow(final_results_2plus)))


# Keep only cell types that converged in BOTH subgroups (for comparability)
converged_both <- intersect(final_results_1cycle$celltype_r2, final_results_2plus$celltype_r2)
removed_1cycle <- setdiff(final_results_1cycle$celltype_r2, converged_both)
removed_2plus <- setdiff(final_results_2plus$celltype_r2, converged_both)
if (length(removed_1cycle) > 0) print(paste("Removed from 1-cycle (not converged in 2+):", paste(removed_1cycle, collapse = ", ")))
if (length(removed_2plus) > 0) print(paste("Removed from 2+ (not converged in 1-cycle):", paste(removed_2plus, collapse = ", ")))
final_results_1cycle <- final_results_1cycle |> filter(celltype_r2 %in% converged_both)
final_results_2plus <- final_results_2plus |> filter(celltype_r2 %in% converged_both)
final_cohort_results_1cycle <- final_cohort_results_1cycle |> filter(celltype_r2 %in% converged_both)
final_cohort_results_2plus <- final_cohort_results_2plus |> filter(celltype_r2 %in% converged_both)
print(paste("Cell types retained in both:", length(converged_both)))


# Plotting 
min_cohorts_1cycle <- 3
min_cohorts_2plus <- 2
effect_clip <- 0.7
color_baseline <- "#04a3bd"
color_treated <- "#f0be3d"
plot_widths <- c(0.5, 0.5)
w_forest <- 6  # width of forest panels in combined summary plot

# xlim will be calculated dynamically based on CIs in each plot

# Define master cohort order
master_cohort_order <- c("SKCM_Mono", "SKCM_Combo",
                         "BRCA_Bassez1", "BRCA_Bassez2",
                         "TNBC_Bassez1", "TNBC_Bassez2", "TNBC_Shiao",
                         "HNSCC_Franken_Mono", "HNSCC_Franken_Combo", "CRC_Chen",
                         "TNBC_Zhang", "HNSCC_Luoma", "HNSCC_vanderLeun",
                         "CRC_Chen_2cycles", "CRC_Chen_3cycles+", "CRC_Li",
                         "NSCLC_Yan", "HCC_Guo", "PCa_Hawley")

# Helper function to compute common row order for a lineage across both subgroups
# Merges cell types from both groups and sorts by pooled effect (1-cycle preferred, fallback to 2+)
compute_common_roworder <- function(celltypes_subset,
                                     results_1cycle, results_2plus,
                                     min_cohorts_1cycle, min_cohorts_2plus) {
  ct_1cycle <- results_1cycle |>
    filter(celltype_r2 %in% celltypes_subset, n_cohorts >= min_cohorts_1cycle) |>
    pull(celltype_r2)
  ct_2plus <- results_2plus |>
    filter(celltype_r2 %in% celltypes_subset, n_cohorts >= min_cohorts_2plus) |>
    pull(celltype_r2)
  common_ct <- union(ct_1cycle, ct_2plus)
  if (length(common_ct) == 0) return(character(0))

  # Build a combined table: use 1-cycle effect when available, fallback to 2+
  eff_1 <- results_1cycle |>
    filter(celltype_r2 %in% common_ct) |>
    select(celltype_r2, pooled_effect)
  eff_2 <- results_2plus |>
    filter(celltype_r2 %in% common_ct) |>
    select(celltype_r2, pooled_effect)
  combined <- eff_1 |>
    full_join(eff_2, by = "celltype_r2", suffix = c("_1cycle", "_2plus")) |>
    mutate(sort_effect = ifelse(!is.na(pooled_effect_1cycle), pooled_effect_1cycle, pooled_effect_2plus)) |>
    arrange(sort_effect)

  return(unique(combined$celltype_r2))
}

# Helper function to calculate plot dimensions
calc_dimensions <- function(n_celltypes, n_cohorts) {
  height <- 1.9 + n_celltypes * 0.25
  width <- 3 + n_cohorts * 0.33
  return(list(width = width, height = height))
}

# Helper to compute shared xlim across two subgroups for a set of cell types
compute_shared_xlim <- function(results_1, results_2, common_roworder,
                                 min_cohorts_1, min_cohorts_2) {
  ci_1 <- results_1 |> filter(celltype_r2 %in% common_roworder, n_cohorts >= min_cohorts_1)
  ci_2 <- results_2 |> filter(celltype_r2 %in% common_roworder, n_cohorts >= min_cohorts_2)
  ci_min <- min(c(ci_1$ci_lower, ci_2$ci_lower), na.rm = TRUE)
  ci_max <- max(c(ci_1$ci_upper, ci_2$ci_upper), na.rm = TRUE)
  ci_range <- ci_max - ci_min
  c(min(ci_min - 0.05 * ci_range, -0.05),
    max(ci_max + 0.05 * ci_range, 0.05))
}

# Helper to compose a combined side-by-side plot from two subgroup component lists
compose_combined <- function(result_1, result_2, title, w_forest, h_annot = 0.5) {
  if (is.null(result_1$forest) || is.null(result_2$forest)) return(NULL)

  n_cohorts_1 <- result_1$n_cohorts
  n_cohorts_2 <- result_2$n_cohorts
  w_left <- w_forest + n_cohorts_1
  w_right <- w_forest + n_cohorts_2

  # Title strips for each subgroup
  title_1 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = result_1$title_label,
                                  size = 4.5, fontface = "bold") +
    theme_void() + coord_cartesian(clip = "off")
  title_2 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = result_2$title_label,
                                  size = 4.5, fontface = "bold") +
    theme_void() + coord_cartesian(clip = "off")

  design <- "AABB\nCDEF\nGHIJ"

  # Fixed absolute sizes (inches) for title and annotation rows
  h_title_in <- 0.3
  h_annot_in <- h_annot
  n_celltypes <- result_1$n_celltypes
  h_content_in <- n_celltypes * 0.25 + 1.1
  height <- h_title_in + h_annot_in + h_content_in
  width <- 3 + (n_cohorts_1 + n_cohorts_2) * 0.33 + 4

# Helper to build a forest+dotplot+annotation for one subgroup
build_summary_plot <- function(final_results_sub, final_cohort_results_sub, df_sub,
                                common_roworder, title_label,
                                min_cohorts_sub, shared_xlim = NULL) {

  plot_data_raw <- final_results_sub |>
    filter(celltype_r2 %in% common_roworder, n_cohorts >= min_cohorts_sub) |>
    arrange(pooled_effect)

  if (nrow(plot_data_raw) == 0) return(NULL)

  # Use shared xlim if provided, otherwise calculate dynamically
  if (!is.null(shared_xlim)) {
    xlim_min <- shared_xlim[1]
    xlim_max <- shared_xlim[2]
  } else {
    ci_min <- min(plot_data_raw$ci_lower, na.rm = TRUE)
    ci_max <- max(plot_data_raw$ci_upper, na.rm = TRUE)
    ci_range <- ci_max - ci_min
    xlim_min <- ci_min - 0.05 * ci_range
    xlim_max <- ci_max + 0.05 * ci_range
    xlim_min <- min(xlim_min, -0.05)
    xlim_max <- max(xlim_max, 0.05)
  }

  # Dynamic breaks that adapt to the data range, always including 0
  x_breaks <- pretty(c(xlim_min, xlim_max), n = 3)
  if (!0 %in% x_breaks) x_breaks <- sort(unique(c(0, x_breaks)))
  x_breaks <- x_breaks[x_breaks >= xlim_min & x_breaks <= xlim_max]

  plot_data <- plot_data_raw |>
    mutate(
      sig_label = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      ci_lower_clipped = pmax(ci_lower, xlim_min),
      ci_upper_clipped = pmin(ci_upper, xlim_max),
      needs_left_arrow = ci_lower < xlim_min,
      needs_right_arrow = ci_upper > xlim_max,
      celltype_r2 = factor(celltype_r2, levels = common_roworder),
      significant = p_value < 0.05
    )

  forest_plot <- ggplot(plot_data, aes(x = pooled_effect, y = celltype_r2)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(
      aes(xmin = ci_lower_clipped, xmax = ci_upper_clipped),
      height = 0.2, color = "black", linewidth = 0.4
    ) +
    geom_segment(
      data = subset(plot_data, needs_left_arrow),
      aes(x = xlim_min, xend = xlim_min - 0.02, y = celltype_r2, yend = celltype_r2),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black"
    ) +
    geom_segment(
      data = subset(plot_data, needs_right_arrow),
      aes(x = xlim_max, xend = xlim_max + 0.02, y = celltype_r2, yend = celltype_r2),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"), color = "black"
    ) +
    geom_point(
      data = subset(plot_data, significant),
      aes(color = ifelse(pooled_effect > 0, "Treated", "Baseline")),
      size = 3, shape = 16
    ) +
    geom_point(
      data = subset(plot_data, !significant),
      aes(color = ifelse(pooled_effect > 0, "Treated", "Baseline")),
      size = 3, shape = 1, stroke = 1
    ) +
    geom_text(
      aes(x = xlim_max + 0.06, label = sig_label),
      hjust = 0, size = 6, color = "black"
    ) +
    scale_color_manual(
      values = c("Baseline" = color_baseline, "Treated" = color_treated),
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = function(x) sprintf("%.1f", x)
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Pooled Effect Size \n(rank-biserial r)", y = NULL, title = "") +
    coord_cartesian(xlim = c(xlim_min, xlim_max), clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      plot.margin = margin(t = 5, r = 30, b = 25, l = 5),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 11, color = "black", margin = margin(t = 8)),
      axis.text.y = element_text(size = 11, color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}
  # Dotplot
  cohorts_in_data <- final_cohort_results_sub |>
    filter(celltype_r2 %in% common_roworder) |>
    pull(cohort) |>
    unique()
  cohort_order_filtered <- master_cohort_order[master_cohort_order %in% cohorts_in_data]

  cohort_treatment <- data.frame(
    cohort = cohort_order_filtered,
    treatment = df_sub$treatment[match(cohort_order_filtered, df_sub$cohort)]
  )

  dotplot_data <- final_cohort_results_sub |>
    filter(celltype_r2 %in% common_roworder) |>
    left_join(cohort_treatment, by = "cohort") |>
    mutate(
      celltype_r2 = factor(celltype_r2, levels = common_roworder),
      cohort = factor(cohort, levels = cohort_order_filtered),
      effect_clipped = pmax(-effect_clip, pmin(effect_clip, effect_size)),
      weight = pmax(1, pmin(6, 1 / ifelse(is.na(effect_se), 1e+10, effect_se))),
      significant = p_value < 0.05
    )

  dotplot <- ggplot(dotplot_data, aes(x = cohort, y = celltype_r2)) +
    geom_point(aes(color = effect_clipped, size = weight), na.rm = TRUE) +
    geom_point(data = subset(dotplot_data, significant),
               aes(size = weight), shape = 21, color = "black", fill = NA, stroke = 0.5) +
    scale_color_gradient2(
      low = color_baseline, mid = "white", high = color_treated,
      midpoint = 0, limits = c(-effect_clip, effect_clip),
      na.value = "lightgray", name = "Effect size"
    ) +
    scale_size_continuous(range = c(0, 5), name = "Confidence",
                          limits = c(0, 6), breaks = c(2, 4, 6), labels = c("2", "4", "6")) +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9, color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(t = 5, r = 0, b = 25, l = 0)
    )

  # Annotation bars
  annotation_data <- data.frame(
    cohort = factor(cohort_order_filtered, levels = cohort_order_filtered),
    treatment = df_sub$treatment[match(cohort_order_filtered, df_sub$cohort)],
    subtype = df_sub$subtype[match(cohort_order_filtered, df_sub$cohort)]
  )

  combined_annotation <- ggplot(annotation_data, aes(x = cohort)) +
    geom_tile(aes(y = 1, fill = subtype), color = "white", linewidth = 0.5) +
    scale_fill_manual(values = cancerytpe_colors, name = "Cancer type") +
    new_scale_fill() +
    geom_tile(aes(y = 2, fill = treatment), color = "white", linewidth = 0.5) +
    scale_fill_manual(values = treatment_colors, name = "Treatment") +
    annotate("text", x = 0.35, y = 1, label = "Cancer type", hjust = 1, size = 2.8) +
    annotate("text", x = 0.35, y = 2, label = "Treatment", hjust = 1, size = 2.8) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_line(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 55),
      legend.position = "right",
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  return(list(forest = forest_plot, dotplot = dotplot, annotation = combined_annotation,
              title_label = title_label,
              n_celltypes = length(common_roworder), n_cohorts = length(cohort_order_filtered)))
}

# --- Load pooled (overall) meta-analysis results ---
pooled_meta <- read.csv('tables/wicox_r2_per_cohort_tx_status_meta.csv')
min_cohorts_pooled <- 5

# --- Significant cell types from either cycle group OR pooled ---
sig_1cycle <- final_results_1cycle |>
  filter(n_cohorts >= min_cohorts_1cycle, p_value < 0.05) |> pull(celltype_r2)
sig_2plus <- final_results_2plus |>
  filter(n_cohorts >= min_cohorts_2plus, p_value < 0.05) |> pull(celltype_r2)
sig_pooled <- pooled_meta |>
  filter(n_cohorts >= min_cohorts_pooled, p_value < 0.05) |> pull(celltype_r2)
ov2_sig_ct <- union(sig_1cycle, sig_2plus) |> union(sig_pooled)

# --- Merge all three groups ---
merged_ov2 <- pooled_meta |>
  filter(celltype_r2 %in% ov2_sig_ct, n_cohorts >= min_cohorts_pooled) |>
  select(celltype_r2, pooled_effect, ci_lower, ci_upper, p_value) |>
  mutate(cycle = "Pooled") |>
  bind_rows(
    final_results_1cycle |>
      filter(celltype_r2 %in% ov2_sig_ct, n_cohorts >= min_cohorts_1cycle) |>
      select(celltype_r2, pooled_effect, ci_lower, ci_upper, p_value) |>
      mutate(cycle = "1 cycle"),
    final_results_2plus |>
      filter(celltype_r2 %in% ov2_sig_ct, n_cohorts >= min_cohorts_2plus) |>
      select(celltype_r2, pooled_effect, ci_lower, ci_upper, p_value) |>
      mutate(cycle = "2+ cycles")
  )

# --- Sort by pooled effect (fallback to 1c, then 2+c) ---
ov2_sort <- pooled_meta |>
  filter(celltype_r2 %in% ov2_sig_ct) |>
  select(celltype_r2, eff_pooled = pooled_effect) |>
  full_join(
    final_results_1cycle |> filter(celltype_r2 %in% ov2_sig_ct) |>
      select(celltype_r2, eff_1c = pooled_effect), by = "celltype_r2") |>
  full_join(
    final_results_2plus |> filter(celltype_r2 %in% ov2_sig_ct) |>
      select(celltype_r2, eff_2p = pooled_effect), by = "celltype_r2") |>
  mutate(sort_eff = coalesce(eff_pooled, eff_1c, eff_2p)) |>
  arrange(sort_eff)
ov2_roworder <- ov2_sort$celltype_r2

# --- Prepare merged data ---
# Dodge order (bottom to top): 2+ cycles, Pooled, 1 cycle
merged_ov2 <- merged_ov2 |>
  mutate(
    celltype_r2 = factor(celltype_r2, levels = ov2_roworder),
    cycle = factor(cycle, levels = c("2+ cycles", "Pooled", "1 cycle")),
    significant = p_value < 0.05,
    direction = ifelse(pooled_effect > 0, "Treated > BL", "BL > Treated"),
    pt_alpha = ifelse(significant, 1, 0.35)
  )

# --- xlim ---
ov2_ci_min <- min(merged_ov2$ci_lower, na.rm = TRUE)
ov2_ci_max <- max(merged_ov2$ci_upper, na.rm = TRUE)
ov2_ci_rng <- ov2_ci_max - ov2_ci_min
ov2_xlim <- c(min(ov2_ci_min - 0.05 * ov2_ci_rng, -0.1),
              max(ov2_ci_max + 0.05 * ov2_ci_rng, 0.1))

merged_ov2 <- merged_ov2 |>
  mutate(
    ci_lower_clip = pmax(ci_lower, ov2_xlim[1]),
    ci_upper_clip = pmin(ci_upper, ov2_xlim[2])
  )

# Simple shapes: triangle=1c, square=2+c, circle=Pooled
# Significance shown via alpha (opaque = p<0.05, transparent = n.s.)
ov2_shape_vals <- c("1 cycle" = 17, "2+ cycles" = 15, "Pooled" = 16)
ov2_dodge <- 0.55

# --- Striped background data (ggforestplot style) ---
ov2_n_ct <- length(ov2_roworder)
stripe_df <- data.frame(
  ymin = seq(0.5, ov2_n_ct + 0.5, by = 2),
  ymax = pmin(seq(1.5, ov2_n_ct + 1.5, by = 2), ov2_n_ct + 0.5)
)

# --- Forest plot with striped background ---
p_ov2_forest <- ggplot(merged_ov2, aes(x = pooled_effect, y = celltype_r2,
                                        group = cycle)) +
  # Striped background
  geom_rect(data = stripe_df,
            aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
            fill = "#ededed", color = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # Black CI segments (Pooled thicker than cycle subgroups)
  geom_errorbarh(
    aes(xmin = ci_lower_clip, xmax = ci_upper_clip, linewidth = cycle),
    height = 0.08, color = "black",
    position = position_dodge(width = ov2_dodge)
  ) +
  scale_linewidth_manual(
    values = c("1 cycle" = 0.2, "2+ cycles" = 0.2, "Pooled" = 0.4),
    guide = "none"
  ) +
  # Points: shape = cycle, color = direction, alpha = significance
  geom_point(aes(shape = cycle, color = direction, alpha = pt_alpha),
             size = 2.5,
             position = position_dodge(width = ov2_dodge)) +
  scale_color_manual(
    values = c("BL > Treated" = color_baseline, "Treated > BL" = color_treated),
    name = "Direction"
  ) +
  scale_shape_manual(
    values = ov2_shape_vals, name = "Analysis",
    guide = guide_legend(reverse = TRUE, override.aes = list(color = "black", alpha = 1, size = 2.5))
  ) +
  scale_alpha_identity() +
  scale_x_continuous(
    breaks = pretty(ov2_xlim, n = 5),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_y_discrete(drop = FALSE) +
  coord_cartesian(xlim = ov2_xlim, clip = "off") +
  labs(x = "Pooled effect size\n(rank-biserial r)", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.border     = element_rect(fill = NA, color = "grey70"),
    axis.text.y      = element_text(size = 11, color = "black"),
    axis.text.x      = element_text(size = 10, color = "black"),
    axis.title.x     = element_text(size = 11, color = "black", margin = margin(t = 8)),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 25, l = 5)
  )

# --- Combined dotplot: all cohorts ordered by cancer type, then cycle ---
combined_cohort_res <- bind_rows(
  final_cohort_results_1cycle |> mutate(cycle_group = "1 cycle"),
  final_cohort_results_2plus |> mutate(cycle_group = "2+ cycles")
)

# Cohort order: cancer type first, then cycle within cancer type
combined_cohort_order <- c(
  "SKCM_Mono", "SKCM_Combo",
  "BRCA_Bassez1", "BRCA_Bassez2",
  "TNBC_Bassez1", "TNBC_Bassez2", "TNBC_Shiao", "TNBC_Zhang",
  "HNSCC_Franken_Mono", "HNSCC_Franken_Combo", "HNSCC_Luoma", "HNSCC_vanderLeun",
  "CRC_Chen", "CRC_Chen_2cycles", "CRC_Chen_3cycles+", "CRC_Li",
  "NSCLC_Yan", "HCC_Guo", "PCa_Hawley"
)
cohorts_present <- combined_cohort_res |> pull(cohort) |> unique()
combined_cohort_order <- combined_cohort_order[combined_cohort_order %in% cohorts_present]

# Get cohort metadata from both subgroups
cohort_meta_ov2 <- bind_rows(
  df_1cycle |> select(cohort, subtype, treatment) |> distinct(),
  df_2plus |> select(cohort, subtype, treatment) |> distinct()
) |> distinct(cohort, .keep_all = TRUE)

# Add cycle_group
cycle_grp_map <- combined_cohort_res |> select(cohort, cycle_group) |> distinct()
cohort_meta_ov2 <- cohort_meta_ov2 |> left_join(cycle_grp_map, by = "cohort")

# Prepare dotplot data
dot_data_ov2 <- combined_cohort_res |>
  filter(celltype_r2 %in% ov2_roworder) |>
  left_join(cohort_meta_ov2 |> select(cohort, treatment), by = "cohort") |>
  mutate(
    celltype_r2 = factor(celltype_r2, levels = ov2_roworder),
    cohort = factor(cohort, levels = combined_cohort_order),
    effect_clipped = pmax(-effect_clip, pmin(effect_clip, effect_size)),
    weight = pmax(1, pmin(6, 1 / ifelse(is.na(effect_se), 1e10, effect_se))),
    significant = p_value < 0.05
  )

p_ov2_dotplot <- ggplot(dot_data_ov2, aes(x = cohort, y = celltype_r2)) +
  geom_point(aes(color = effect_clipped, size = weight), na.rm = TRUE) +
  geom_point(data = subset(dot_data_ov2, significant),
             aes(size = weight), shape = 21, color = "black", fill = NA, stroke = 0.5) +
  scale_color_gradient2(
    low = color_baseline, mid = "white", high = color_treated,
    midpoint = 0, limits = c(-effect_clip, effect_clip),
    na.value = "lightgray", name = "Effect size"
  ) +
  scale_size_continuous(range = c(0, 5), name = "Confidence",
                        limits = c(0, 6), breaks = c(2, 4, 6)) +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_line(color = "black"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 0, b = 25, l = 0)
  )

# --- Combined annotation: cancer type + treatment + cycle ---
ov2_all_subtypes <- sort(unique(c(df_1cycle$subtype, df_2plus$subtype)))
ov2_all_treatments <- names(treatment_colors)[names(treatment_colors) %in%
                        unique(c(df_1cycle$treatment, df_2plus$treatment))]
hokusai2_cols <- MetBrewer::met.brewer('Hokusai2', 6)
cycle_colors <- c("1 cycle" = hokusai2_cols[1], "2+ cycles" = hokusai2_cols[4])

ann_data_ov2 <- data.frame(
  cohort = factor(combined_cohort_order, levels = combined_cohort_order),
  subtype = cohort_meta_ov2$subtype[match(combined_cohort_order, cohort_meta_ov2$cohort)],
  treatment = cohort_meta_ov2$treatment[match(combined_cohort_order, cohort_meta_ov2$cohort)],
  cycle_group = cohort_meta_ov2$cycle_group[match(combined_cohort_order, cohort_meta_ov2$cohort)]
)

p_ov2_annot <- ggplot(ann_data_ov2, aes(x = cohort)) +
  geom_tile(aes(y = 1, fill = subtype), color = "white", linewidth = 0.5) +
  scale_fill_manual(values = cancerytpe_colors, name = "Cancer type",
                    limits = ov2_all_subtypes) +
  new_scale_fill() +
  geom_tile(aes(y = 2, fill = treatment), color = "white", linewidth = 0.5) +
  scale_fill_manual(values = treatment_colors, name = "Treatment",
                    limits = ov2_all_treatments) +
  new_scale_fill() +
  geom_tile(aes(y = 3, fill = cycle_group), color = "white", linewidth = 0.5) +
  scale_fill_manual(values = cycle_colors, name = "Tx cycle") +
  annotate("text", x = 0.35, y = 1, label = "Cancer type", hjust = 1, size = 2.5) +
  annotate("text", x = 0.35, y = 2, label = "Treatment", hjust = 1, size = 2.5) +
  annotate("text", x = 0.35, y = 3, label = "Treatment cycle", hjust = 1, size = 2.5) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(breaks = NULL) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 11) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), panel.border = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 55)
  )

# --- Compose with patchwork ---
# Layout:  row 1: [empty] [annotation]
#          row 2: [forest] [dotplot]
ov2_nc <- length(combined_cohort_order)
ov2_w_forest <- 7.5
ov2_h_annot <- 0.6
ov2_h_content <- ov2_n_ct * 0.25 + 1.1

ov2_design <- "
#A
BC
"

p_ov2_combined <- p_ov2_annot +
  p_ov2_forest + p_ov2_dotplot +
  plot_layout(
    design = ov2_design,
    heights = c(ov2_h_annot, ov2_h_content),
    widths = c(ov2_w_forest, ov2_nc * 0.33),
    guides = "collect"
  ) &
  theme(legend.position = "right")


ggsave('figures/Fig2/Temporal_r2_meta_overlay_summary.pdf',
       p_ov2_combined, width = 12, height = 9)


# per-patient fold-change dynamics
temporal_cohorts <- c('SKCM_this study', 'SKCM_Plozniak',
                      'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao', 'TNBC_Zhang',
                      'HNSCC_Franken', 'HNSCC_Luoma', 'HNSCC_vanderLeun',
                      'CRC_Chen', 'CRC_Li','NSCLC_Yan', 'HCC_Guo')

unwanted_celltypes <- c('Melanocytes(CNA-)', 'Melanocytes(CNA+)',
                         'Epithelial(CNA+)', 'Epithelial(CNA-)', 'Malignant(CNA+)')
unmatched_pt <- c("SKCM_this study_Patient2", "SKCM_this study_Patient3")

df <- metadata |>
  filter(
    cohort %in% temporal_cohorts,
    count_immune >= 50,
    non_malignant_count >= 100,
    !patient %in% unmatched_pt,
    tx_status %in% c("Baseline", "Treated"),
    subset %in% c('TME', 'CD45+sorted'),
    !celltype_r2 %in% unwanted_celltypes,
    tx_cycle != "Not specified"
  ) |>
  distinct(sample, celltype_r2, freq_r2_comp, .keep_all = TRUE) |>
  mutate(
    cohort = case_when(
      cohort == "BRCA_Bassez1" & subtype == "TNBC" ~ "TNBC_Bassez1",
      cohort == "BRCA_Bassez2" & subtype == "TNBC" ~ "TNBC_Bassez2",
      TRUE ~ cohort
    ),
    cancer_type = subtype,
    proportion = freq_r2_comp,
    tx_status = factor(tx_status, levels = c("Baseline", "Treated"))
  )

# Keep only matched patients (both Baseline and Treated)
matched_patients <- df |>
  group_by(cancer_type, cohort, patient) |>
  summarise(has_both = n_distinct(tx_status) == 2, .groups = "drop") |>
  filter(has_both) |>
  pull(patient) |>
  unique()

df_matched <- df |> filter(patient %in% matched_patients)


cancers <- c("SKCM", "TNBC", "HNSCC", "CRC", "NSCLC", "HCC")
pseudocount <- 0.001
min_n <- 3  # minimum samples per group to display

# Colors consistent with r2_overview_split.R (Austria palette, order from metadata$subtype)
cancer_cols <- c("SKCM" = "#A40000", "BCC" = "#452053", "TNBC" = "#0E4A63", "BRCA" = "#007E2F",
                 "HNSCC" = "#A9B21B", "CRC" = "#E7A83C", "NSCLC" = "#B86092", "PCa" = "#89325A",
                 "RCC" = "#4C4F61", "HCC" = "#00B7A7")


# compute per-patient fold change and summarise

compute_fc <- function(df_sub) {
  df_wide <- df_sub %>%
    pivot_wider(names_from = tx_status, values_from = proportion) %>%
    filter(!is.na(Baseline), !is.na(Treated)) %>%
    mutate(
      fc = (Treated + pseudocount) / (Baseline + pseudocount),
      timepoint = case_when(
        cancer_type == "TNBC" & grepl("Bassez", cohort) ~ "1 cycle\n(early)",
        cancer_type == "SKCM"                            ~ "1 cycle\n(early)",
        cohort == "CRC_Chen" & tx_cycle == "3+ cycles"  ~ "3 cycles",
        cohort == "CRC_Li"                               ~ "4 cycles",
        cohort == "NSCLC_Yan"                              ~ "3 cycles",
        cohort == "HCC_Guo"                                ~ "3 cycles",
        TRUE ~ as.character(tx_cycle)
      ),
      timepoint = factor(timepoint,
        levels = c("1 cycle\n(early)", "1 cycle", "2 cycles", "3 cycles", "4 cycles"))
    )

  summary_df <- df_wide %>%
    group_by(cancer_type, timepoint) %>%
    summarise(
      median_fc = median(fc, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n > min_n)

  # Add baseline rows
  keep_cancers <- unique(summary_df$cancer_type)
  baseline_rows <- data.frame(
    cancer_type = keep_cancers,
    timepoint = factor("Baseline", levels = c("Baseline", levels(summary_df$timepoint))),
    median_fc = 1, n = NA
  )
  summary_df$timepoint <- factor(summary_df$timepoint,
    levels = c("Baseline", levels(summary_df$timepoint)))
  summary_df <- bind_rows(baseline_rows, summary_df)

  return(summary_df)
}

# plot fold change dynamics

plot_fc <- function(summary_df, title) {
  # Dynamic y-axis
  ymin <- min(summary_df$median_fc, na.rm = TRUE) * 0.9
  ymax <- max(summary_df$median_fc, na.rm = TRUE) * 1.1
  breaks <- c(0.125, 0.25, 0.5, 1, 2, 4)
  breaks <- breaks[breaks >= ymin & breaks <= ymax]
  if (!1 %in% breaks) breaks <- sort(c(breaks, 1))

  label_df <- summary_df %>% filter(!is.na(n))

  p <- ggplot(summary_df, aes(x = timepoint, y = median_fc,
                               color = cancer_type, group = cancer_type)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_text_repel(data = label_df,
      aes(label = paste0("n=", n)),
      size = 4, show.legend = FALSE, color = 'black',
      box.padding = 0.4, point.padding = 0.3,
      min.segment.length = 0.2,
      segment.color = "grey60", segment.size = 0.3,
      direction = "both", max.overlaps = 20) +
    scale_color_manual(values = cancer_cols) +
    scale_y_continuous(trans = "log2",
      limits = c(ymin, ymax),
      breaks = breaks,
      labels = as.character(breaks)) +
    labs(
      title = title,
      x = "",
      y = "Log2(fold change)",
      color = "Cancer Type"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 13, hjust = 0.5),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )

  return(p)
}


ct <- "CD4_Tfh"

df_sub <- df_matched %>%
  filter(cancer_type %in% cancers, celltype_r2 == ct) %>%
  select(cancer_type, cohort, patient, tx_status, tx_cycle, proportion)

summary_df <- compute_fc(df_sub)
p1 <- plot_fc(summary_df, ct);p1


cts <- c("CD4_T-ISG", "CD8_T-ISG")

df_sub <- df_matched %>%
  filter(cancer_type %in% cancers, celltype_r2 %in% cts) %>%
  filter(!(cohort == "CRC_Chen" & tx_cycle == "2 cycles")) %>%
  group_by(cancer_type, cohort, patient, tx_status, tx_cycle) %>%
  summarise(proportion = sum(proportion, na.rm = TRUE), .groups = "drop")

summary_df <- compute_fc(df_sub)
p2 <- plot_fc(summary_df, "T-ISG");p2


ct3 <- "Macro-ISG"

df_sub3 <- df_matched %>%
  filter(cancer_type %in% cancers, celltype_r2 == ct3) %>%
  select(cancer_type, cohort, patient, tx_status, tx_cycle, proportion)

summary_df3 <- compute_fc(df_sub3)
p3 <- plot_fc(summary_df3, ct3)


ct4 <- "CD4_Treg"

df_sub4 <- df_matched %>%
  filter(cancer_type %in% cancers, celltype_r2 == ct4) %>%
  select(cancer_type, cohort, patient, tx_status, tx_cycle, proportion)

summary_df4 <- compute_fc(df_sub4)
p4 <- plot_fc(summary_df4, ct4)

p1 + p2 + p3 + p4 + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "bottom")
ggsave('figures/dynamics_fc.pdf', height = 5, width = 20)


