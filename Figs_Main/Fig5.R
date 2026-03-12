pkgs <- c('qs2', 'tidyr', 'dplyr', 'stringr', 'tibble', 'ggplot2',
          'GSVA', 'pROC', 'survival', 'ranger',
          'foreach', 'doParallel', 'metafor', 'scales', 'patchwork')
invisible(lapply(pkgs, function(x) {
  require(package = x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
}))
options(warn = -1)

pb_mat <- qs_read("data/pseudobulk_TIME.qs2")
sample_info <- read.csv("tables/sample_info_time_updated.csv")

paired_bl <- sample_info |>
  dplyr::filter(paired == "Yes", tx_status == "Baseline") |>
  dplyr::select(patient, sample, group, response, cohort)
paired_tx <- sample_info |>
  dplyr::filter(paired == "Yes", tx_status == "Treated") |>
  dplyr::select(patient, group) |>
  dplyr::rename(group_treated = group)

transitions <- dplyr::left_join(paired_bl, paired_tx, by = "patient")
transitions$sample_hyphen <- str_replace_all(transitions$sample, "_", "-")

# Transition labels per TIME group
# I:   favorable = I->B/Q, unfavorable = I->Mye (I->I excluded)
# Mye: favorable = Mye->Others, unfavorable = Mye->Mye
transitions$label <- dplyr::case_when(
  transitions$group == "Immune Inflamed" &
    transitions$group_treated %in% c("Immune Quiescent", "B Cell-enriched") ~ "favorable",
  transitions$group == "Immune Inflamed" &
    transitions$group_treated == "Myeloid-enriched" ~ "unfavorable",
  transitions$group == "Myeloid-enriched" &
    transitions$group_treated != "Myeloid-enriched" ~ "favorable",
  transitions$group == "Myeloid-enriched" &
    transitions$group_treated == "Myeloid-enriched" ~ "unfavorable",
  TRUE ~ "other"
)

i_trans <- transitions |> dplyr::filter(grepl("Inflamed", group), label != "other")
mye_trans <- transitions |> dplyr::filter(grepl("Myeloid", group), label != "other")

# --- TIME-I ---
i_labeled <- transitions |>
  dplyr::filter(grepl("Inflamed", group),
                label %in% c("favorable", "unfavorable"),
                sample_hyphen %in% colnames(pb_mat))

pb_expr_I <- pb_mat[, i_labeled$sample_hyphen, drop = FALSE]
label_I <- i_labeled$label

i_resp <- i_labeled |>
  dplyr::filter(response %in% c("R", "NR")) |>
  dplyr::group_by(label) |>
  dplyr::summarise(n = dplyr::n(), n_R = sum(response == "R"),
                   pct_R = round(100 * n_R / n, 1), .groups = "drop")

# --- TIME-Mye ---
mye_labeled <- transitions |>
  dplyr::filter(grepl("Myeloid", group),
                label %in% c("favorable", "unfavorable"),
                sample_hyphen %in% colnames(pb_mat))

pb_expr_Mye <- pb_mat[, mye_labeled$sample_hyphen, drop = FALSE]
label_Mye <- mye_labeled$label

mye_resp <- mye_labeled |>
  dplyr::filter(response %in% c("R", "NR")) |>
  dplyr::group_by(label) |>
  dplyr::summarise(n = dplyr::n(), n_R = sum(response == "R"),
                   pct_R = round(100 * n_R / n, 1), .groups = "drop")

baseline_clin <- clin_info |>
  dplyr::filter(tx_status == "Baseline", response %in% c("R", "NR"))

# Gene pools (kept for reference)
gene_pools <- list()
gene_pools_filtered <- list()

# RF stability selection
n_cores <- min(parallel::detectCores() - 1, 40)
n_hvg   <- 4000
n_iter  <- 200
top_n   <- 100

run_stability_rf <- function(seed, x, y, wts, top_n, ntree = 1000) {
  set.seed(seed)
  mtry_val <- max(1, floor(sqrt(ncol(x))))
  train_df <- as.data.frame(x, check.names = FALSE)
  train_df$transition <- y
  rf_fit <- tryCatch(
    ranger::ranger(transition ~ ., data = train_df, num.trees = ntree,
                   mtry = mtry_val, importance = "permutation",
                   case.weights = wts, num.threads = 1),
    error = function(e) NULL)
  if (is.null(rf_fit)) return(NULL)
  imp <- sort(ranger::importance(rf_fit), decreasing = TRUE)
  list(top_genes = names(imp)[seq_len(min(top_n, length(imp)))],
       all_importance = imp, oob_error = rf_fit$prediction.error)
}

prepare_pool_matrix <- function(genes, pb_expr) {
  x <- as.matrix(t(pb_expr[genes, , drop = FALSE]))
  for (j in seq_len(ncol(x))) {
    nas <- is.na(x[, j]) | is.nan(x[, j]) | is.infinite(x[, j])
    if (any(nas)) x[nas, j] <- median(x[!nas, j], na.rm = TRUE)
  }
  v <- apply(x, 2, var, na.rm = TRUE)
  keep <- !is.na(v) & v > 0
  x <- x[, keep, drop = FALSE]
  orig <- colnames(x)
  safe_names <- make.names(orig, unique = TRUE)
  colnames(x) <- safe_names
  list(x = x, name_map = setNames(orig, safe_names))
}

run_pool_stability <- function(pool_name, x, name_map, y, wts) {
  cat(sprintf("  Running %d RF iterations on %d cores...\n", n_iter, n_cores))
  pool_runs <- parallel::mclapply(seq_len(n_iter), function(i) {
    tryCatch(
      run_stability_rf(seed = i * 7 + abs(digest::digest2int(pool_name, seed = 42)),
                       x = x, y = y, wts = wts, top_n = top_n),
      error = function(e) NULL)
  }, mc.cores = n_cores)
  pool_runs <- pool_runs[!sapply(pool_runs, is.null)]
  if (length(pool_runs) < 50) {
    cat("  WARNING: only", length(pool_runs), "successful runs\n")
    return(NULL)
  }
  oob_errors <- sapply(pool_runs, function(r) r$oob_error)
  cat("  Successful:", length(pool_runs), "/", n_iter, "\n")
  cat("  Mean OOB:", round(mean(oob_errors), 4), "\n")
  gene_counts <- table(unlist(lapply(pool_runs, function(r) r$top_genes)))
  gene_freq <- sort(gene_counts / length(pool_runs), decreasing = TRUE)
  all_genes_u <- unique(unlist(lapply(pool_runs, function(r) names(r$all_importance))))
  mean_imp <- setNames(rep(0, length(all_genes_u)), all_genes_u)
  for (r in pool_runs) {
    imp <- r$all_importance
    mean_imp[names(imp)] <- mean_imp[names(imp)] + imp
  }
  mean_imp <- sort(mean_imp / length(pool_runs), decreasing = TRUE)
  list(n_runs = length(pool_runs), mean_oob = mean(oob_errors),
       gene_freq = gene_freq, mean_imp = mean_imp, name_map = name_map)
}

find_elbow <- function(freq_vec, min_genes = 15, max_genes = 80) {
  n <- length(freq_vec)
  if (n < min_genes) return(n)
  idx_range <- min_genes:min(max_genes, n)
  vals <- as.numeric(freq_vec[idx_range])
  x <- seq(0, 1, length.out = length(idx_range))
  y <- (vals - min(vals)) / (max(vals) - min(vals) + 1e-10)
  x1 <- x[1]; y1 <- y[1]; x2 <- x[length(x)]; y2 <- y[length(y)]
  dists <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) /
    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  idx_range[which.max(dists)]
}

#' Extract directional signature: keep ALL elbow genes with up/down labels
extract_directional_signature <- function(stab, elbow_n, x_pool, fav_idx, unfav_idx) {
  gene_freq <- stab$gene_freq
  mean_imp  <- stab$mean_imp
  name_map  <- stab$name_map
  candidates <- names(gene_freq)[1:min(elbow_n, length(gene_freq))]
  sig_genes <- unname(name_map[candidates])
  sig_importance <- unname(mean_imp[candidates]); names(sig_importance) <- sig_genes
  sig_freq <- as.numeric(gene_freq[candidates]); names(sig_freq) <- sig_genes
  
  gene_dir <- vapply(candidates, function(g) {
    mf <- mean(x_pool[fav_idx, g], na.rm = TRUE)
    mu <- mean(x_pool[unfav_idx, g], na.rm = TRUE)
    if (mf >= mu) "up" else "down"
  }, character(1))
  names(gene_dir) <- sig_genes
  
  up   <- names(gene_dir[gene_dir == "up"])
  down <- names(gene_dir[gene_dir == "down"])
  cat(sprintf("    Elbow: %d genes (%d up, %d down)\n", length(sig_genes), length(up), length(down)))
  
  list(size = length(sig_genes), genes = sig_genes,
       up_genes = up, down_genes = down,
       importance = sig_importance, frequency = sig_freq, direction = gene_dir)
}

# Run RF for TIME-I and TIME-Mye 
run_group <- function(pb_expr, labels, group_name) {
  cat(sprintf("\n########## %s (HVG = %d) ##########\n", group_name, n_hvg))
  cat("  Samples:", ncol(pb_expr), "\n")
  cat("  favorable:", sum(labels == "favorable"),
      " unfavorable:", sum(labels == "unfavorable"), "\n")
  
  y <- factor(labels, levels = c("unfavorable", "favorable"))
  wts <- ifelse(y == "unfavorable", sum(y == "favorable") / sum(y == "unfavorable"), 1)
  fav_idx   <- which(y == "favorable")
  unfav_idx <- which(y == "unfavorable")
  
  x_all <- as.matrix(t(pb_expr))
  gene_vars <- apply(x_all, 2, var, na.rm = TRUE)
  gene_vars <- gene_vars[!is.na(gene_vars) & gene_vars > 0]
  top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(n_hvg, length(gene_vars))]
  cat("  Genes after variance filter:", length(top_var_genes), "\n")
  
  prep <- prepare_pool_matrix(top_var_genes, pb_expr)
  stab <- run_pool_stability(paste0(group_name, "_HVG_", n_hvg), prep$x, prep$name_map, y, wts)
  if (is.null(stab)) stop(paste("RF failed for", group_name))
  
  elbow_n <- find_elbow(stab$gene_freq)
  cat(sprintf("  Elbow at: %d genes (freq threshold: %.3f)\n",
              elbow_n, as.numeric(stab$gene_freq[elbow_n])))
  
  sig <- extract_directional_signature(stab, elbow_n, prep$x, fav_idx, unfav_idx)
  list(stab = stab, prep = prep, elbow_n = elbow_n, sig = sig)
}

I_result   <- run_group(pb_expr_I,   label_I,   "I")
Mye_result <- run_group(pb_expr_Mye, label_Mye, "Mye")

i_sig <- I_result$sig
m_sig <- Mye_result$sig
all_genes <- union(i_sig$genes, m_sig$genes)

# Resolve direction: use importance-weighted direction when both have the gene
gene_direction <- vapply(all_genes, function(g) {
  in_i <- g %in% i_sig$genes; in_m <- g %in% m_sig$genes
  if (in_i && in_m) {
    if (i_sig$importance[g] >= m_sig$importance[g]) i_sig$direction[g] else m_sig$direction[g]
  } else if (in_i) { i_sig$direction[g] } else { m_sig$direction[g] }
}, character(1))

gene_importance <- vapply(all_genes, function(g) {
  imp_i <- if (g %in% names(i_sig$importance)) unname(i_sig$importance[g]) else 0
  imp_m <- if (g %in% names(m_sig$importance)) unname(m_sig$importance[g]) else 0
  max(imp_i, imp_m)
}, numeric(1))

gene_source <- vapply(all_genes, function(g) {
  in_i <- g %in% i_sig$genes; in_m <- g %in% m_sig$genes
  if (in_i && in_m) "both" else if (in_i) "I" else "Mye"
}, character(1))

up_genes   <- names(gene_direction[gene_direction == "up"])
down_genes <- names(gene_direction[gene_direction == "down"])

# Elbow diagnostic plot
pdf(file.path(fig_dir, "elbow_analysis.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

freq_I <- as.numeric(I_result$stab$gene_freq)
plot(1:min(100, length(freq_I)), freq_I[1:min(100, length(freq_I))],
     type = "l", lwd = 2, col = "steelblue",
     xlab = "Gene rank", ylab = "Selection frequency",
     main = sprintf("TIME-I (elbow=%d)", I_result$elbow_n))
abline(v = I_result$elbow_n, col = "red", lty = 2, lwd = 2)
grid()

freq_M <- as.numeric(Mye_result$stab$gene_freq)
plot(1:min(100, length(freq_M)), freq_M[1:min(100, length(freq_M))],
     type = "l", lwd = 2, col = "darkorange",
     xlab = "Gene rank", ylab = "Selection frequency",
     main = sprintf("TIME-Mye (elbow=%d)", Mye_result$elbow_n))
abline(v = Mye_result$elbow_n, col = "red", lty = 2, lwd = 2)
grid()

dev.off()

# Gene table
sig_df <- data.frame(
  gene = all_genes,
  importance = gene_importance[all_genes],
  direction = gene_direction[all_genes],
  source = gene_source[all_genes],
  stringsAsFactors = FALSE
) |> arrange(desc(importance))
write.csv(sig_df, file.path(tab_dir, "transition_signature_genes.csv"), row.names = FALSE)

# Heatmap of all transition signature genes (up + down)
all_genes <- transition_sig$genes
all_genes_present <- intersect(all_genes, rownames(pb_mat))

# Order: direction (up first, then down), within each by source (I, Mye), then by importance
gene_order_df <- data.frame(
  gene      = all_genes_present,
  direction = transition_sig$direction[all_genes_present],
  source    = transition_sig$source[all_genes_present],
  imp       = transition_sig$importance[all_genes_present],
  stringsAsFactors = FALSE
) |> arrange(source, direction, desc(imp))
genes_ordered <- gene_order_df$gene

# Order samples: I-favorable, I-unfavorable, Mye-favorable, Mye-unfavorable
plot_samples <- plot_samples |> arrange(time_group, label, cohort)

expr_mat <- as.matrix(pb_mat[genes_ordered, plot_samples$sample_hyphen, drop = FALSE])
if (max(expr_mat, na.rm = TRUE) > 50) {
  cat("  Applying log2(x+1)\n")
  expr_mat <- log2(expr_mat + 1)
}

# Row-scale, cap at +/-2.5
expr_scaled <- t(scale(t(expr_mat)))
expr_scaled[is.na(expr_scaled)] <- 0
cap <- 2.5
expr_scaled[expr_scaled >  cap] <-  cap
expr_scaled[expr_scaled < -cap] <- -cap

# Column annotations
col_anno <- HeatmapAnnotation(
  `Subtype` = plot_samples$time_group,
  Transition   = plot_samples$label,
  Response     = plot_samples$response,
  col = list(
    `Subtype` = c('TIME-I' = "#208cc0", 'TIME-Mye' = "#cf5e4e"),
    Transition   = c("favorable" = "#009E73", "unfavorable" = "#E69F00"),
    Response     = c("R" = "#CC0C00FF", "NR" = "#5C88DAFF",
                     "NE" = "grey70", "SD" = "grey85")
  ),
  annotation_name_side = "right",
  gap = unit(1, "mm")
)

# Left annotation: RF importance barplot (bars extend leftward, away from heatmap)
row_imp <- gene_order_df$imp[match(genes_ordered, gene_order_df$gene)]
left_anno <- rowAnnotation(
  Importance = anno_barplot(
    row_imp,
    width = unit(2, "cm"),
    gp = gpar(fill = "grey65", col = NA),
    baseline = max(row_imp),
    axis_param = list(gp = gpar(fontsize = 6))
  ),
  annotation_name_rot = 90
)

# Right annotation: source + direction
row_source_raw <- transition_sig$source[genes_ordered]
row_source <- ifelse(row_source_raw == "I", "TIME-I",
                     ifelse(row_source_raw == "Mye", "TIME-Mye", row_source_raw))
row_dir    <- transition_sig$direction[genes_ordered]
right_anno <- rowAnnotation(
  Subtype = row_source,
  col = list(
    Subtype = c('TIME-I' = "#208cc0", 'TIME-Mye' = "#cf5e4e", "both" = "#984EA3")
  ),
  annotation_width = unit(4, "mm")
)

col_split <- factor(
  paste0(plot_samples$time_group, "\n", plot_samples$label),
  levels = c("TIME-I\nfavorable", "TIME-I\nunfavorable", "TIME-Mye\nfavorable", "TIME-Mye\nunfavorable")
)

# Row split by TIME subtype, then direction
row_split <- data.frame(
  Subtype    = factor(row_source, levels = c("TIME-I", "TIME-Mye", "both")),
  Direction = factor(ifelse(row_dir == "up", "Up", "Down"), levels = c("Up", "Down"))
)

col_fun <- colorRamp2(
  c(-cap, -cap/2, 0, cap/2, cap),
  c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")
)

ht <- Heatmap(
  expr_scaled,
  name = "Z-score",
  col  = col_fun,
  column_split = col_split,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  column_title_gp = gpar(fontsize = 10),
  column_gap = unit(2, "mm"),
  row_split = row_split,
  cluster_row_slices = FALSE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  row_title_gp = gpar(fontsize = 10),
  row_gap = unit(2, "mm"),
  top_annotation = col_anno,
  left_annotation = left_anno,
  right_annotation = right_anno,
  border = TRUE,
  use_raster = TRUE,
  raster_quality = 3,
  heatmap_legend_param = list(
    title = "Scaled\nExpression",
    legend_height = unit(4, "cm"),
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(ht,
     column_title = "Transition Signature Genes (n=59)",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     padding = unit(c(2, 2, 2, 15), "mm"),
     merge_legend = TRUE)

# Compute transition score per sample via directional ssGSEA on all pseudobulk
up_g   <- transition_sig$up_genes
down_g <- transition_sig$down_genes
deg.list <- list(up = up_g, down = down_g)

score_fn <- function(expr_mat, gs_list) {
  ssgseaPar <- ssgseaParam(as.matrix(expr_mat), gs_list, normalize = TRUE)
  sc <- gsva(ssgseaPar) |> t()
  return(sc)
}

scores <- score_fn(pb_mat, deg.list)
scores <- data.frame(scores, check.names = FALSE)
scores$sample_hyphen <- rownames(scores)
scores$score <- scores$up - scores$down

all_samples <- sample_info |>
  mutate(sample_hyphen = str_replace_all(sample, "_", "-")) |>
  filter(sample_hyphen %in% scores$sample_hyphen)
all_samples$score <- scores$score[match(all_samples$sample_hyphen, scores$sample_hyphen)]
all_samples <- all_samples |> filter(!is.na(score))

# Correlation bar — transition score vs cell fractions (Baseline)
# Cell fractions from metadata (TME subset, all cell types including non-immune)
mtx_cells <- metadata |>
  filter(sample %in% all_samples$sample,
         non_malignant_count >= 100,
         subset == 'TME',
         !celltype_r2 %in% c('Melanocytes(CNA-)')) |>
  mutate(celltype_r2 = case_when(
    celltype_r2 %in% c('iCAF_MMP1', 'iCAF_IL6') ~ 'iCAF',
    .default = celltype_r2
  )) |>
  dplyr::select(sample, freq_r2_comp, celltype_r2) |>
  distinct(sample, celltype_r2, .keep_all = TRUE) |>
  pivot_wider(values_from = freq_r2_comp, names_from = celltype_r2, values_fill = 0) |>
  column_to_rownames(var = 'sample')

# Baseline samples
bl_samples <- all_samples |> filter(tx_status == "Baseline")
common_samples <- intersect(rownames(mtx_cells), bl_samples$sample)
cat("  Baseline samples with cell fractions:", length(common_samples), "\n")

mtx <- cbind(mtx_cells[common_samples, , drop = FALSE],
             Score = bl_samples$score[match(common_samples, bl_samples$sample)])

cor.res <- cor(mtx, use = "pairwise.complete.obs")
testRes <- corrplot::cor.mtest(mtx, conf.level = 0.95, method = "spearman")
cor.test.res <- testRes$p
diag(cor.test.res) <- NA

df.res <- data.frame(
  cor_eff = cor.res[, 'Score'],
  pvalue  = cor.test.res[, 'Score']
)
df.res <- df.res[rownames(df.res) != 'Score', ]
df.res$cell_type <- rownames(df.res)

cat("  Significant (p<0.05):\n")
print(df.res |> filter(pvalue < 0.05) |> arrange(desc(cor_eff)))

custom_colors <- c("Positive" = "#CC0C00FF", "Negative" = "#5C88DAFF")

# Keep only significant cell types
df_plot <- df.res |>
  filter(pvalue < 0.05) |>
  mutate(
    color_group = ifelse(cor_eff > 0, "Positive", "Negative"),
    cell_type = factor(cell_type, levels = cell_type[order(cor_eff)])
  )

p_cor <- ggplot(df_plot, aes(x = cell_type, y = cor_eff, fill = color_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = custom_colors, name = "") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "", x = "", y = "Correlation Coefficient") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text.y = element_text(size = 8, color = 'black', hjust = 1),
    axis.text.x = element_text(size = 8, color = 'black'),
    legend.position = "none"
  )

n_bars <- nrow(df_plot)

# Boxplot — transition score by response (Baseline only)
df_box <- all_samples |>
  filter(response %in% c("R", "NR"), tx_status == "Baseline") |>
  mutate(response = factor(response, levels = c("R", "NR")))

# Wilcoxon test
wt <- wilcox.test(score ~ response, data = df_box)
cat("  Wilcoxon p =", format(wt$p.value, digits = 3), "\n")

p_box <- ggplot(df_box, aes(x = response, y = score, fill = response)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.5) +
  scale_fill_manual(values = c("R" = "#CC0C00FF", "NR" = "#5C88DAFF")) +
  labs(x = "", y = "Transition Score\n(ssGSEA up - down)",
       title = "Baseline") +
  annotate("text", x = 1.5, y = max(df_box$score) * 1.05,
           label = paste0("p = ", format(wt$p.value, digits = 2)),
           size = 4.5) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = 'black'),
    legend.position = "none"
  )

# Treatment Response
# Directional ssGSEA: score = ssGSEA(up) - ssGSEA(down)
compute_directional_ssgsea <- function(expr_list, genes_up, genes_down) {
  gs_list <- list(favorable = genes_up, unfavorable = genes_down)
  score_list <- lapply(expr_list, function(expr_mat) {
    gu <- intersect(genes_up, rownames(expr_mat))
    gd <- intersect(genes_down, rownames(expr_mat))
    if (length(gu) < 3 || length(gd) < 3) return(NULL)
    gs_present <- list(favorable = gu, unfavorable = gd)
    ssgseaPar <- ssgseaParam(as.matrix(expr_mat), gs_present, normalize = FALSE)
    sc <- gsva(ssgseaPar) |> t() |> as.data.frame()
    sc$sample <- sub("^[^.]+\\.", "", rownames(sc))
    sc$score <- sc$favorable - sc$unfavorable
    data.frame(sample = sc$sample, score = sc$score, stringsAsFactors = FALSE)
  })
  scores <- do.call(rbind, score_list[!sapply(score_list, is.null)])
  rownames(scores) <- NULL
  scores
}

# Within-cohort z-scoring
scale_within_cohort <- function(df) {
  df |>
    group_by(cohort) |>
    mutate(score_z = as.numeric(scale(score))) |>
    ungroup() |>
    mutate(score_z = ifelse(is.na(score_z), 0, score_z))
}

# Logistic regression AUROC (cohort-adjusted, stratified 10-fold CV)
compute_response_auroc <- function(df, nfolds = 10) {
  df <- df |> filter(response %in% c("R", "NR"), !is.na(score_z))
  df$resp_bin <- ifelse(df$response == "R", 1, 0)
  if (nrow(df) < 30 || length(unique(df$resp_bin)) < 2) return(list(auc = NA, lower = NA, upper = NA))
  
  # Stratified CV by cohort composition
  set.seed(42)
  folds <- sample(rep(1:nfolds, length.out = nrow(df)))
  pred_prob <- rep(NA, nrow(df))
  
  for (k in 1:nfolds) {
    train_idx <- folds != k
    test_idx <- folds == k
    tryCatch({
      fit <- glm(resp_bin ~ score_z + cohort, data = df[train_idx, ], family = binomial)
      pred_prob[test_idx] <- predict(fit, newdata = df[test_idx, ], type = "response")
    }, error = function(e) {
      # Fallback: score-only model if cohort causes issues
      tryCatch({
        fit <- glm(resp_bin ~ score_z, data = df[train_idx, ], family = binomial)
        pred_prob[test_idx] <<- predict(fit, newdata = df[test_idx, ], type = "response")
      }, error = function(e2) NULL)
    })
  }
  
  valid <- !is.na(pred_prob)
  if (sum(valid) < 30) return(list(auc = NA, lower = NA, upper = NA))
  
  roc_obj <- tryCatch(
    pROC::roc(df$resp_bin[valid], pred_prob[valid], quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(roc_obj)) return(list(auc = NA, lower = NA, upper = NA))
  
  ci <- tryCatch(
    pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 1000, quiet = TRUE),
    error = function(e) c(NA, as.numeric(roc_obj$auc), NA)
  )
  
  list(auc = as.numeric(roc_obj$auc), lower = as.numeric(ci[1]), upper = as.numeric(ci[3]))
}

load(file.path(result_dir, "prepared_data.RData"))
load(file.path(result_dir, "feature_selection.RData"))
best_genes     <- transition_sig$genes
best_importance <- transition_sig$importance
gene_direction <- transition_sig$direction
genes_up       <- transition_sig$up_genes
genes_down     <- transition_sig$down_genes

transition_scores <- compute_directional_ssgsea(expr_list, genes_up, genes_down)
df_all <- merge_score_clin(transition_scores, clin_info, baseline_only = TRUE)
df_all <- scale_within_cohort(df_all)

# Response
cohorts_response <- df_all |>
  filter(response %in% c("R", "NR")) |>
  group_by(cohort) |>
  summarise(n = n(), n_r = sum(response == "R"), .groups = "drop") |>
  filter(n >= 15, n_r >= 3, (n - n_r) >= 3) |>
  pull(cohort)

response_forest <- lapply(cohorts_response, function(coh) {
  d <- df_all |> filter(cohort == coh, response %in% c("R", "NR"))
  d$resp_bin <- ifelse(d$response == "R", 1, 0)
  fit <- tryCatch(glm(resp_bin ~ score_z, data = d, family = binomial), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  ce <- coef(summary(fit))["score_z", ]
  data.frame(cohort = coh, n = nrow(d), n_r = sum(d$resp_bin),
             OR = exp(ce["Estimate"]),
             ci_lower = exp(ce["Estimate"] - 1.96 * ce["Std. Error"]),
             ci_upper = exp(ce["Estimate"] + 1.96 * ce["Std. Error"]),
             p = ce["Pr(>|z|)"], stringsAsFactors = FALSE)
})
response_forest_df <- do.call(rbind, response_forest[!sapply(response_forest, is.null)])

df_resp <- df_all |> filter(response %in% c("R", "NR"))
df_resp$resp_bin <- ifelse(df_resp$response == "R", 1, 0)
pooled_resp <- tryCatch({
  fit <- glm(resp_bin ~ score_z + cohort, data = df_resp, family = binomial)
  ce <- coef(summary(fit))["score_z", ]
  est <- ce[1]; se <- ce[2]; pval <- ce[4]
  data.frame(cohort = "Pooled", n = nrow(df_resp), n_r = sum(df_resp$resp_bin),
             OR = exp(est),
             ci_lower = exp(est - 1.96 * se),
             ci_upper = exp(est + 1.96 * se),
             p = pval)
}, error = function(e) NULL)

# PFS
cat("\n=== PFS forest plot ===\n")
df_pfs <- df_all |>
  mutate(pfs_event = ifelse(pfs %in% c(0, "Alive"), 0, 1)) |>
  filter(!is.na(time_pfs), !is.na(pfs_event))

cohorts_pfs <- df_pfs |>
  group_by(cohort) |>
  summarise(n = n(), n_events = sum(pfs_event), .groups = "drop") |>
  filter(n >= 15, n_events >= 5) |> pull(cohort)

pfs_forest <- lapply(cohorts_pfs, function(coh) {
  d <- df_pfs |> filter(cohort == coh)
  fit <- tryCatch(coxph(Surv(time_pfs, pfs_event) ~ score_z, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(cohort = coh, n = nrow(d), n_events = sum(d$pfs_event),
             HR = s$conf.int[1, 1], ci_lower = s$conf.int[1, 3],
             ci_upper = s$conf.int[1, 4], p = s$coefficients[1, 5], stringsAsFactors = FALSE)
})
pfs_forest_df <- do.call(rbind, pfs_forest[!sapply(pfs_forest, is.null)])
pooled_pfs <- tryCatch({
  fit <- coxph(Surv(time_pfs, pfs_event) ~ score_z + cohort, data = df_pfs)
  s <- summary(fit)
  data.frame(cohort = "Pooled", n = nrow(df_pfs), n_events = sum(df_pfs$pfs_event),
             HR = s$conf.int["score_z", 1], ci_lower = s$conf.int["score_z", 3],
             ci_upper = s$conf.int["score_z", 4], p = s$coefficients["score_z", 5])
}, error = function(e) NULL)
if (!is.null(pooled_pfs)) pfs_forest_df <- rbind(pfs_forest_df, pooled_pfs)

# OS
cat("\n=== OS forest plot ===\n")
df_os <- df_all |>
  mutate(os_event = ifelse(os %in% c(0, "Alive"), 0, 1)) |>
  filter(!is.na(time_os), !is.na(os_event))

cohorts_os <- df_os |>
  group_by(cohort) |>
  summarise(n = n(), n_events = sum(os_event), .groups = "drop") |>
  filter(n >= 15, n_events >= 5) |> pull(cohort)

os_forest <- lapply(cohorts_os, function(coh) {
  d <- df_os |> filter(cohort == coh)
  fit <- tryCatch(coxph(Surv(time_os, os_event) ~ score_z, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(cohort = coh, n = nrow(d), n_events = sum(d$os_event),
             HR = s$conf.int[1, 1], ci_lower = s$conf.int[1, 3],
             ci_upper = s$conf.int[1, 4], p = s$coefficients[1, 5], stringsAsFactors = FALSE)
})
os_forest_df <- do.call(rbind, os_forest[!sapply(os_forest, is.null)])
pooled_os <- tryCatch({
  fit <- coxph(Surv(time_os, os_event) ~ score_z + cohort, data = df_os)
  s <- summary(fit)
  data.frame(cohort = "Pooled", n = nrow(df_os), n_events = sum(df_os$os_event),
             HR = s$conf.int["score_z", 1], ci_lower = s$conf.int["score_z", 3],
             ci_upper = s$conf.int["score_z", 4], p = s$coefficients["score_z", 5])
}, error = function(e) NULL)
if (!is.null(pooled_os)) os_forest_df <- rbind(os_forest_df, pooled_os)

# Benchmark
all_sig_scores <- list()
all_sig_scores[["Transition score"]] <- transition_scores

for (sig_name in names(benchmark_sigs)) {
  cat("  Scoring:", sig_name, "\n")
  all_sig_scores[[sig_name]] <- compute_ssgsea(expr_list, benchmark_sigs[[sig_name]])
}

benchmark_sigs <- list(
  'Interferon Gamma' = c('IDO1', 'CXCL10', 'CXCL9', 'HLA-DRA', 'STAT1', 'IFNG'),
  'Expanded Immune Gene Signature' = c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK",
                                       "CD2", "HLA-DRA", "CXCL13", "IL2RG", "NKG7",
                                       "HLA-E", "CXCR6", "LAG3", "TAGAP", "CXCL10",
                                       "STAT1", "GZMB"),
  'T Cell Reactivity' = c("PDCD1", "CTLA4", "LAG3", "GZMB", "PRF1",
                          "ITGAE", "TNFRSF9", "TNFRSF18"),
  'Hot-tumor Signature' = c("CXCL9", "CXCL10", "CXCL11", "CXCR3", "CD3", "CD4",
                            "CD8A", "CD8B", "CD274", "PDCD1", "CXCR4", "CCL5"),
  'Tumor Inflamation Signature' = c("PSMB10", "HLA-DQA1", "HLA-DRB1", "CMKLR1", "HLA-E", "NKG7", "CD8A", 
                                    "CCL5", "CXCL9", "CD27", "CXCR6", "IDO1", "STAT1", "TIGIT", "LAG3", 
                                    "CD274", "PDCD1LG2", "CD276"),
  'INCITE-12' = c("KLRD1", "KLRK1", "FASLG", "GZMK", "KLRC1", "KLRC4", "CD8A", "PRF1", "APOBEC3G", "APOL3", "CD244", "FCRL6"),
  'IMPRES' = c("CD40", "CD27", "TNFRSF14", "CD276", "HAVCR2", "CD200", "VSIR",
               "PDCD1", "CD274", "CTLA4", "CD28", "CD80", "CD86", "TNFSF4")
)

compute_lr_auroc <- function(df, nfolds = 10) {
  df <- df |> filter(response %in% c("R", "NR"), !is.na(score_z))
  df$resp_bin <- ifelse(df$response == "R", 1, 0)
  if (nrow(df) < 30 || length(unique(df$resp_bin)) < 2)
    return(data.frame(auc = NA, lower = NA, upper = NA))
  
  set.seed(42)
  folds <- sample(rep(1:nfolds, length.out = nrow(df)))
  pred_prob <- rep(NA, nrow(df))
  
  for (k in 1:nfolds) {
    train <- folds != k; test <- folds == k
    tryCatch({
      fit <- glm(resp_bin ~ score_z + cohort, data = df[train, ], family = binomial)
      pred_prob[test] <- predict(fit, newdata = df[test, ], type = "response")
    }, error = function(e) {
      tryCatch({
        fit <- glm(resp_bin ~ score_z, data = df[train, ], family = binomial)
        pred_prob[test] <<- predict(fit, newdata = df[test, ], type = "response")
      }, error = function(e2) NULL)
    })
  }
  
  valid <- !is.na(pred_prob)
  if (sum(valid) < 30) return(data.frame(auc = NA, lower = NA, upper = NA))
  
  roc_obj <- tryCatch(pROC::roc(df$resp_bin[valid], pred_prob[valid], quiet = TRUE),
                      error = function(e) NULL)
  if (is.null(roc_obj)) return(data.frame(auc = NA, lower = NA, upper = NA))
  
  ci <- tryCatch(pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 1000, quiet = TRUE),
                 error = function(e) c(NA, as.numeric(roc_obj$auc), NA))
  
  data.frame(auc = round(as.numeric(roc_obj$auc), 3),
             lower = round(as.numeric(ci[1]), 3),
             upper = round(as.numeric(ci[3]), 3))
}

all_scores <- list()


all_scores[["Transition score"]] <- compute_directional_ssgsea(
  expr_list, transition_sig$up_genes, transition_sig$down_genes)

for (nm in names(benchmark_sigs)) {
  cat("  ", nm, "\n")
  all_scores[[nm]] <- compute_ssgsea(expr_list, benchmark_sigs[[nm]])
}

for (nm in names(xcell2_sigs)) {
  cat("  ", nm, "\n")
  all_scores[[nm]] <- xcell2_sigs[[nm]]
}

benchmark_list <- lapply(names(all_scores), function(sig_name) {
  cat("  ", sig_name, "... ")
  score_df <- all_scores[[sig_name]]
  df <- merge(clin_info, score_df, by = "sample", all.x = FALSE) |>
    filter(tx_status == "Baseline") |>
    scale_within_cohort()
  res <- compute_lr_auroc(df)
  res$Score <- sig_name
  cat("AUROC =", res$auc, "\n")
  res
})

benchmark_df <- do.call(rbind, benchmark_list) |>
  select(Score, auc, lower, upper) |>
  arrange(desc(auc))
rownames(benchmark_df) <- NULL

# Barplot — Response LR AUROC
df_bm <- benchmark_df |>
  mutate(Score = factor(Score, levels = Score[order(auc)]))

p_bm <- ggplot(df_bm, aes(x = Score, y = auc, fill = auc)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, linewidth = 0.4) +
  coord_flip(ylim = c(0.5, 0.75)) +
  scale_fill_gradientn(colors = c("#C6DBEF", "#3182BD", "#08306B"),
                       limits = c(min(df_bm$auc, na.rm = TRUE),
                                  max(df_bm$auc, na.rm = TRUE)),
                       guide = "none") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  labs(x = "", y = "AUROC \n(10-fold CV)",
       title = "Prediction Benchmark \n(response)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "right"
  )











