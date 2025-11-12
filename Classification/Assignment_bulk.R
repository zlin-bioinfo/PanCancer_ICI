pkgs <- c('qs2','tidyr','dplyr','plyr','stringr','tibble','janitor','Seurat', 'clusterProfiler','enrichplot',
          'ggplot2','RColorBrewer','ComplexHeatmap','MetBrewer','GSVA','glmnet','xCell2','survminer','survival',
          'progeny','viridis','GEOquery','AnnotationDbi','org.Hs.eg.db')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1) 
seu <- qs_read('data/pseudobulk_TME.qs2')
gs_list <- readxl::read_xlsx('tables/Signature.xlsx') |>
  lapply(function(x){return(x[!is.na(x)])})
ssgseaPar <- ssgseaParam(GetAssayData(seu, layer = 'data'), gs_list)
score_pb <- gsva(ssgseaPar)
gs_sub <- c("T cells","Treg cells", "Effector T cells","NK cells", "Th1","Th2","Th17", 
            "B cells", "Plasma cells", "MHC-I", "MHC-II", "TLS formation", 
            "DC", "Pan-macrophage", "M1-like", "M2-like",
            "Co-activation molecules", "Checkpoint molecules", "Interferon response",
            "Endothelium", "Angiogenesis", "CAF", "Matrix remodeling","Stromal immiune exclusion","Stromal suppression")
sample_info <- read.csv('tables/sample_info_time_updated.csv')
score_pb <- score_pb |> t() |> data.frame(check.names = F)
score_pb$group <- sample_info$group[match(rownames(score_pb), str_replace_all(sample_info$sample, '_', '-'))]
score_median <- aggregate(. ~ group, data = score_pb, FUN = median) |> column_to_rownames(var = 'group') |> apply(1, as.numeric)
rownames(score_median) <- names(gs_list)
pdf('figures/NMF/sig_score_median.pdf', height = 4, width = 10)
Heatmap(score_median[gs_sub,] |> t() |> scale(), 
        col = rev(brewer.pal(9, 'RdBu')), 
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 12),
        clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2',
        rect_gp = gpar(col = "white", lwd = 0.5),
        heatmap_legend_param = list(title = 'Enrichment Score \n(median)', 
                                    at = c(-1.5, 0, 1.5), labels = c("Min", "", "Max"), border = "black",
                                    title_position = "topcenter",
                                    legend_direction = "horizontal"),
        width = ncol(score_median[gs_sub,])*unit(35, "mm"),
        height = nrow(score_median[gs_sub,])*unit(1, "mm"))
dev.off()

expr_list <- qs_read('data/bulk_datasets/expr_list.qs2')
clin_list <- qs_read('data/bulk_datasets/clin_list.qs2')
score <- function(expr_mat){
  ssgseaPar <- ssgseaParam(as.matrix(expr_mat), gs_list)
  score <- gsva(ssgseaPar)
  return(score)
}
score_list <- lapply(expr_list, score)

all_sample_assignments <- lapply(score_list, function(score_matrix) {
  score_matrix_scaled <- score_matrix[gs_sub,] |> t() |> scale() |> t() |> data.frame(check.names = F)
  mtx_combined <- cbind(score_median[gs_sub,] |> t() |> scale() |> t(), score_matrix_scaled)
  M_iter <- cor(mtx_combined, method = 'spearman')
  mtx_raw_cor <- M_iter[colnames(score_median), colnames(score_matrix)] |>
    data.frame(check.names = F)
  mtx_normalized_score <- t(mtx_raw_cor) |>
    # apply(2, minmax_normalize) |>
    scale() |>
    data.frame(check.names = F)
  df_raw_cor_long <- mtx_raw_cor |>
    rownames_to_column(var = "Subtype") |>
    pivot_longer(
      cols = -Subtype,
      names_to = "sample",
      values_to = "Actual_Correlation"
    )
  df_norm_score_long <- mtx_normalized_score |>
    rownames_to_column(var = "sample") |>
    pivot_longer(
      cols = -sample,
      names_to = "Subtype",
      values_to = "Normalized_Score"
    )
  df_combined <- df_raw_cor_long |>
    inner_join(df_norm_score_long, by = c("sample", "Subtype"))
  sample_assignment_iter <- df_combined |>
    group_by(sample) |>
    slice_max(Actual_Correlation, n = 1, with_ties = FALSE) |>
    ungroup()
  return(sample_assignment_iter)
})
sample_assignment <- do.call(rbind, all_sample_assignments)

# Add clinical info
clin_info <- do.call(rbind, clin_list)
clin_info <- clin_info[clin_info$sample %in% sample_assignment$sample,]
clin_info$patient <- paste0(clin_info$cohort, '_',clin_info$patient)
clin_info$cancertype[clin_info$cancertype == 'mUC'] <- 'UC'
clin_info$cancertype[clin_info$cancertype == 'Esophageal'] <- 'EC'
clin_info$cohort <- paste0(clin_info$cancertype, '_', clin_info$cohort)
clin_info$response[is.na(clin_info$response) | clin_info$response %in% c('Unknown','Not Evaluable')] <- 'NE'
clin_info$response[clin_info$response %in% c('CR','MR','PR','Responder','CRPR','1','Complete Response','Partial Response')] <- 'R'
clin_info$response[clin_info$response %in% c('NR','PD','SD','Non-Responder','0','Progressive Disease','Stable Disease','N','Progressor','Stable')] <- 'NR'
clin_info$tx_status <- factor(clin_info$tx_status, levels = c('Baseline','Treated','Not specified'))
clin_info$group <- sample_assignment$Subtype[match(clin_info$sample, sample_assignment$sample)]
clin_info$cor <- sample_assignment$Actual_Correlation[match(clin_info$sample, sample_assignment$sample)]
# filter out samples
clin_info <- clin_info |> filter(cor > 0.25)
clin_info$group <- factor(clin_info$group, levels = c('Immune Quiescent', 'Immune Inflamed', 'B Cell-enriched', 'Myeloid-enriched'))
clin_info <- clin_info |> 
  filter(! response == '') |> 
  mutate(response = factor(response, levels = c('R','NR','NE'))) |> 
  arrange(group, tx_status, response, cohort, cancertype)
write.csv(clin_info, 'tables/sample_info_time_bulk.csv', row.names = F)
clin_info <- read.csv('tables/sample_info_time_bulk.csv')
clin_info <- clin_info |> 
  mutate(group = case_when(group == 'Immune Quiescent' ~ 'TIME-Q',
                           group == 'Immune Inflamed' ~ 'TIME-I',
                           group == 'B Cell-enriched' ~ 'TIME-B',
                           group == 'Myeloid-enriched' ~ 'TIME-Mye')) |> 
  mutate(group = factor(group, levels = c('TIME-Q','TIME-I','TIME-B','TIME-Mye')))

col_ha = HeatmapAnnotation(
  `Cancer Type` = clin_info$cancertype,
  Cohort = clin_info$cohort,
  Response = clin_info$response,
  `Treatment Status` = factor(clin_info$tx_status, levels = c('Baseline','Treated','Not specified')),
  `TIME Subtype` = clin_info$group,
  col = list(`Cancer Type` = structure(names = as.character(unique(clin_info$cancertype)), met.brewer("Juarez", length(unique(clin_info$cancertype)))),
             `Cohort` = structure(names = as.character(unique(clin_info$cohort)), met.brewer("Juarez", length(unique(clin_info$cohort)))),
             Response = c('R' = '#CC0C00FF','NR' = '#5C88DAFF','NE' = '#84BD00FF'),
             `Treatment Status` = c('Baseline' = "#04a3bd", 'Treated' = "#f0be3d", 'Not specified' = "#931e18"),
             `TIME Subtype` = structure(names = as.character(unique(clin_info$group)), met.brewer("Juarez",  length(unique(clin_info$group))))
  )
)
scaled_score_list <- lapply(score_list, function(mat) {
  t(scale(t(mat)))
})
scores <- do.call(cbind, scaled_score_list)
# row split
celltype <- c(
  "T cells","Treg cells", "Effector T cells","NK cells", "Th1","Th2","Th17", 
  "B cells", "Plasma cells", "MHC-I", "MHC-II", "TLS formation", 
  "DC", "Pan-macrophage", "M1-like", "M2-like",
  "Co-activation molecules", "Checkpoint molecules", "Interferon response",
  "Endothelium", "Angiogenesis", "CAF", "Matrix remodeling","Stromal immiune exclusion","Stromal suppression"
)
celltype.group <- factor(c(
  rep("T/NK",7),
  rep("B/plasma",5),
  rep('myeloid',4),
  rep("Signature",3),
  rep("Stromal",6)))
df.row<- data.frame(group = celltype.group); rownames(df.row) <- celltype
pdf('figures/NMF/RNA-seq/ht_sample_bulk.pdf', height = 6, width = 14)
ht <- Heatmap(scores[gs_sub, clin_info$sample], name = "Enrichment Score\n(scaled)", 
              column_title = paste0('Pan-cancer Bulk Datasets (n=', nrow(clin_info),')'),
              clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2',
              column_split = clin_info$group, 
              row_split = df.row,
              cluster_rows = T, cluster_columns = F,
              show_row_names = T, show_column_names = F, 
              show_row_dend = T, show_column_dend = F,
              col = circlize::colorRamp2(c(-2, 0, 2), c("#154999", "white", "#CF0034")),
              # col = rev(brewer.pal(10, 'RdBu')),
              na_col = 'lightgray',
              heatmap_legend_param = list(
                title_position = "topcenter",
                at = c(-2, 0, 2),
                labels = c("Min", "", "Max"),
                legend_direction = "horizontal", border = T),
              top_annotation = col_ha, 
              row_names_gp = gpar(fontsize = 10), use_raster=T); ht
dev.off()
# OS
sample_info_os <- clin_info |> 
  filter(!is.na(os), !is.na(time_os), tx_status == 'Baseline') |> 
  mutate(os = case_when(os == 'Alive' ~ 0,
                        os == 'Dead' ~ 1,
                        TRUE ~ as.numeric(os)),
         time_os = as.numeric(time_os)) |> 
  distinct(patient, .keep_all = T)
fit <- survfit(Surv(time_os, os) ~ group, data = sample_info_os)
p <- ggsurvplot(fit, data = sample_info_os,conf.int = F, pval = T, 
                legend.title = "",
                legend.labs = c('TIME-Q','TIME-I','TIME-B','TIME-Mye'),
                surv.median.line = "hv",
                risk.table = T,
                risk.table.col = "strata",
                risk.table.height = 0.3,
                ggtheme = theme_classic(),
                palette =  met.brewer("Juarez",4),
                xlab = 'Overall Survival(month)', title='', 
                break.time.by = 12, pval.coord = c(48, 0.5)); p
pdf('figures/Fig3/surv_os_ici.pdf', height = 5, width = 4.5)
print(p, newpage = FALSE)
dev.off()

# PFS
sample_info_pfs <- clin_info |> 
  filter(!is.na(pfs), !is.na(time_pfs), tx_status == 'Baseline') |> 
  mutate(pfs = case_when(pfs == 'Alive' ~ 0,
                         pfs == 'Dead' ~ 1,
                         TRUE ~ as.numeric(pfs)),
         time_pfs = as.numeric(time_pfs)) |> 
  distinct(patient, .keep_all = T)
fit <- survfit(Surv(time_pfs, pfs) ~ group, data = sample_info_pfs)
p <- ggsurvplot(fit, data = sample_info_pfs,conf.int = F, pval = T, 
                legend.title = "",
                legend.labs = c('TIME-Q','TIME-I','TIME-B','TIME-Mye'),
                surv.median.line = "hv",
                risk.table = T,
                risk.table.col = "strata",
                risk.table.height = 0.3,
                ggtheme = theme_classic(),
                palette =  met.brewer("Juarez",4),
                xlab = 'Prognosis-free Survival(month)', title='', 
                break.time.by = 12, pval.coord = c(24, 0.5)); p
pdf('figures/Fig3/surv_pfs_ici.pdf', height = 5, width = 4.5)
print(p, newpage = FALSE)
dev.off()