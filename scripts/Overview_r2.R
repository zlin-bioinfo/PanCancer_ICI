rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','tibble','pheatmap','qs2','MetBrewer','viridis','forcats','grid','corrplot','ComplexHeatmap','colorRamp2','corrr','igraph','ggraph','tidygraph','graphlayouts','ggforce','ggnetwork','ggnewscale')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

t <- c('CD4_T-naive','CD4_Tcm','CD4_Treg','CD4_T-ISG','CD4_Tfh','CD4_Tstr','CD4_Tctl','CD4_Th17',
           'CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
           "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
           "CD8_NK-like",'MAIT','gdT')
nk <- c('NK_CD56loCD16hi','NK_CD56hiCD16lo', 'Cycling T/NK')
bplasma <- c("B-naive", "B-ISG", "B-HSP", "B_MT2A", 
             "ACB_EGR1", "ACB_NR4A2", "ACB_CCR7", "B-memory", "B-AtM", 
             "GCB-pre", "GCB-DZ_SUGCT", "GCB-LZ_LMO2",
             "GCB-cycling", "PC-cycling",
             "PC-early_RGS13", "PC-early_LTB", "PC_IGHG", "PC_IGHA")
mye <- c('Mast','pDC','cDC1', 
         'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'mregDC', 'MoDC', 
         'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
         'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro-ISG', 
         'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2', 'Cycling myeloids')
nonimmune <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein",
               "Pericytes","SMC", "Myofibroblasts", "CAF_SFRP2", 
               "CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap")
# meta_int <- read.csv('tables/meta_all.csv') 
# meta_int |> 
#   distinct(celltype_r2, sample, .keep_all = T) |> 
#   select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |> 
#   filter(celltype_r2 %in% t[grepl('CD4', t)]) |> 
#   ggplot(aes(x = reorder(celltype_r2, freq_r2_comp, FUN = median, decreasing = T), y = freq_r2_comp)) + 
#   geom_jitter(size = 0.1) + geom_boxplot(aes(fill = celltype_r2), alpha = 0.4, show.legend = F) + 
#   theme_minimal() + theme(plot.title = element_text(hjust = 0.5),
#                           axis.text.x = element_text(size = 8, colour = 'black')) + 
#   ggtitle('CD4+T Cells') +
#   Seurat::RotatedAxis() + xlab('') + ylab('Relative Frequency')
# ggsave('figures/Abundance/cd4t.pdf', width = 8, height = 5)
# 
# meta_int |> 
#   distinct(celltype_r2, sample, .keep_all = T) |> 
#   select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |> 
#   filter(celltype_r2 %in% t[!grepl('CD4', t)]) |> 
#   ggplot(aes(x = reorder(celltype_r2, freq_r2_comp, FUN = median, decreasing = T), y = freq_r2_comp)) + 
#   geom_jitter(size = 0.1) + geom_boxplot(aes(fill = celltype_r2), alpha = 0.4, show.legend = F) + 
#   theme_minimal() + theme(plot.title = element_text(hjust = 0.5),
#                           axis.text.x = element_text(size = 8, colour = 'black')) + 
#   ggtitle('CD8+T&NK Cells') +
#   Seurat::RotatedAxis() + xlab('') + ylab('Relative Frequency')
# ggsave('figures/Abundance/cd8t.pdf', width = 7.5, height = 5)
# 
# meta_int |> 
#   distinct(celltype_r2, sample, .keep_all = T) |> 
#   select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |> 
#   filter(celltype_r2 %in% mye) |> 
#   ggplot(aes(x = reorder(celltype_r2, freq_r2_comp, FUN = median, decreasing = T), y = freq_r2_comp)) + 
#   geom_jitter(size = 0.1) + geom_boxplot(aes(fill = celltype_r2), alpha = 0.4, show.legend = F) + 
#   theme_minimal() + theme(plot.title = element_text(hjust = 0.5),
#                           axis.text.x = element_text(size = 8, colour = 'black')) + 
#   ggtitle('Myeloid Cells') +
#   Seurat::RotatedAxis() + xlab('') + ylab('Relative Frequency')
# ggsave('figures/Abundance/mye.pdf', width = 8, height = 5)
# 
# meta_int |> 
#   distinct(celltype_r2, sample, .keep_all = T) |> 
#   select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |> 
#   filter(celltype_r2 %in% bplasma) |> 
#   ggplot(aes(x = reorder(celltype_r2, freq_r2_comp, FUN = median, decreasing = T), y = freq_r2_comp)) + 
#   geom_jitter(size = 0.1) + geom_boxplot(aes(fill = celltype_r2), alpha = 0.4, show.legend = F) + 
#   theme_minimal() + theme(plot.title = element_text(hjust = 0.5),
#                           axis.text.x = element_text(size = 8, colour = 'black')) + 
#   ggtitle('B Lymphocytes') +
#   Seurat::RotatedAxis() + xlab('') + ylab('Relative Frequency')
# ggsave('figures/Abundance/bplasma.pdf', width = 8, height = 5)
# 
# meta_int |> 
#   distinct(celltype_r2, sample, .keep_all = T) |> 
#   select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |> 
#   filter(celltype_r2 %in% nonimmune) |> 
#   ggplot(aes(x = reorder(celltype_r2, freq_r2_comp, FUN = median, decreasing = T), y = freq_r2_comp)) + 
#   geom_jitter(size = 0.1) + geom_boxplot(aes(fill = celltype_r2), alpha = 0.4, show.legend = F) + 
#   theme_minimal() + theme(plot.title = element_text(hjust = 0.5),
#                           axis.text.x = element_text(size = 8, colour = 'black')) + 
#   ggtitle('Non-immune') +
#   Seurat::RotatedAxis() + xlab('') + ylab('Relative Frequency')
# ggsave('figures/Abundance/non-immune.pdf', width = 5, height = 3.5)

# Existance across dataset
list_metadata <- qs_read('tables/meta_list.qs2')
meta_int <- read.csv('tables/meta_all.csv') 
list_metadata$RCC_Bi <- NULL
meta_int <- meta_int |> filter(cohort != 'RCC_Bi')
dist_fun <- function(list_metadata, cellstates, .time_point){
  list_dist <- lapply(list_metadata, function(metadata){
    print(unique(metadata$cohort))
    # metadata <- metadata |> 
    #   group_by(sample, celltype_r2) |> 
    #   mutate(count_subset = n()) |>  
    #   ungroup() 
    metadata <- filter(metadata, time_point == .time_point)
    subtype_dist <- table(metadata$sample, metadata$celltype_r2)
    rowname <- rownames(subtype_dist)
    subtype_dist <- apply(subtype_dist, 2, function(x) as.numeric(x>=3)) 
    rownames(subtype_dist) <- rowname
    subtype_dist <- apply(subtype_dist, 2, sum) |> 
      sapply(function(x){return(ifelse(x >= 1, 1, ifelse(length(unique(metadata$sample)) <= 2 & x >= 1, 1, 0)))}) |> 
      data.frame() |> 
      rownames_to_column(var = 'subtype')
    names(subtype_dist)[2] <- unique(metadata$cohort)
    return(subtype_dist)
  })
  df_dist <- Reduce(function(x, y) merge(x, y, by = "subtype", all = TRUE), list_dist) |> column_to_rownames(var = 'subtype')
  df_dist[is.na(df_dist)] <- 0
  df_dist <- df_dist[cellstates,]
  df_dist <- t(df_dist) |> data.frame(check.names = F) |> rownames_to_column(var = 'cohort')
  return(df_dist)
}

dist_df <- lapply(c('Pre','On'), function(timepoint){
  dist_t <- dist_fun(list_metadata = list_metadata, 
                     cellstates = t, .time_point = timepoint)
  dist_nk <- dist_fun(list_metadata = list_metadata[c('SKCM_Becker', 'SKCM_Plozniak', 'BCC_Yost', 
                                                      'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 
                                                      'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
                                                      'NSCLC_Yan', 
                                                      'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo')], 
                      cellstates = nk, .time_point = timepoint)
  
  dist_bplasma <- dist_fun(list_metadata = list_metadata[c('SKCM_Becker', 'SKCM_Plozniak', 'BCC_Yost', 
                                                           'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 
                                                           'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
                                                           'NSCLC_Yan', 
                                                           'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo')], 
                           cellstates =  bplasma, .time_point = timepoint)
  dist_mye <- dist_fun(list_metadata = list_metadata[c('SKCM_Becker', 'SKCM_Plozniak', 'BCC_Yost', 
                                                       'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 
                                                       'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
                                                       'NSCLC_Yan', 
                                                       'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo')], 
                       cellstates =  mye, .time_point = timepoint)
  dist_nonimmune <- dist_fun(list_metadata = list_metadata[c('SKCM_Becker', 'SKCM_Plozniak', 'BCC_Yost', 
                                                             'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao', 
                                                             'HNSC_Franken', 
                                                             'NSCLC_Yan', 
                                                             'CRC_Li', 'CRC_Chen', 'PCa_Hawley')], 
                             cellstates =  nonimmune, .time_point = timepoint)
  
  df_dist <- Reduce(function(x, y) merge(x, y, by = "cohort", all = TRUE), list(dist_t, dist_nk, dist_bplasma, dist_mye, dist_nonimmune)) |> 
    column_to_rownames(var = 'cohort') |> 
    t() |> data.frame(check.names = F)
  df_dist[is.na(df_dist)] <- 0
  return(df_dist)
})

combined_dist <- (dist_df[[1]] + 2 * dist_df[[2]])

# Interpret the combined values
result <- apply(combined_dist, c(1, 2), function(x) {
  if (x == 0) "None"
  else if (x == 1) "Pre-Tx only"
  else if (x == 2) "On-Tx only"
  else if (x == 3) "Both"
})

result <- result[,c("SKCM_this study", "SKCM_Plozniak",   "BCC_Yost", "SCC_Yost", 
                    "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Zhang", "TNBC_Shiao", 
                    "HNSC_Franken", "HNSC_vanderLeun", "HNSC_Luoma" , 
                    "CRC_Li", "CRC_Chen", "NSCLC_Yan", "NSCLC_Liu", "PCa_Hawley",'HCC_Guo')]
df_anno <- meta_int |> 
  distinct(celltype_r2, component) |> 
  mutate(component = case_when(component == 'T_NK' ~ 'T&NK',
                               component == 'Bplasma' ~ 'B&plasma',
                               .default = component)) |> 
  mutate(component = factor(component, levels = c('T&NK', 'B&plasma', 'Myeloids', 'Non-immune'))) 

pdf('figures/dist_existing.pdf', height = 9, width = 4.5)
Heatmap(result, 
        show_column_names = T, show_row_names = T,
        cluster_rows = F, cluster_columns = F, 
        col = c("#1F78B4","lightgray", "#c5484f",'#308240'),
        heatmap_legend_param = list(
          at = c("Both", "None", "Pre-Tx only", "On-Tx only"),
          labels = c("Both", "None", "Pre-Tx only", "On-Tx only"),
          title = "Existing"
        ),
        column_names_gp = gpar(fontsize = 6),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 6),
        row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
        rect_gp = gpar(col = "white", lwd = 1),
        width = ncol(result)*unit(3, "mm"), 
        height = nrow(result)*unit(2, "mm"),
        row_split = df_anno$component[match(rownames(result), df_anno$celltype_r2)])
dev.off()

# R(o/e)
library(Startrac)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
metadata <- read.csv('tables/meta_all.csv')
# unwanted_celltypes <- c('Melanocytes(CNA-)', 'Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Epithelial(CNA-)', 'Malignant(CNA+)', 
#                         'Cycling T/NK', "GCB-cycling", "PC-cycling", 'Cycling_myeloids')
unwanted_celltypes <- c('Melanocytes(CNA-)', 'Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Epithelial(CNA-)', 'Malignant(CNA+)')
col_fun <- colorRamp2(c(0, 1, 2), colors = c("white", "#ffe7d4", "#FE5B14"))
metadata <- metadata |> 
  mutate(subtype = case_when(subtype %in% c('BCC','SCC') ~ 'BCC&SCC',
                             TRUE ~ subtype))
subtypes <- unique(metadata$subtype)
subtypes <- subtypes[-9]
ht_list_t <- lapply(subtypes, function(.subtype){
  print(.subtype)
  roie <- metadata |> 
    mutate(time_point = factor(time_point, levels = c('Pre','On'))) |> 
    filter(celltype_r2 %in% c(t,nk),
           subtype == .subtype) |> 
    calTissueDist(byPatient = F,
                  colname.cluster = "celltype_r2",
                  colname.patient = "patient",
                  colname.tissue = "time_point",
                  method = "chisq", 
                  min.rowSum = 10) |> 
    as.matrix()
  missing_rows <- setdiff(c(t,nk), rownames(roie))
  if (length(missing_rows) > 0) {
    extra_rows <- matrix(NA, nrow = length(missing_rows), ncol = ncol(roie))
    rownames(extra_rows) <- missing_rows
    roie <- rbind(roie, extra_rows)
  }
  
  # Reorder rows to match desired_row_order
  roie <- roie[c(t,nk), , drop = FALSE]
  print(dim(roie))
  ht <- Heatmap(roie, 
          column_title = .subtype,
          # column_title_rot = 45,
          column_title_gp = gpar(fontsize = 10),
          column_title_rot = 0,
          column_names_rot = 30,
          show_heatmap_legend = T, 
          cluster_rows = F, 
          cluster_columns = F,
          row_names_side = 'right', 
          show_column_names = T,
          show_row_names = T,
          col = col_fun,
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          na_col = "grey",
          heatmap_legend_param = list(
            title = "Ro/e",
            at =  seq(0, 2, by = 0.5),
            labels =  seq(0, 2, by = 0.5),
            legend_gp = gpar(fill = col_fun(2))
          ),
          cell_fun = function(j, i, x, y, width, height, fill) {
            value <- roie[i, j]
            text_label <- ifelse(is.na(value), "",  # No text for NA values
                                 ifelse(value > 2, "+++",
                                        ifelse(value >= 1.5, "++",
                                               ifelse(value >= 0.8, "+",
                                                      ifelse(value > 0, "+/-", "-")))))
            grid.text(text_label, x, y, gp = gpar(fontsize = 8))
          }
  )
  return(ht)
})
ht_list <- NULL
for (i in 1:(length(ht_list_t)-1)) {
  ht_list = ht_list + ht_list_t[[i]]
}

# RCC_Bi
.subtype <- 'RCC'
roie <- metadata |>
  mutate(time_point = factor(time_point, levels = c('NoICI','ICI_exposed'))) |> 
  filter(celltype_r2 %in% c(t,nk),
         subtype == .subtype) |> 
  calTissueDist(byPatient = F,
                colname.cluster = "celltype_r2",
                colname.patient = "patient",
                colname.tissue = "time_point",
                method = "chisq", 
                min.rowSum = 10) |> 
  as.matrix()
missing_rows <- setdiff(c(t,nk), rownames(roie))
if (length(missing_rows) > 0) {
  extra_rows <- matrix(NA, nrow = length(missing_rows), ncol = ncol(roie))
  rownames(extra_rows) <- missing_rows
  roie <- rbind(roie, extra_rows)
}

# Reorder rows to match desired_row_order
roie <- roie[c(t,nk), , drop = FALSE]
print(dim(roie))
ht <- Heatmap(roie, 
              column_title = .subtype,
              # column_title_rot = 45,
              column_title_gp = gpar(fontsize = 10),
              column_title_rot = 0,
              column_names_rot = 30,
              show_heatmap_legend = T, 
              cluster_rows = F, 
              cluster_columns = F,
              row_names_side = 'right', 
              show_column_names = T,
              show_row_names = T,
              col = col_fun,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              na_col = "grey",
              heatmap_legend_param = list(
                title = "Ro/e",
                at =  seq(0, 2, by = 0.5),
                labels =  seq(0, 2, by = 0.5),
                legend_gp = gpar(fill = col_fun(2))
              ),
              cell_fun = function(j, i, x, y, width, height, fill) {
                value <- roie[i, j]
                text_label <- ifelse(is.na(value), "",  # No text for NA values
                                     ifelse(value > 2, "+++",
                                            ifelse(value >= 1.5, "++",
                                                   ifelse(value >= 0.8, "+",
                                                          ifelse(value > 0, "+/-", "-")))))
                grid.text(text_label, x, y, gp = gpar(fontsize = 8))
              }
)
ht_list <- ht_list+ht
pdf('figures/Roie/t.pdf', height = 5, width = 9)
draw(ht_list, ht_gap = unit(0.4, "cm"))
dev.off()

# myeloids
subtypes <- unique(metadata$subtype)
subtypes <- subtypes[-9]
ht_list_t <- lapply(subtypes, function(.subtype){
  print(.subtype)
  roie <- metadata |> 
    mutate(time_point = factor(time_point, levels = c('Pre','On'))) |> 
    filter(celltype_r2 %in% mye,
           subtype == .subtype) |> 
    calTissueDist(byPatient = F,
                  colname.cluster = "celltype_r2",
                  colname.patient = "patient",
                  colname.tissue = "time_point",
                  method = "chisq", 
                  min.rowSum = 10) |> 
    as.matrix()
  missing_rows <- setdiff(mye, rownames(roie))
  if (length(missing_rows) > 0) {
    extra_rows <- matrix(NA, nrow = length(missing_rows), ncol = ncol(roie))
    rownames(extra_rows) <- missing_rows
    roie <- rbind(roie, extra_rows)
  }
  
  # Reorder rows to match desired_row_order
  roie <- roie[mye, , drop = FALSE]
  print(dim(roie))
  ht <- Heatmap(roie, 
                column_title = .subtype,
                # column_title_rot = 45,
                column_title_gp = gpar(fontsize = 10),
                column_title_rot = 0,
                column_names_rot = 30,
                show_heatmap_legend = T, 
                cluster_rows = F, 
                cluster_columns = F,
                row_names_side = 'right', 
                show_column_names = T,
                show_row_names = T,
                col = col_fun,
                row_names_gp = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                na_col = "grey",
                heatmap_legend_param = list(
                  title = "Ro/e",
                  at =  seq(0, 2, by = 0.5),
                  labels =  seq(0, 2, by = 0.5),
                  legend_gp = gpar(fill = col_fun(2))
                ),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  value <- roie[i, j]
                  text_label <- ifelse(is.na(value), "",  # No text for NA values
                                       ifelse(value > 2, "+++",
                                              ifelse(value >= 1.5, "++",
                                                     ifelse(value >= 0.8, "+",
                                                            ifelse(value > 0, "+/-", "-")))))
                  grid.text(text_label, x, y, gp = gpar(fontsize = 8))
                }
  )
  return(ht)
})
ht_list <- NULL
for (i in 1:(length(ht_list_t)-1)) {
  ht_list = ht_list + ht_list_t[[i]]
}

# RCC_Bi
.subtype <- 'RCC'
roie <- metadata |>
  mutate(time_point = factor(time_point, levels = c('NoICI','ICI_exposed'))) |> 
  filter(celltype_r2 %in% mye,
         subtype == .subtype) |> 
  calTissueDist(byPatient = F,
                colname.cluster = "celltype_r2",
                colname.patient = "patient",
                colname.tissue = "time_point",
                method = "chisq", 
                min.rowSum = 10) |> 
  as.matrix()
missing_rows <- setdiff(mye, rownames(roie))
if (length(missing_rows) > 0) {
  extra_rows <- matrix(NA, nrow = length(missing_rows), ncol = ncol(roie))
  rownames(extra_rows) <- missing_rows
  roie <- rbind(roie, extra_rows)
}

# Reorder rows to match desired_row_order
roie <- roie[mye, , drop = FALSE]
print(dim(roie))
ht <- Heatmap(roie, 
              column_title = .subtype,
              # column_title_rot = 45,
              column_title_gp = gpar(fontsize = 10),
              column_title_rot = 0,
              column_names_rot = 30,
              show_heatmap_legend = T, 
              cluster_rows = F, 
              cluster_columns = F,
              row_names_side = 'right', 
              show_column_names = T,
              show_row_names = T,
              col = col_fun,
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 8),
              na_col = "grey",
              heatmap_legend_param = list(
                title = "Ro/e",
                at =  seq(0, 2, by = 0.5),
                labels =  seq(0, 2, by = 0.5),
                legend_gp = gpar(fill = col_fun(2))
              ),
              cell_fun = function(j, i, x, y, width, height, fill) {
                value <- roie[i, j]
                text_label <- ifelse(is.na(value), "",  # No text for NA values
                                     ifelse(value > 2, "+++",
                                            ifelse(value >= 1.5, "++",
                                                   ifelse(value >= 0.8, "+",
                                                          ifelse(value > 0, "+/-", "-")))))
                grid.text(text_label, x, y, gp = gpar(fontsize = 8))
              }
)
ht_list <- ht_list+ht
pdf('figures/Roie/myeloids.pdf', height = 5, width = 9)
draw(ht_list, ht_gap = unit(0.4, "cm"))
dev.off()











roie <- metadata |> 
  mutate(time_point = factor(time_point, levels = c('Pre','On'))) |> 
  filter(!celltype_r2 %in% unwanted_celltypes,
         subtype == 'SKCM') |> 
  calTissueDist(byPatient = F,
                colname.cluster = "celltype_r2",
                colname.patient = "patient",
                colname.tissue = "time_point",
                method = "chisq", 
                min.rowSum = 0) |> 
  as.matrix() |> t() |> 
  apply(1, function(x) pmin(x, 2))
Heatmap(roie, 
        column_title = 'SKCM',
        column_title_rot = 45,
        show_heatmap_legend = T, 
        cluster_rows = T, 
        cluster_columns = F,
        row_names_side = 'right', 
        show_column_names = T,
        show_row_names = T,
        col = col_fun,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "R(o/e)",
          at =  seq(0, 2, by = 0.5),
          labels =  seq(0, 2, by = 0.5),
          legend_gp = gpar(fill = col_fun(2))
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", roie[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)

# Heatmap
library(Startrac)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
meta_int <- read.csv('tables/meta_all.csv') 
unwanted_celltypes <- c('Melanocytes(CNA-)', 'Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Epithelial(CNA-)', 'Malignant(CNA+)')
# logit_transform <- function(x, epsilon = 1e-6) {
#   x <- ifelse(x == 0, epsilon, x)  # Replace 0 with a small value
#   x <- ifelse(x == 1, 1 - epsilon, x)  # Replace 1 with (1 - small value)
#   log(x / (1 - x))  # Logit transformation
# }

meta_df <- meta_int |> 
  distinct(celltype_r2, sample, .keep_all = TRUE) |> 
  filter(!celltype_r2 %in% unwanted_celltypes) 
df_sample <- meta_df |> distinct(sample, .keep_all = T)
meta_matrix <- meta_df |> 
  select(sample, celltype_r2, freq_r2_comp) |> 
  pivot_wider(names_from = celltype_r2, values_from = freq_r2_comp, values_fill = list(freq_r2_comp = 0)) |> 
  column_to_rownames(var = 'sample') |> 
  as.matrix() |> t()
row_names <- rownames(meta_matrix)
# Apply logit transformation while keeping NA values
meta_matrix <- apply(meta_matrix, 2, function(col) scale(log(col + 1e-6)))

# Generate the heatmap
col_ha <- HeatmapAnnotation(dataset = df_sample$cohort, time_point = df_sample$time_point, response = df_sample$response)
celltype_class <- c(
  setNames(rep("T/NK", length(c(t,nk))), c(t,nk)),
  setNames(rep("B/plasma", length(bplasma)), bplasma),
  setNames(rep("Myeloid", length(mye)), mye),
  setNames(rep("Non-immune", length(nonimmune)), nonimmune)
)
# Match row names to their categories
celltype_annotation <- celltype_class[row_names] # your_data_matrix is the data for heatmap

# Define colors for annotations
category_colors <- c("T/NK" = "#E41A1C", 
                     "B/plasma" = "#377EB8", 
                     "Myeloid" = "#4DAF4A", 
                     "Non-immune" = "#984EA3")

# Create row annotation
row_ha <- rowAnnotation(
  Classification = celltype_annotation, 
  col = list(Classification = category_colors)
)
Heatmap(meta_matrix, show_column_names = F, top_annotation = col_ha, left_annotation = row_ha)
  










# Heatmap for change
meta_int <- read.csv('tables/meta_int.csv') 
pt_df <- read.csv('tables/meta_patient.csv')
# freq_filt <- meta_int |> 
#   select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
#   distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
#   pivot_wider(values_from = count_r2, names_from = time_point, values_fill = 0)
# freq_filt$pt_r2 <- paste0(freq_filt$patient, '_', freq_filt$celltype_r2)
freq_wide <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(change = log2((On + 0.01)/(Pre + 0.01)), diff = (On - Pre))
freq_wide |> 
  select(patient, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, dataset, Pre, On, change) |> 
  write.csv(file = 'data/df_patient_celltype.csv', row.names = F)
mat_change <- freq_wide |> 
  select(patient, diff, celltype_r2) |> 
  pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
  column_to_rownames(var = 'celltype_r2')
celltype_order_r2 <- c('CD4_Naive','CD4_Tm_CREM-','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM', 'CD4_Tm_CCL5', 
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
                       'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2',
                       'Endo_lymphatic','Endo_artery','Endo_capillary','Endo_tip','Endo_vein',
                       'EndMT','CAF_inflammatory', 'CAF_adipogenic', 'CAF_PN', 'CAF_AP', 'Myofibroblast')
mat_change <- mat_change[celltype_order_r2,]
mat_change[which(rownames(mat_change) %in% 'NK_CD56loCD16hi'):which(rownames(mat_change) %in% 'Myofibroblast'),
         grepl('SCC|NSCLC', colnames(mat_change))] <- NA
mat_change[which(rownames(mat_change) %in% 'Endo_lymphatic'):which(rownames(mat_change) %in% 'Myofibroblast'),
         grepl('TNBC_Shiao|TNBC_Zhang|HNSC_IMCISION|HNSC_Luoma', colnames(mat_change))] <- NA

pt_df <- pt_df[match(colnames(mat_change), pt_df$patient),]
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
             Modality = c('Mono' = "#9cc184", 'Dual' = "#3c7c3d"),
             Response = structure(names = as.character(unique(pt_df$response)), pal_startrek()(length(unique(pt_df$response)))),
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

df_anno <- meta_int |> 
  distinct(celltype_r2, component) |> 
  mutate(component = case_when(component == 'T_NK' ~ 'T&NK',
                               component == 'Bplasma' ~ 'B&plasma',
                               .default = component)) |> 
  mutate(component = factor(component, levels = c('T&NK', 'B&plasma', 'Myeloids', 'Non-immune')))  
  
pdf('figures/Change/ht_dynamics.pdf', height = 10, width = 9)
Heatmap(t(scale(t(mat_change))), na_col = 'lightgray', name = 'ΔRF (z-score)',
        row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
        show_column_names = F, show_row_names = T,
        cluster_rows = F, cluster_columns = F,
        col = colorRamp2(c(-1, 0, 1), c("#1F78B4", "white", "#E31A1C")),
        column_names_gp = gpar(fontsize = 8),
        heatmap_legend_param = list(
          legend_direction = "horizontal", 
          legend_width = unit(2, "cm"), at = c(-1, 0, 1),
          legend_side = 'bottom',
          title_position = "topcenter"),
        row_names_gp = gpar(fontsize = 6),
        top_annotation = col_ha,
        width = ncol(mat_change)*unit(0.8, "mm"), 
        height = nrow(mat_change)*unit(2.5, "mm"),
        row_split = df_anno$component[match(rownames(mat_change), df_anno$celltype_r2)])
dev.off()
# # Mean
# mat_change_mean <- freq_wide |> 
#   select(dataset, diff, celltype_r2) |> 
#   group_by(dataset,celltype_r2) |> 
#   dplyr::summarize(mean_change = mean(diff)) |> 
#   pivot_wider(values_from = mean_change, names_from = dataset, values_fill = 0) |> 
#   column_to_rownames(var = 'celltype_r2') 
# mat_change_mean <- mat_change_mean[celltype_order_r2, 
#                                    c("SKCM_Becker", "BCC_Yost", "SCC_Yost", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "HNSC_Franken", "HNSC_IMCISION", "HNSC_Luoma", "CRC_Li", "PCa_Hawley", "NSCLC_Liu")]
# mat_change_mean[which(rownames(mat_change_mean) %in% 'NK_CD56loCD16hi'):which(rownames(mat_change_mean) %in% 'Myofibroblast'),
#            grepl('SCC|NSCLC', colnames(mat_change_mean))] <- NA
# mat_change_mean[which(rownames(mat_change_mean) %in% 'EC_lymphatic'):which(rownames(mat_change_mean) %in% 'Myofibroblast'),
#            grepl('TNBC_Shiao|TNBC_Zhang|HNSC_IMCISION|HNSC_Luoma', colnames(mat_change_mean))] <- NA
# pdf('figures/Change/ht_dynamics_mean.pdf', height = 8, width = 4)
# Heatmap(t(scale(t(mat_change_mean))), na_col = 'lightgray', name = 'Mean ΔRF \n(z-score)',
#         row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
#         show_column_names = T, show_row_names = T,
#         cluster_rows = F, cluster_columns = F,
#         col = colorRamp2(c(-1, 0, 1), c("#1F78B4", "white", "#E31A1C")),
#         heatmap_legend_param = list(
#           legend_direction = "horizontal", 
#           legend_width = unit(2, "cm"), at = c(-1, 0, 1),
#           legend_side = 'bottom',
#           title_position = "topcenter"),
#         column_names_gp = gpar(fontsize = 6),
#         row_names_gp = gpar(fontsize = 6),
#         # top_annotation = col_ha,
#         width = ncol(mat_change)*unit(0.3, "mm"), 
#         height = nrow(mat_change)*unit(2.5, "mm"),
#         row_split = df_anno$component[match(rownames(mat_change_mean), df_anno$celltype_r2)])
# dev.off()

# Median
mat_change_median <- freq_wide |> 
  select(dataset, diff, celltype_r2) |> 
  group_by(dataset,celltype_r2) |> 
  dplyr::summarize(median_change = median(diff)) |> 
  pivot_wider(values_from = median_change, names_from = dataset, values_fill = 0) |> 
  column_to_rownames(var = 'celltype_r2') 
mat_change_median <- mat_change_median[celltype_order_r2, 
                                   c("SKCM_Becker", "BCC_Yost", "SCC_Yost", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "HNSC_Franken", "HNSC_IMCISION", "HNSC_Luoma", "CRC_Li", "CRC_Chen", "PCa_Hawley", "NSCLC_Liu")]
mat_change_median[which(rownames(mat_change_median) %in% 'NK_CD56loCD16hi'):which(rownames(mat_change_median) %in% 'Myofibroblast'),
                grepl('SCC|NSCLC', colnames(mat_change_median))] <- NA
mat_change_median[which(rownames(mat_change_median) %in% 'Endo_lymphatic'):which(rownames(mat_change_median) %in% 'Myofibroblast'),
                grepl('TNBC_Shiao|TNBC_Zhang|HNSC_IMCISION|HNSC_Luoma', colnames(mat_change_median))] <- NA
pdf('figures/Change/ht_dynamics_median.pdf', height = 9, width = 4.5)
Heatmap(t(scale(t(mat_change_median))), na_col = 'lightgray', name = 'Median ΔRF \n(z-score)',
        row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
        show_column_names = T, show_row_names = T,
        cluster_rows = F, cluster_columns = F,
        col = colorRamp2(c(-1, 0, 1), c("#1F78B4", "white", "#E31A1C")),
        heatmap_legend_param = list(
          legend_direction = "horizontal", 
          legend_width = unit(2, "cm"), at = c(-1, 0, 1),
          legend_side = 'bottom',
          title_position = "topcenter"),
        column_names_gp = gpar(fontsize = 7),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 7),
        # top_annotation = col_ha,
        width = ncol(mat_change_median)*unit(3.2, "mm"), 
        height = nrow(mat_change_median)*unit(2.5, "mm"),
        row_split = df_anno$component[match(rownames(mat_change_median), df_anno$celltype_r2)])
dev.off()

# Pre Post median
pdf('figures/Change/ht_median.pdf', height = 11, width = 6)
mat_pre_median <- freq_wide |> 
  select(dataset, Pre, celltype_r2) |> 
  group_by(dataset,celltype_r2) |> 
  dplyr::summarize(median_pre = median(Pre)) |> 
  pivot_wider(values_from = median_pre, names_from = dataset, values_fill = 0) |> 
  column_to_rownames(var = 'celltype_r2') 
mat_pre_median <- mat_pre_median[celltype_order_r2, 
                                 c("SKCM_Becker", "BCC_Yost", "SCC_Yost", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "HNSC_Franken", "HNSC_IMCISION", "HNSC_Luoma", "CRC_Li", "CRC_Chen", "PCa_Hawley", "NSCLC_Liu")]
mat_pre_median[which(rownames(mat_pre_median) %in% 'NK_CD56loCD16hi'):which(rownames(mat_pre_median) %in% 'Myofibroblast'),
               grepl('SCC|NSCLC', colnames(mat_pre_median))] <- NA
mat_pre_median[which(rownames(mat_pre_median) %in% 'Endo_lymphatic'):which(rownames(mat_pre_median) %in% 'Myofibroblast'),
               grepl('TNBC_Shiao|TNBC_Zhang|HNSC_IMCISION|HNSC_Luoma', colnames(mat_pre_median))] <- NA
ht_pre <- Heatmap(t(scale(t(mat_pre_median))), na_col = 'lightgray', column_title = 'Pre',
                  column_title_gp = gpar(fontface = 'bold'),
                  row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
                  show_column_names = T, show_row_names = F,
                  cluster_rows = F, cluster_columns = F,
                  col = colorRamp2(c(-1, 0, 1), c("#1F78B4", "white", "#E31A1C")),
                  show_heatmap_legend = F,
                  column_names_gp = gpar(fontsize = 8),
                  row_names_gp = gpar(fontsize = 8),
                  # top_annotation = col_ha,
                  width = ncol(mat_pre_median)*unit(3, "mm"), 
                  height = nrow(mat_pre_median)*unit(3, "mm"),
                  row_split = df_anno$component[match(rownames(mat_pre_median), df_anno$celltype_r2)])
mat_post_median <- freq_wide |> 
  select(dataset, On, celltype_r2) |> 
  group_by(dataset,celltype_r2) |> 
  dplyr::summarize(median_post = median(On)) |> 
  pivot_wider(values_from = median_post, names_from = dataset, values_fill = 0) |> 
  column_to_rownames(var = 'celltype_r2') 
mat_post_median <- mat_post_median[celltype_order_r2, 
                                 c("SKCM_Becker", "BCC_Yost", "SCC_Yost", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "HNSC_Franken", "HNSC_IMCISION", "HNSC_Luoma", "CRC_Li", "CRC_Chen", "PCa_Hawley", "NSCLC_Liu")]
mat_post_median[which(rownames(mat_post_median) %in% 'NK_CD56loCD16hi'):which(rownames(mat_post_median) %in% 'Myofibroblast'),
               grepl('SCC|NSCLC', colnames(mat_post_median))] <- NA
mat_post_median[which(rownames(mat_post_median) %in% 'Endo_lymphatic'):which(rownames(mat_post_median) %in% 'Myofibroblast'),
               grepl('TNBC_Shiao|TNBC_Zhang|HNSC_IMCISION|HNSC_Luoma', colnames(mat_post_median))] <- NA
ht_post <- Heatmap(t(scale(t(mat_post_median))), na_col = 'lightgray', name = 'Median RF \n(z-score)', column_title = 'Post',
                   column_title_gp = gpar(fontface = 'bold'),
                   row_title_gp = gpar(fontsize = 8, fontface = 'bold'),
                   show_column_names = T, show_row_names = T,
                   cluster_rows = F, cluster_columns = F,
                   col = colorRamp2(c(-1, 0, 1), c("#1F78B4", "white", "#E31A1C")),
                   show_heatmap_legend = T,
                   column_names_gp = gpar(fontsize = 8),
                   row_names_gp = gpar(fontsize = 8),
                   heatmap_legend_param = list( 
                     legend_width = unit(2, "cm"), at = c(-1, 0, 1),
                     legend_direction = "horizontal", 
                     legend_side = 'bottom',
                     title_position = "topcenter"),
                   width = ncol(mat_post_median)*unit(3, "mm"), 
                   height = nrow(mat_post_median)*unit(3, "mm"),
                   row_split = df_anno$component[match(rownames(mat_post_median), df_anno$celltype_r2)])
ht_pre + ht_post
dev.off()


# Correlation of dynamic changes (co-regulatiion)
meta_int <- read.csv('tables/meta_int.csv')
pt_df <- read.csv('tables/meta_patient.csv')
head(meta_int)
range(meta_int$freq_r2_comp)
# Immune cells
immune <- c('CD4_Naive','CD4_Tm_CREM-','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM', 'CD4_Tm_CCL5', 
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

freq_wide <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  # filter(abs(Pre-Post) >= 3, (Pre >= 3 | Post >= 3)) |> 
  mutate(change = log2((On + 0.01)/(Pre + 0.01)), diff = (On - Pre))
# freq_wide$pt_r2 <- paste0(freq_wide$patient, '_', freq_wide$celltype_r2)
# freq_wide <- filter(freq_wide, pt_r2 %in% freq_filt$pt_r2)

mat_change <- freq_wide |> 
  filter(!dataset %in% c('NSCLC_Liu'), 
         !patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su011', 'BCC/SCC_Yost_su012', 'BCC/SCC_Yost_su014'),
         celltype_r2 %in% immune) |> 
  select(patient, diff, celltype_r2) |> 
  pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
M <- cor(t(mat_change))
testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
pdf('figures/Co-regulating/freq_immune.pdf', height = 10, width = 9)
cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
# corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
#          tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
#          mar = c(0,0,0.7,0),
#          sig.level = c(0.01, 0.05), insig = 'label_sig', title = 'Immune Cells', diag = F, method = 'square', type = 'upper') |> 
#   corrRect(c(1,8)) |>
#   corrRect(c(20,31)) |>
#   corrRect(c(41,43)) |> 
#   corrRect(c(53,62)) 
corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
         tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
         mar = c(0,0,0.7,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = 'Immune Cells', diag = F, method = 'square', type = 'upper')
dev.off()

# Create a tidy data frame of correlations
tidy_cors <- mat_change |> 
  t() |> 
  correlate() |> 
  stretch()
df_p <- testRes$p
diag(df_p) <- NA
df <- df_p |> as.data.frame() |> 
  rownames_to_column(var = 'x') |> 
  pivot_longer(values_to = 'p', names_to = 'y', cols = -x) |> 
  mutate(FDR = p.adjust(p, method = "fdr")) |> 
  merge(tidy_cors, by = c('x','y'))
# Convert correlations stronger than some value to an undirected graph object
graph_cors <- df |> 
  filter(FDR<0.05, r>0.3) |> 
  graph_from_data_frame(directed = FALSE)
V(graph_cors)$`Cell type` <- ifelse(V(graph_cors)$name %in% t, 'T lymphocytes',
                                 ifelse(V(graph_cors)$name %in% nk, 'NK cells',
                                        ifelse(V(graph_cors)$name %in% mye, 'Myeloid cells',
                                               ifelse(V(graph_cors)$name %in% bplasma, 'B lymphocytes', 'Non-immune'))))
V(graph_cors)$`Cell type` <- factor(V(graph_cors)$`Cell type`, levels = c("T lymphocytes", "NK cells", "B lymphocytes", "Myeloid cells"))
communities <- cluster_leiden(graph_cors, resolution = 0.01)
communities$nb_clusters
V(graph_cors)$Cluster <- as.factor(communities$membership)
df_cluster <- data.frame(cluster = V(graph_cors)$Cluster, celltype = V(graph_cors)$name, cellgroup = V(graph_cors)$`Cell type`)
cluster_to_mark <- df_cluster |> 
  distinct(cluster, cellgroup, .keep_all = T) |> 
  group_by(cluster) |> 
  summarize(count = n()) |> 
  filter(count > 1) |> 
  pull(cluster) 
graph_cors |> 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = r, size = -log10(p), alpha = abs(r))) +
  scale_size(name = 'Edge Width\n -log(p)', breaks = c(5,10,15)) +
  scale_colour_gradientn(
    name = "Pearson's Correlation", limits = c(0, 1), breaks = c(-1,0,1), 
    colors = colorRampPalette(c("white", "#CC0000"))(50),
    guide = guide_colorbar(direction = "horizontal", 
                           title.position = "top", 
                           title.hjust = 0.5)) +
  new_scale_color() + 
  geom_nodes(aes(color = `Cell type`), size = 5) +
  scale_color_manual(values=met.brewer("Egypt", 4), name = 'Node Color') +
  geom_mark_hull(
    aes(x, y, group = Cluster, filter = Cluster %in% cluster_to_mark),
    concavity = 5,
    expand = unit(2.5, "mm"), 
    color = 'darkgray') +
  geom_nodetext_repel(aes(label = name), size = 3) +
  guides(alpha = "none", size = guide_legend(nrow = 1), color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme_blank() + 
  theme(legend.position="bottom") +
  ggtitle('(FDR<0.05 & rho>0.3)')
ggsave('figures/Co-regulating/graph_immune_pos.pdf', height = 8, width = 9)

# Convert correlations stronger than some value to an undirected graph object
graph_cors <- df |> 
  filter(FDR<0.05, r< -0.3) |> 
  graph_from_data_frame(directed = FALSE)
V(graph_cors)$`Cell type` <- ifelse(V(graph_cors)$name %in% t, 'T lymphocytes',
                                    ifelse(V(graph_cors)$name %in% nk, 'NK cells',
                                           ifelse(V(graph_cors)$name %in% mye, 'Myeloid cells',
                                                  ifelse(V(graph_cors)$name %in% bplasma, 'B lymphocytes', 'Non-immune'))))
# fg <- cluster_fast_greedy(graph_cors)
# communities <- cluster_louvain(graph_cors, weights = E(graph_cors)$weight, resolution = 1.2)
communities <- cluster_leiden(graph_cors, resolution = 0.1)
communities$nb_clusters
V(graph_cors)$Cluster <- as.factor(communities$membership)
# df <- data.frame(cluster = V(graph_cors)$Cluster, celltype = V(graph_cors)$name, cellgroup = V(graph_cors)$`Cell type`)
# cluster_to_mark <- df |> 
#   distinct(cluster, cellgroup, .keep_all = T) |> 
#   group_by(cluster) |> 
#   summarize(count = n()) |> 
#   filter(count > 1) |> 
#   pull(cluster) 
graph_cors |> 
  as_tbl_graph() |> 
  mutate(degree = centrality_degree()) |> 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = r, size = -log10(p), alpha = abs(r))) +
  scale_size(name = 'Edge Width\n -log(p)') +
  scale_colour_gradientn(
    name = "Pearson's Correlation", limits = c(-1, 0), breaks = c(-1,0,1), 
    colors = colorRampPalette(c("#364f77", "white"))(50),
    guide = guide_colorbar(direction = "horizontal", 
                           title.position = "top", 
                           title.hjust = 0.5)) +
  new_scale_color() + 
  geom_nodes(aes(color = `Cell type`), size = 5) +
  scale_color_manual(values=met.brewer("Egypt", 4), name = 'Node Color') +
  # geom_mark_hull(
  #   aes(x, y, group = Cluster, filter = Cluster %in% cluster_to_mark),
  #   concavity = 5,
  #   expand = unit(2.5, "mm"), 
  #   color = 'darkgray') +
  geom_nodetext_repel(aes(label = name), size = 3) +
  guides(alpha = "none") +
  theme_blank() + 
  ggtitle('Immune Cells \n(FDR<0.05 & rho<-0.3)')
ggsave('figures/Co-regulating/graph_immune_neg.pdf', height = 6, width = 8)

# # adding coordiantes
# g <- igraph::simplify(graph_cors)
# bb <- layout_as_backbone(g)
# E(graph_cors)$col <- F
# E(graph_cors)$col[bb$backbone] <- T
# # Plot
# p_neg <- ggraph(graph_cors,
#        layout = "manual",
#        x = bb$xy[, 1],
#        y = bb$xy[, 2]) +
#   geom_edge_link(aes(edge_alpha = abs(r), edge_width = -log10(p), color = r, col = col)) +
#   guides(edge_alpha = "none") +
#   scale_edge_colour_gradientn(name = "Pearson's rho", limits = c(-1, 0), 
#                               colors = colorRampPalette(c("#336699", "white"))(50), breaks = c(-1,0)) +
#   geom_node_point(aes(color = Cluster, size = centrality_pagerank())) +
#   scale_size(name = 'Centrality\n(PageRank)') +
#   geom_node_text(aes(label = name), repel = TRUE, size = 3) +
#   geom_mark_hull(
#     aes(x, y, group = Cluster),
#     concavity = 5,
#     expand = unit(1, "mm"),
#     alpha = 0.1
#   ) +
#   scale_color_manual(name = 'Cluster\n(Leiden)', values = brewer.pal(communities$nb_clusters, 'Paired')) +
#   guides(color = guide_legend(override.aes = list(size = 4))) +
#   theme_void() +
#   theme(legend.position = "bottom",
#         legend.text = element_text(size=10)) +
#   ggtitle('Immune Cells\np<0.001, rho<-0.3'); p_neg
# ggsave('figures/Co-regulating/graph_immune_neg.pdf', height = 8, width = 14)
# p_pos + p_neg + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

df <- freq_wide |> 
  filter(celltype_r2 %in% immune, !dataset %in% c('NSCLC_Liu'), 
         !patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su011', 'BCC/SCC_Yost_su012', 'BCC/SCC_Yost_su014'),) |> 
  select(celltype_r2, diff, patient) |> 
  pivot_wider(values_from = 'diff', names_from = celltype_r2, values_fill = 0) |> 
  column_to_rownames(var = 'patient') 

pdf('figures/Co-regulating/ht_isg.pdf', height = 4, width = 14)
a <- apply(df, 1, scales::rescale, to = c(0,1))
a <- a[c('B_ISG15', 'CD4_Th_ISG15', 'CD4_Treg_ISG15', 'CD8_Tex_ISG15', 'cDC2_ISG15', 'Macro_ISG15','cDC2_CXCL9'),]
score_pos <- apply(a, 2, sum)
a <- a[,order(score_pos)]
pt_df <- pt_df[match(colnames(a), pt_df$patient),]
col_ha = HeatmapAnnotation(
    Dataset = pt_df$dataset,
    # `Cancer Type` = pt_df$cancertype,
    # Treatment = pt_df$treatment,
    # Modality = pt_df$modality,
    # `Response Metrics` = pt_df$res_metric,
    Response = pt_df$response,
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
             # Treatment = structure(names = unique(pt_df$treatment), pal_npg()(length(unique(pt_df$treatment)))),
             # Modality = c('Mono' = "#9cc184", 'Dual' = "#3c7c3d"),
             # `Response Metrics` = structure(names = sort(unique(pt_df$res_metric)), met.brewer("Juarez", length(unique(pt_df$res_metric)))),
             Response = structure(names = c('RE','NR'), pal_startrek()(length(unique(pt_df$response))))
             ))
Heatmap(a, cluster_columns = F, cluster_rows = T, 
        show_column_names = F, show_row_names = T, name = 'Z-score', 
        top_annotation = col_ha, 
        width = ncol(a)*unit(1, "mm"), 
        height = nrow(a)*unit(5, "mm"), 
        col = rev(colorRampPalette(brewer.pal(6,'RdBu'))(100)),
        heatmap_legend_param = list(
          legend_direction = "horizontal", 
          legend_width = unit(2, "cm"), at = c(-1, 0, 1),
          title_position = "topcenter"))
dev.off()

score_pos <- apply(df[,c('B_ISG15', 'CD4_Th_ISG15', 'CD4_Treg_ISG15', 'CD8_Tex_ISG15', 'cDC2_ISG15', 'Macro_ISG15', 'cDC2_CXCL9')], 1, sum)
pt_df$isg <- score_pos[match(pt_df$patient, names(score_pos))]
pt_df |> 
  ggplot(aes(x = response, y = isg)) +
  geom_boxplot() + geom_point(aes(color = dataset)) + 
  facet_wrap(.~res_metric, scales = 'free_y') + 
  ggpubr::stat_compare_means(method = 't.test') + ylab('Overall Dynamics(ISG15+)') + theme_bw()
pt_df |>
  filter(res_metric %in% c('T-cell expansion','T-cell expansion+MRI')) |> 
  mutate(expansion = case_when(response == 'RE' ~ 'E',
                               response == 'NR' ~ 'NE')) |> 
  ggplot(aes(x = factor(expansion, level = c('NE','E')), y = isg)) + 
  geom_violin(aes(color = expansion)) +
  scale_color_manual(name = 'Expansion', values = c('#CF0034','#154999')) +
  ggnewscale::new_scale_color() +
  geom_jitter(aes(color = dataset), width = 0.1) +
  scale_color_jama(name = 'Cohort') +
  geom_boxplot(width = 0.5, alpha = 0.5) +
  # facet_wrap(.~dataset) +
  ggpubr::stat_compare_means(method = 't.test') + 
  ylab('Overall Dynamics \n(ISG15+)') + xlab('') + ggtitle('T Cell Expansion') +
  theme_bw()
ggsave('figures/Co-regulating/expansion_isg.pdf', height = 3.5, width = 5)

# pt_df |>
#   filter(res_metric %in% c('T-cell expansion','T-cell expansion+MRI')) |> 
#   mutate(expansion = case_when(response == 'RE' ~ 'E',
#                                response == 'NR' ~ 'NE')) |> 
#   ggplot(aes(x = factor(expansion, level = c('NE','E')), y = isg)) + 
#   geom_violin(aes(color = expansion)) +
#   scale_color_manual(name = 'Expansion', values = c('#CF0034','#154999')) +
#   ggnewscale::new_scale_color() +
#   # geom_jitter(aes(color = dataset), width = 0.1) +
#   # scale_color_jama(name = 'Cohort') +
#   geom_boxplot(width = 0.5, alpha = 0.5) +
#   facet_wrap(.~dataset) +
#   ggpubr::stat_compare_means(method = 't.test') + 
#   ylab('Overall Dynamics \n(ISG15+)') + xlab('') + ggtitle('T Cell Expansion') +
#   theme_bw()
# ggsave('figures/Co-regulating/expansion_isg_cohort.pdf', height = 3.5, width = 8)

pdf('figures/Co-regulating/ht_mac_naive_plasmsa.pdf', height = 5, width = 14)
a <- apply(df, 1, scales::rescale, to = c(0,1)) 
a <- a[c('B_Naive','CD4_Naive', 'CD8_Naive', 'CD4_Tm_AREG','CD4_Tm_CREM','B_Memory', 'B_AtM','NK_CD56loCD16hi', 'CD8_Temra_CX3CR1', 'CD4_Temra_CX3CR1','CD4_Treg_TNFRSF9', 'Macro_C1QC', 'Plasma_cell'),]
score_pos <- apply(a[c('B_Naive','CD4_Naive', 'CD8_Naive', 'CD4_Tm_AREG','CD4_Tm_CREM','B_Memory', 'B_AtM','NK_CD56loCD16hi', 'CD8_Temra_CX3CR1', 'CD4_Temra_CX3CR1'),], 2, sum)
score_neg <- apply(a[c('CD4_Treg_TNFRSF9', 'Macro_C1QC', 'Plasma_cell'),], 2, sum)
a <- a[,order(score_neg-score_pos)]
pt_df <- pt_df[match(colnames(a), pt_df$patient),]
col_ha = HeatmapAnnotation(
  Dataset = pt_df$dataset,
  # `Cancer Type` = pt_df$cancertype,
  # Treatment = pt_df$treatment,
  # Modality = pt_df$modality,
  # `Response Metrics` = pt_df$res_metric,
  # Response = pt_df$response,
  col = list(Dataset = structure(names = unique(pt_df$dataset), dittoColors()[1:length(unique(pt_df$dataset))])
             # `Cancer Type` = c("PCa" = "#a82203",
             #                   "CRC" = "#208cc0",
             #                   "NSCLC" = "#f1af3a",
             #                   "HNSC" = "#cf5e4e",
             #                   "SCC" = "#3B7546",
             #                   "BCC" = "#0092A5",
             #                   "SKCM" = "#000000",
             #                   "TNBC" = "#FF0196",
             #                   "ER+BC" = "#F073EA",
             #                   "HER2+BC" = "#FFB0FF"),
             # Treatment = structure(names = unique(pt_df$treatment), pal_npg()(length(unique(pt_df$treatment)))),
             # Modality = c('Mono' = "#9cc184", 'Dual' = "#3c7c3d"),
             # `Response Metrics` = structure(names = sort(unique(pt_df$res_metric)), met.brewer("Juarez", length(unique(pt_df$res_metric)))),
             # Response = structure(names = c('RE','NR'), pal_startrek()(length(unique(pt_df$response))))
             ))
Heatmap(a, cluster_columns = F, cluster_rows = F, row_names_side = 'left',
        show_column_names = F, show_row_names = T, name = 'Z-score', 
        top_annotation = col_ha, 
        width = ncol(a)*unit(1, "mm"), 
        height = nrow(a)*unit(4, "mm"), 
        col = rev(colorRampPalette(brewer.pal(6,'RdBu'))(100)),
        heatmap_legend_param = list(
          legend_direction = "horizontal", 
          legend_width = unit(2, "cm"), at = c(-1, 0, 1),
          title_position = "topcenter"))
dev.off()

a <- apply(df, 1, scales::rescale, to = c(-1,1))
a <- a[c('B_ISG15','CD4_Th_ISG15', 'CD4_Treg_ISG15','CD8_Tex_ISG15', 'cDC2_ISG15', 'Macro_ISG15','CD4_TfhTh1_IFNG'),]
score <- apply(a[1:7,], 2, sum)
a <- a[,order(score)]
pt_df <- pt_df[match(colnames(a), pt_df$patient),]
# col_ha = HeatmapAnnotation(
#   Dataset = pt_df$dataset,
#   `Cancer Type` = pt_df$cancertype,
#   Treatment = pt_df$treatment,
#   Modality = pt_df$modality,
#   `Response Metrics` = pt_df$res_metric,
#   Response = pt_df$response,
#   Cluster = pt_df$cluster,
#   col = list(Dataset = structure(names = unique(pt_df$dataset), dittoColors()[1:length(unique(pt_df$dataset))]),
#              `Cancer Type` = c("PCa" = "#a82203", 
#                                "CRC" = "#208cc0",     
#                                "NSCLC" = "#f1af3a",  
#                                "HNSC" = "#cf5e4e",
#                                "SCC" = "#3B7546",     
#                                "BCC" = "#0092A5", 
#                                "SKCM" = "#000000", 
#                                "TNBC" = "#FF0196", 
#                                "ER+BC" = "#F073EA",    
#                                "HER2+BC" = "#FFB0FF"),
#              Treatment = structure(names = unique(pt_df$treatment), pal_npg()(length(unique(pt_df$treatment)))),
#              Modality = c('Mono' = "#9cc184", 'Dual' = "#3c7c3d"),
#              Response = structure(names = c('RE','NR'), pal_startrek()(length(unique(pt_df$response)))),
#              `Response Metrics` = structure(names = sort(unique(pt_df$res_metric)), met.brewer("Juarez", length(unique(pt_df$res_metric)))))
#   )
col_ha = HeatmapAnnotation(
  Dataset = pt_df$dataset,
  `Cancer Type` = pt_df$cancertype,
  Treatment = pt_df$treatment,
  Modality = pt_df$modality,
  `Response Metrics` = pt_df$res_metric,
  Response = pt_df$response,
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
             Modality = c('Mono' = "#9cc184", 'Dual' = "#3c7c3d"),
             Response = structure(names = c('RE','NR'), pal_startrek()(length(unique(pt_df$response)))),
             `Response Metrics` = structure(names = sort(unique(pt_df$res_metric)), met.brewer("Juarez", length(unique(pt_df$res_metric))))))
Heatmap(a, cluster_columns = F, cluster_rows = F, 
        show_column_names = F, show_row_names = T, name = 'Z-score', 
        top_annotation = col_ha, 
        width = ncol(a)*unit(1, "mm"), 
        height = nrow(a)*unit(6, "mm"), 
        col = rev(colorRampPalette(brewer.pal(6,'RdBu'))(100)),
        heatmap_legend_param = list(
          legend_direction = "horizontal", 
          legend_width = unit(2, "cm"), at = c(-1, 0, 1),
          title_position = "topcenter"))


df <- freq_wide |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley', 'HNSC_Franken'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014")) |> 
  select(celltype_r2, diff, patient) |> 
  pivot_wider(values_from = 'diff', names_from = celltype_r2, values_fill = 0) |> 
  column_to_rownames(var = 'patient') 
a <- apply(df, 1, scales::rescale, to = c(-1,1))
a <- a[c('B_Naive','B_Memory', 'B_AtM', 'CD4_Naive', 'CD8_Naive', 'Mono_CD14', 'Mono_CD16', 'gdT', 'NK_CD56loCD16hi', 'CD8_Temra_CX3CR1', 'CD4_Tm_CREM-', 'CD8_Trm_ZNF683','CD4_Treg_TNFRSF9', 'Macro_C1QC', 'Plasma_cell', 'CAF_inflammatory'),]
score_pos <- apply(a[c('B_Naive','B_Memory', 'B_AtM', 'CD4_Naive', 'CD8_Naive', 'Mono_CD14', 'Mono_CD16', 'gdT', 'NK_CD56loCD16hi', 'CD8_Temra_CX3CR1'),], 2, sum)
score_neg <- apply(a[c('CD4_Tm_CREM-', 'CD8_Trm_ZNF683', 'CD4_Treg_TNFRSF9', 'Macro_C1QC',  'CAF_inflammatory','Plasma_cell'),], 2, sum)
# score_neg <- a[13,]
a <- a[,order(score_pos-score_neg)]
pt_df <- pt_df[match(colnames(a), pt_df$patient),]
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
             Modality = c('Mono' = "#9cc184", 'Dual' = "#3c7c3d"),
             Response = structure(names = c('RE','NR'), pal_startrek()(length(unique(pt_df$response)))),
             `Response Metrics` = structure(names = sort(unique(pt_df$res_metric)), met.brewer("Juarez", length(unique(pt_df$res_metric)))))
)
pdf('figures/Co-regulating/ht_mac_naive_plasmsa_caf.pdf', height = 7, width = 12)
Heatmap(a, cluster_columns = F, cluster_rows = F, 
        show_column_names = F, show_row_names = T, name = 'Z-score', 
        top_annotation = col_ha, 
        width = ncol(a)*unit(2, "mm"), 
        height = nrow(a)*unit(6, "mm"), col = rev(colorRampPalette(brewer.pal(6,'RdBu'))(100)))
dev.off()

a <- apply(df, 1, scales::rescale, to = c(-1,1))
a <- a[c('CD4_Tm_AREG','CD4_Tm_CREM', 'cDC1', 'cDC2_IL1B', 'CD8_Tcm_IL7R', 'CD4_Tm_CREM-', 'CD8_Trm_ZNF683', 'MigrDC', 'CD4_Th17_IL26'),]
a <- a[,order(a[1,]-a[nrow(a),])]
pt_df <- pt_df[match(colnames(a), pt_df$patient),]
col_ha = HeatmapAnnotation(
  Response = pt_df$response,
  col = list(
    Response = structure(names = c('RE','NR'), pal_startrek()(length(unique(pt_df$response)))))
)
Heatmap(a, cluster_columns = T, cluster_rows = T, 
        show_column_names = F, show_row_names = T, name = 'Z-score', 
        top_annotation = col_ha, 
        width = ncol(a)*unit(1, "mm"), 
        height = nrow(a)*unit(2.5, "mm"), col = rev(colorRampPalette(brewer.pal(6,'RdBu'))(100)))
# 
# # modality
# pdf('figures/Co-regulating/freq_immune_modality.pdf', height = 10, width = 18)
# par(mfrow=c(1,2))
# cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
# mat_change <- freq_wide |> 
#   filter(dataset != 'NSCLC_Liu', 
#          modality == 'Mono',
#          celltype_r2 %in% immune) |> 
#   select(patient, diff, celltype_r2) |> 
#   pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
#   column_to_rownames(var = 'celltype_r2')
# M <- cor(t(mat_change))
# testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
# corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
#          tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
#          mar = c(0,0,0.7,0),
#          sig.level = c(0.01, 0.05), insig = 'label_sig', title = 'Mono', diag = F, method = 'square', type = 'upper') 
# mat_change <- freq_wide |> 
#   filter(dataset != 'NSCLC_Liu', 
#          modality == 'Dual',
#          celltype_r2 %in% immune) |> 
#   select(patient, diff, celltype_r2) |> 
#   pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
#   column_to_rownames(var = 'celltype_r2')
# M <- cor(t(mat_change))
# testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
# corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
#          tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
#          mar = c(0,0,0.7,0),
#          sig.level = c(0.01, 0.05), insig = 'label_sig', title = 'Dual', diag = F, method = 'square', type = 'upper') 
# dev.off()
# 
# # response
# pdf('figures/Co-regulating/freq_immune_response.pdf', height = 10, width = 18)
# par(mfrow=c(1,2))
# cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
# mat_change <- freq_wide |> 
#   filter(dataset != 'NSCLC_Liu', 
#          response == 'RE',
#          celltype_r2 %in% immune) |> 
#   select(patient, diff, celltype_r2) |> 
#   pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
#   column_to_rownames(var = 'celltype_r2')
# M <- cor(t(mat_change))
# testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
# corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
#          tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
#          mar = c(0,0,0.7,0),
#          sig.level = c(0.01, 0.05), insig = 'label_sig', title = 'RE', diag = F, method = 'square', type = 'upper') 
# mat_change <- freq_wide |> 
#   filter(dataset != 'NSCLC_Liu', 
#          response == 'NR',
#          celltype_r2 %in% immune) |> 
#   select(patient, diff, celltype_r2) |> 
#   pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
#   column_to_rownames(var = 'celltype_r2')
# M <- cor(t(mat_change))
# testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
# corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
#          tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
#          mar = c(0,0,0.7,0),
#          sig.level = c(0.01, 0.05), insig = 'label_sig', title = 'NR', diag = F, method = 'square', type = 'upper') 
# dev.off()

# TME
mat_change <- freq_wide |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley', 'HNSC_Franken', 'CRC_Chen', 'CRC_Li'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014")) |> 
  select(patient, diff, celltype_r2) |> 
  pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
  column_to_rownames(var = 'celltype_r2') 
M <- cor(t(mat_change))
testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
pdf('figures/Co-regulating/freq_TME.pdf', height = 10, width = 10)
cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
corrplot(M, p.mat = testRes$p, method = 'square', col = cols(100), tl.srt = 45, 
         tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2", mar = c(0,0,1,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = 'TME', diag = F) 
dev.off()

# Create a tidy data frame of correlations
tidy_cors <- mat_change |> 
  t() |> 
  correlate() |> 
  stretch()
df_p <- testRes$p
diag(df_p) <- NA
df <- df_p |> as.data.frame() |> 
  rownames_to_column(var = 'x') |> 
  pivot_longer(values_to = 'p', names_to = 'y', cols = -x) |> 
  mutate(FDR = p.adjust(p, method = "fdr")) |> 
  merge(tidy_cors, by = c('x','y'))
graph_cors <- df |> 
  filter(FDR<0.05, r>0.3) |> 
  graph_from_data_frame(directed = FALSE)
V(graph_cors)$`Cell type` <- ifelse(V(graph_cors)$name %in% t, 'T lymphocytes',
                                    ifelse(V(graph_cors)$name %in% nk, 'NK cells',
                                           ifelse(V(graph_cors)$name %in% mye, 'Myeloid cells',
                                                  ifelse(V(graph_cors)$name %in% bplasma, 'B lymphocytes', 'Non-immune'))))
V(graph_cors)$`Cell type` <- factor(V(graph_cors)$`Cell type`, levels = c("T lymphocytes", "NK cells", "B lymphocytes", "Myeloid cells","Non-immune"))
communities <- cluster_leiden(graph_cors, resolution = 0.1)
communities$nb_clusters
V(graph_cors)$Cluster <- as.factor(communities$membership)
df_cluster <- data.frame(cluster = V(graph_cors)$Cluster, celltype = V(graph_cors)$name, cellgroup = V(graph_cors)$`Cell type`)
# cluster_to_mark <- df_cluster |> 
#   distinct(cluster, cellgroup, .keep_all = T) |> 
#   group_by(cluster) |> 
#   summarize(count = n()) |> 
#   filter(count > 1) |> 
#   pull(cluster) 
graph_cors |> 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = r, size = -log10(p), alpha = abs(r))) +
  scale_size(name = 'Edge Width\n -log(p)') +
  scale_colour_gradientn(
    name = "Pearson's Correlation", limits = c(0, 1), breaks = c(-1,0,1), 
    colors = colorRampPalette(c("white", "#CC0000"))(50),
    guide = guide_colorbar(direction = "horizontal", 
                           title.position = "top", 
                           title.hjust = 0.5)) +
  new_scale_color() + 
  geom_nodes(aes(color = `Cell type`), size = 5) +
  scale_color_manual(values=met.brewer("Egypt", 5), name = 'Node Color') +
  # geom_mark_hull(
  #   aes(x, y, group = Cluster, filter = Cluster %in% cluster_to_mark),
  #   concavity = 3,
  #   expand = unit(2.5, "mm"),
  #   color = 'darkgray') +
  geom_nodetext_repel(aes(label = name), size = 3) +
  guides(alpha = "none") +
  theme_blank() + 
  ggtitle('TME \n(FDR<0.05 & rho>0.3)')
ggsave('figures/Co-regulating/graph_TME_pos.pdf', height = 8, width = 9)

graph_cors <- df |> 
  filter(FDR<0.05, r< -0.3) |> 
  graph_from_data_frame(directed = FALSE)
V(graph_cors)$`Cell type` <- ifelse(V(graph_cors)$name %in% t, 'T lymphocytes',
                                    ifelse(V(graph_cors)$name %in% nk, 'NK cells',
                                           ifelse(V(graph_cors)$name %in% mye, 'Myeloid cells',
                                                  ifelse(V(graph_cors)$name %in% bplasma, 'B lymphocytes', 'Non-immune'))))
# fg <- cluster_fast_greedy(graph_cors)
# communities <- cluster_louvain(graph_cors, weights = E(graph_cors)$weight, resolution = 1.2)
communities <- cluster_leiden(graph_cors, resolution = 0.1)
communities$nb_clusters
V(graph_cors)$Cluster <- as.factor(communities$membership)
# df <- data.frame(cluster = V(graph_cors)$Cluster, celltype = V(graph_cors)$name, cellgroup = V(graph_cors)$`Cell type`)
# cluster_to_mark <- df |> 
#   distinct(cluster, cellgroup, .keep_all = T) |> 
#   group_by(cluster) |> 
#   summarize(count = n()) |> 
#   filter(count > 1) |> 
#   pull(cluster) 
graph_cors |> 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = r, size = -log10(p), alpha = abs(r))) +
  scale_size(name = 'Edge Width\n -log(p)') +
  scale_colour_gradientn(
    name = "Pearson's Correlation", limits = c(-1, 0), breaks = c(-1,0,1), 
    colors = colorRampPalette(c("#336699", "white"))(50),
    guide = guide_colorbar(direction = "horizontal", 
                           title.position = "top", 
                           title.hjust = 0.5)) +
  new_scale_color() + 
  geom_nodes(aes(color = `Cell type`), size = 5) +
  scale_color_manual(values=met.brewer("Egypt", 5), name = 'Node Color') +
  # geom_mark_hull(
  #   aes(x, y, group = Cluster, filter = Cluster %in% cluster_to_mark),
  #   concavity = 5,
  #   expand = unit(2.5, "mm"), 
  #   color = 'darkgray') +
  geom_nodetext_repel(aes(label = name), size = 3) +
  guides(alpha = "none") +
  theme_blank() + 
  ggtitle('TME \n(FDR<0.05, rho< -0.3)')
ggsave('figures/Co-regulating/graph_TME_neg.pdf', height = 8, width = 9)

freq_wide |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley', 'HNSC_Franken'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014")) |> 
  select(patient, diff, celltype_r2) |> 
  pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
  filter(celltype_r2 %in% c('CAF_inflammatory','CD4_Treg_TNFRSF9')) |> 
  column_to_rownames(var = 'celltype_r2') |> 
  t() |> data.frame() |> 
  rownames_to_column(var = 'celltypes') |> 
  dplyr::mutate(`Cancer type` = str_split(celltypes, '_', simplify = T)[,1]) |> 
  ggplot(aes(CAF_inflammatory, CD4_Treg_TNFRSF9)) + 
  geom_point(aes(color = `Cancer type`)) + 
  scale_color_manual(values = brewer.pal(4, 'Set1')) +
  geom_smooth( method = lm) + 
  ggpubr::stat_cor() + 
  theme_minimal()
ggsave('figures/Co-regulating/scatter_treg_caf.pdf', height = 4, width = 5)

freq_wide |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley', 'HNSC_Franken', 'CRC_Li', 'CRC_Chen'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014")) |> 
  select(patient, diff, celltype_r2) |> 
  pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
  filter(celltype_r2 %in% c('CAF_inflammatory','CD4_Treg_TNFRSF9')) |> 
  column_to_rownames(var = 'celltype_r2') |> 
  t() |> data.frame() |> 
  rownames_to_column(var = 'celltypes') |> 
  dplyr::mutate(`Cancer type` = str_split(celltypes, '_', simplify = T)[,1]) |> 
  ggplot(aes(CAF_inflammatory, CD4_Treg_TNFRSF9)) + 
  geom_point(aes(color = `Cancer type`)) + 
  scale_color_manual(values = brewer.pal(5, 'Set1')) +
  geom_smooth( method = lm) + 
  ggpubr::stat_cor() + 
  theme_minimal() 
ggsave('figures/Co-regulating/scatter_caf_treg.pdf', height = 4, width = 5)

# Pre vs diff
list_res <- lapply(unique(freq_wide$celltype_r2), function(subtype){
  cor_res <- c()
  df <- freq_wide |> filter(celltype_r2 == subtype) 
  summary_cor <- cor.test(df$Pre, df$Post)
  cor_res <- c(subtype, summary_cor$p.value, summary_cor$estimate, sqrt(var(df$Pre)))
  return(cor_res)
})
res <- do.call(rbind, list_res) |> data.frame()
colnames(res) <- c('subtype', 'pvalue', 'rho', 'std_Pre')
res$pvalue <- as.numeric(res$pvalue)
res$rho <- as.numeric(res$rho)
res$std_Pre <- as.numeric(res$std_Pre)
res$fdr <- p.adjust(res$pvalue, method = 'fdr', n = nrow(res))
res |> 
  ggplot(aes(x = rho, y = log10(fdr), label = subtype)) +
  geom_point() +
  geom_hline(yintercept=log10(0.05), linetype="dashed", color = "red") +
  theme_classic() + geom_text_repel() + 
  xlab("Correlation Coefficient \n(Pearson's r)") + ylab('-log10(FDR)') +
  ggtitle('Pre vs Change') + theme(plot.title = element_text(hjust = 0.5))

library(ggrepel)


