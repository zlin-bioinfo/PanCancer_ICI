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
pdf('figures/Roie/t.pdf', height = 6, width = 8)
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
pdf('figures/Roie/myeloids.pdf', height = 6, width = 8)
draw(ht_list, ht_gap = unit(0.4, "cm"))
dev.off()

# Bplasma
subtypes <- unique(metadata$subtype)
subtypes <- subtypes[-9]
ht_list_t <- lapply(subtypes, function(.subtype){
  print(.subtype)
  roie <- metadata |> 
    mutate(time_point = factor(time_point, levels = c('Pre','On'))) |> 
    filter(celltype_r2 %in% bplasma,
           subtype == .subtype) |> 
    calTissueDist(byPatient = F,
                  colname.cluster = "celltype_r2",
                  colname.patient = "patient",
                  colname.tissue = "time_point",
                  method = "chisq", 
                  min.rowSum = 10) |> 
    as.matrix()
  missing_rows <- setdiff(bplasma, rownames(roie))
  if (length(missing_rows) > 0) {
    extra_rows <- matrix(NA, nrow = length(missing_rows), ncol = ncol(roie))
    rownames(extra_rows) <- missing_rows
    roie <- rbind(roie, extra_rows)
  }
  
  # Reorder rows to match desired_row_order
  roie <- roie[bplasma, , drop = FALSE]
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
  filter(celltype_r2 %in% bplasma,
         subtype == .subtype) |> 
  calTissueDist(byPatient = F,
                colname.cluster = "celltype_r2",
                colname.patient = "patient",
                colname.tissue = "time_point",
                method = "chisq", 
                min.rowSum = 10) |> 
  as.matrix()
missing_rows <- setdiff(bplasma, rownames(roie))
if (length(missing_rows) > 0) {
  extra_rows <- matrix(NA, nrow = length(missing_rows), ncol = ncol(roie))
  rownames(extra_rows) <- missing_rows
  roie <- rbind(roie, extra_rows)
}

# Reorder rows to match desired_row_order
roie <- roie[bplasma, , drop = FALSE]
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
pdf('figures/Roie/bplasma.pdf', height = 4.5, width = 8)
draw(ht_list, ht_gap = unit(0.4, "cm"))
dev.off()

# nonimmune
subtypes <- unique(metadata$subtype)
subtypes <- subtypes[-c(9,10)]
ht_list_t <- lapply(subtypes, function(.subtype){
  print(.subtype)
  roie <- metadata |> 
    mutate(time_point = factor(time_point, levels = c('Pre','On'))) |> 
    filter(celltype_r2 %in% nonimmune,
           subtype == .subtype) |> 
    calTissueDist(byPatient = F,
                  colname.cluster = "celltype_r2",
                  colname.patient = "patient",
                  colname.tissue = "time_point",
                  method = "chisq", 
                  min.rowSum = 10) |> 
    as.matrix()
  missing_rows <- setdiff(nonimmune, rownames(roie))
  if (length(missing_rows) > 0) {
    extra_rows <- matrix(NA, nrow = length(missing_rows), ncol = ncol(roie))
    rownames(extra_rows) <- missing_rows
    roie <- rbind(roie, extra_rows)
  }
  
  # Reorder rows to match desired_row_order
  roie <- roie[nonimmune, , drop = FALSE]
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
  filter(celltype_r2 %in% nonimmune,
         subtype == .subtype) |> 
  calTissueDist(byPatient = F,
                colname.cluster = "celltype_r2",
                colname.patient = "patient",
                colname.tissue = "time_point",
                method = "chisq", 
                min.rowSum = 10) |> 
  as.matrix()
missing_rows <- setdiff(nonimmune, rownames(roie))
if (length(missing_rows) > 0) {
  extra_rows <- matrix(NA, nrow = length(missing_rows), ncol = ncol(roie))
  rownames(extra_rows) <- missing_rows
  roie <- rbind(roie, extra_rows)
}

# Reorder rows to match desired_row_order
roie <- roie[nonimmune, , drop = FALSE]
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
pdf('figures/Roie/nonimmune.pdf', height = 4, width = 7.5)
draw(ht_list, ht_gap = unit(0.4, "cm"))
dev.off()

