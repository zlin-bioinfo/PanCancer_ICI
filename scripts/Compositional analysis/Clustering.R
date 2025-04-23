rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','tibble','qs2','MetBrewer','viridis','grid','ComplexHeatmap','circlize','colorRamp2')
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
# Heatmap
meta_int <- read.csv('tables/meta_all.csv') 
unwanted_celltypes <- c('Melanocytes(CNA-)', 'Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Epithelial(CNA-)', 'Malignant(CNA+)')

meta_df <- meta_int |> 
  distinct(celltype_r2, sample, .keep_all = TRUE) |> 
  filter(!celltype_r2 %in% unwanted_celltypes,
         subset == 'All TME') 
df_sample <- meta_df |> distinct(sample, .keep_all = T)
meta_matrix <- meta_df |> 
  select(sample, celltype_r2, freq_r2_comp) |> 
  pivot_wider(names_from = celltype_r2, values_from = freq_r2_comp, values_fill = list(freq_r2_comp = 0)) |> 
  column_to_rownames(var = 'sample') |> 
  as.matrix() |> t()
row_names <- rownames(meta_matrix)
# Apply logit transformation while keeping NA values
meta_matrix <- apply(meta_matrix, 2, function(col) scale(log(col + 1e-6))) 
rownames(meta_matrix) <- row_names

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
ComplexHeatmap::Heatmap(meta_matrix, show_column_names = F, top_annotation = col_ha, left_annotation = row_ha)






