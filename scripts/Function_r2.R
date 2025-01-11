rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','qs', 'lmerTest','grid','msigdbr','ggplotify','pheatmap','msigdbr','MetBrewer','tibble','ComplexHeatmap','UCell')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

datasets <- c('SKCM_Becker', 'BCC_Yost', 'SCC_Yost', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_IMCISION', 'HNSC_Luoma', 'HNSC_Franken', 'NSCLC_Liu', 'CRC_Li','PCa_Hawley')
# CD4T
gs_list <- readxl::read_xlsx('/bigdata/zlin/PanCancer_ICI/tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 7, skip = 1) |> 
  lapply(function(x){return(x[!is.na(x)])})
gs_list$Exhaustion <- NULL
names(gs_list)[names(gs_list  ) == 'OXPHOS'] = 'Oxidative phosphorylation'
names(gs_list)[names(gs_list) == 'Lipid metabolism'] = 'Fatty acid metabolism'
gs_list$`Fatty acid metabolism`[gs_list$`Fatty acid metabolism` == 'PLA2G16'] <- 'PLAAT3'
gs_list$`Fatty acid metabolism`[gs_list$`Fatty acid metabolism` == 'ARNTL'] <- 'BMAL1'
gs_list$`Oxidative phosphorylation` <- gsub('\\.', '-', gs_list$`Oxidative phosphorylation`)

list_cd4 <- lapply(datasets, function(dataset) {
  print(dataset)
  seu <- qread(paste0('/bigdata/zlin/PanCancer_ICI/data/', dataset, '/seu_r2.qs')) |> 
    subset(subset = celltype_main == 'CD4+T') |> 
    subset(subset = patient %in% c("CRC_Li_P28", "HNSC_Franken_P16", "HNSC_Franken_P20", "SKCM_Becker_P10", "SKCM_Becker_P6"), invert = T) |> 
    NormalizeData() |> 
    AddModuleScore_UCell(features = gs_list)
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  return(seu@meta.data[,c('celltype_r2', 'dataset', names(gs_list))])
})
qsave(list_cd4, '/bigdata/zlin/PanCancer_ICI/tables/score_cd4.qs')
list_cd4 <- qread('/bigdata/zlin/PanCancer_ICI/tables/score_cd4.qs')
df_score <- do.call(rbind, list_cd4)

Differentiation <- c("Naïve", "Activation/Effector function")
Function <- c("TCR signaling", "Cytotoxicity", "Cytokine/Cytokine receptor",
              "Chemokine/Chemokine receptor", "Stress response", "Adhesion",
              "IFN response", "Treg signature", "Costimulatory molecules")
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
col_order <- c('CD4_Naive','CD4_Tm_CREM-','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM', 'CD4_Tm_CCL5', 
               'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
               'CD4_Treg_Early', 'CD4_Treg_TNFRSF9', 'CD4_Treg_ISG15', 
               'CD4_Th_ISG15', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif')
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(df_score$celltype_r2)),
                              nrow = length(MarkerNameVector))
colnames(FunctionScoreMatrix) <- unique(df_score$celltype_r2)
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)) {
  for(ri in 1:nrow(FunctionScoreMatrix)) {
    FunctionVec <- as_tibble(df_score) |> pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[df_score$celltype_r2 == unique(df_score$celltype_r2)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, scales::rescale, to=c(0, 1)))
FunctionScoreMatrix <- FunctionScoreMatrix[, col_order]
cols <- colorRampPalette(brewer.pal('Reds', n=100))(10)
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis))))
signatureType_row$Signature.type <- factor(signatureType_row$Signature.type, levels = c('Differentiation', 'Function','Metabolism','Apoptosis'))
pdf('/bigdata/zlin/PanCancer_ICI/figures/Functional_score_T/ht_func_cd4.pdf', height = 6, width = 8)
Heatmap(FunctionScoreMatrix, name = 'Mean \n(z-score)',
        column_title = 'CD4+ T cells',
        cluster_rows = F, cluster_columns = F,
        column_names_rot = 45,
        col = cols, 
        row_names_gp = gpar(fontsize = 10), 
        column_names_gp = gpar(fontsize = 9), 
        width = ncol(FunctionScoreMatrix)*unit(6, "mm"), 
        height = nrow(FunctionScoreMatrix)*unit(6, "mm"),
        row_title_gp = gpar(fontsize = 7, fontface = 'bold'),
        heatmap_legend_param = list(
          legend_direction = "horizontal",
          legend_width = unit(2, "cm"), at = c(0, 1),
          legend_side = 'bottom',
          title_position = "topcenter"),
        row_split = signatureType_row$Signature.type)
dev.off()

# CD8T
gs_list <- readxl::read_xlsx('/bigdata/zlin/PanCancer_ICI/tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 5, skip = 1) |> 
  lapply(function(x){return(x[!is.na(x)])})
names(gs_list)[names(gs_list) == 'Activation:Effector function'] <- 'Activation/Effector function'
names(gs_list)[names(gs_list) == 'TCR Signaling'] = 'TCR signaling'
names(gs_list)[names(gs_list) == 'IFN Response'] = 'IFN response'
gs_list$`IFN response`[gs_list$`IFN response` == 'DDX58'] <- 'RIGI'
gs_list$`Oxidative phosphorylation` <- gsub('\\.', '-', gs_list$`Oxidative phosphorylation`)
gs_list <- gs_list[-which(names(gs_list)=='Senescence')]

list_cd8 <- lapply(datasets, function(dataset) {
  print(dataset)
  seu <- qread(paste0('/bigdata/zlin/PanCancer_ICI/data/', dataset, '/seu_r2.qs')) |> 
    subset(subset = celltype_main == 'CD8+T') |>
    subset(subset =  (celltype_r2 %in% c('gdT', 'MAIT')), invert = T) |> 
    subset(subset = patient %in% c("CRC_Li_P28", "HNSC_Franken_P16", "HNSC_Franken_P20", "SKCM_Becker_P10", "SKCM_Becker_P6"), invert = T) |> 
    NormalizeData() |> 
    AddModuleScore_UCell(features = gs_list)
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  return(seu@meta.data[,c('celltype_r2', names(gs_list))])
})
qsave(list_cd8, '/bigdata/zlin/PanCancer_ICI/tables/score_cd8.qs')
list_cd8 <- qread('/bigdata/zlin/PanCancer_ICI/tables/score_cd8.qs')
df_score <- do.call(rbind, list_cd8)

Differentiation <- c("Naïve", "Activation/Effector function", "Exhaustion")
Function <- c("TCR signaling", "Cytotoxicity", "Cytokine/Cytokine receptor",
              "Chemokine/Chemokine receptor", "Anergy", "NFKB Signaling", "Stress response", "MAPK Signaling", "Adhesion",
              "IFN response")
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
col_order <- c('CD8_Naive', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
               'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13',
               'CD8_Tex_ISG15', 'CD8_Temra_CX3CR1', 'CD8_NK-like', 'CD8_Prolif')
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(df_score$celltype_r2)),
                              nrow = length(MarkerNameVector))
colnames(FunctionScoreMatrix) <- unique(df_score$celltype_r2)
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)) {
  for(ri in 1:nrow(FunctionScoreMatrix)) {
    FunctionVec <- as_tibble(df_score) |> pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[df_score$celltype_r2 == unique(df_score$celltype_r2)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, scales::rescale, to=c(0, 1)))
FunctionScoreMatrix <- FunctionScoreMatrix[, col_order]
cols <- colorRampPalette(brewer.pal('PuBu', n=100))(10)
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis))))
signatureType_row$Signature.type <- factor(signatureType_row$Signature.type, levels = c('Differentiation', 'Function','Metabolism','Apoptosis'))
pdf('/bigdata/zlin/PanCancer_ICI/figures/Functional_score_T/ht_func_cd8.pdf', height = 6, width = 8)
Heatmap(FunctionScoreMatrix, name = 'Mean \n(z-score)',
        column_title = 'CD8+ T cells',
        cluster_rows = F, cluster_columns = F,
        column_names_rot = 45,
        col = cols, 
        row_names_gp = gpar(fontsize = 10), 
        column_names_gp = gpar(fontsize = 9), 
        width = ncol(FunctionScoreMatrix)*unit(6, "mm"), 
        height = nrow(FunctionScoreMatrix)*unit(6, "mm"),
        row_title_gp = gpar(fontsize = 7, fontface = 'bold'),
        heatmap_legend_param = list(
          legend_direction = "horizontal",
          legend_width = unit(2, "cm"), at = c(0, 1),
          legend_side = 'bottom',
          title_position = "topcenter"),
        row_split = signatureType_row$Signature.type)
dev.off()

# Dataset names
datasets <- c('SKCM_Becker', 'BCC_Yost', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_IMCISION', 'HNSC_Luoma', 'HNSC_Franken', 'CRC_Li','PCa_Hawley')
# Macrophages
ifna <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_INTERFERON_ALPHA_RESPONSE') |>
  pull(gene_symbol) |>
  unique()
ifng <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_INTERFERON_GAMMA_RESPONSE') |>
  pull(gene_symbol) |>
  unique()
ifng[ifng == 'DDX58'] <- 'RIGI'
mhc1 <- msigdbr(species = "Homo sapiens", category = 'C2') |> 
  filter(gs_name  == 'REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION') |> 
  pull(gene_symbol) |> 
  unique()
mhc2 <- msigdbr(species = "Homo sapiens", category = 'C2') |> 
  filter(gs_name  == 'REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION') |> 
  pull(gene_symbol) |> 
  unique()
complement <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_COMPLEMENT') |> 
  pull(gene_symbol) |> 
  unique()
tnfa <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB') |> 
  pull(gene_symbol) |> 
  unique()
tnfa[tnfa == 'DDX58'] <- 'RIGI'
tgfb <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_TGF_BETA_SIGNALING') |> 
  pull(gene_symbol) |> 
  unique()
infla <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_INFLAMMATORY_RESPONSE') |>
  pull(gene_symbol) |>
  unique()
sig_mac <- readxl::read_xlsx('/bigdata/zlin/PanCancer_ICI/tables/1-s2.0-S0092867421000106-mmc5.xlsx', skip = 1)
sig_mac$M1[sig_mac$M1 == 'IL23'] <- 'IL23A'
sig_mac$M2[sig_mac$M2 == 'FASL'] <- 'FASLG'
gs_list <- list(`M1 Signature` = sig_mac$M1[!is.na(sig_mac$M1)],
                `M2 Signature` = sig_mac$M2[!is.na(sig_mac$M2)],
                Angiogenesis = sig_mac$Angiogenesis[!is.na(sig_mac$Angiogenesis)],
                Phagocytosis = sig_mac$Phagocytosis[!is.na(sig_mac$Phagocytosis)],
                Inflammatory = infla,
                `TNF-alpha Signaling` = tnfa,
                `TGF-beta Signaling` = tgfb,
                `IFN-alpha Response` = ifna,
                `IFN-gamma Response` = ifng,
                Complement = complement,
                `MHC I Antigen Proc.&Pres.` = mhc1,
                `MHC II Antigen Presentation` = mhc2)

list_mac <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qread(paste0('/bigdata/zlin/PanCancer_ICI/data/', dataset, '/seu_r2.qs')) |> 
    subset(subset = celltype_main %in% c('Mono','Macro')) |> 
    subset(subset = patient %in% c("CRC_Li_P28", "HNSC_Franken_P16", "HNSC_Franken_P20", "SKCM_Becker_P10", "SKCM_Becker_P6"), invert = T) |> 
    NormalizeData() |> 
    AddModuleScore_UCell(features = gs_list)
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  return(seu@meta.data[,c('celltype_r2', 'dataset', 'sample', names(gs_list))])
})
qsave(list_mac, '/bigdata/zlin/PanCancer_ICI/tables/score_mac.qs')
list_mac <- qread('/bigdata/zlin/PanCancer_ICI/tables/score_mac.qs')
df_score <- do.call(rbind, list_mac)
# mac_order <- c('Macro_IL1B','Macro_INHBA','Macro_SPP1','Macro_FN1','Macro_ISG15','Macro_TNF','Macro_LYVE1','Macro_C1QC','Macro_TREM2')
# func_order <- c("Angiogenesis", "Inflammatory" , "Phagocytosis", "M1 Signature", "M2 Signature", "TNF-alpha Signaling", "TGF-beta Signaling", "IFN-alpha Response", "IFN-gamma Response")
# df_score |> 
#   pivot_longer(names_to = 'function', values_to = 'score', cols = names(gs_list)) |> 
#   filter(`function` %in% func_order) |> 
#   mutate(celltype_r2 = factor(celltype_r2, levels = mac_order),
#          `function` = factor(`function`, levels = func_order)) |> 
#   ggplot(aes(celltype_r2, score, fill = `function`)) + 
#   geom_violin(alpha = 0.5) + 
#   geom_boxplot(alpha = 0.3, width = 0.5, coef = 6) +
#   scale_fill_met_d(name="Cross", n=9) +
#   scale_y_continuous(position = "right") +
#   facet_grid(rows = vars(`function`), scales = 'free', switch = 'y') + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = 'black', size = 8),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 10),
#         strip.text.y.left = element_text(angle = 0),
#         axis.text.y.right = element_text(size = 8),
#         plot.title = element_text(hjust = 0.5),
#         legend.position = "none") + xlab('') + ylab('') + ggtitle('Macrophages')
# ggsave('/bigdata/zlin/PanCancer_ICI/figures/Functional_score_Mac/func_mac.pdf', height = 7, width = 5)

# Heatmap
pdf('/bigdata/zlin/PanCancer_ICI/figures/Functional_score_Mac/ht_func_mac.pdf', height = 6, width = 6)
df_score |> 
  group_by(celltype_r2) |> 
  summarise(across(where(is.numeric), median, na.rm = TRUE)) |> 
  column_to_rownames(var = 'celltype_r2') |> 
  apply(2, scales::rescale, to=c(0, 1)) |> t() |>
  Heatmap(column_names_rot = 45, name = 'Mean\n(z-score)', column_title = 'Macrophages',
          heatmap_legend_param = list(
            legend_direction = "horizontal",
            legend_width = unit(2, "cm"), at = c(0, 1),
            title_position = "topcenter"),
          row_names_gp = gpar(fontsize = 10), 
          column_names_gp = gpar(fontsize = 9), 
          col = colorRampPalette(brewer.pal(10, "YlGn"))(10),
          width = 9*unit(6, "mm"), 
          height = 12*unit(6, "mm"))
dev.off()

# by dataset
funcscore_ht <- function(df, dataset, w = 4, h = 2){
  print(unique(df$dataset))
  subtypes <- df |> 
    filter(dataset == dataset) |> 
    group_by(sample, celltype_r2) |> 
    summarise(count = n(), .groups = 'drop') |> 
    tidyr::spread(celltype_r2, count, fill = 0) |> 
    select(-sample) |> 
    apply(2, function(x) sum(x >= 5)) |> 
    sapply(function(x){return(ifelse(x >= 2, 1, ifelse(length(unique(df$sample)) <= 4 & x>=1, 1, 0)))}) 
  subtypes <- names(subtypes)[subtypes == 1]
  df <- filter(df, celltype_r2 %in% subtypes, dataset == dataset)
  FunctionScoreMatrix <- matrix(0,
                                ncol = length(unique(df$celltype_r2)),
                                nrow = length(names(gs_list)))
  colnames(FunctionScoreMatrix) <- unique(df$celltype_r2)
  rownames(FunctionScoreMatrix) <- names(gs_list)
  for(ci in 1:ncol(FunctionScoreMatrix)) {
    for(ri in 1:nrow(FunctionScoreMatrix)) {
      FunctionVec <- as_tibble(df) |> pull(names(gs_list)[ri])
      fv <- median(FunctionVec[df$celltype_r2 == unique(df$celltype_r2)[ci]])
      FunctionScoreMatrix[ri, ci] <- fv
    }
  }
  FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, scales::rescale, to=c(-1, 1)))
  my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
  ht_median <- as.ggplot(pheatmap(FunctionScoreMatrix, main = unique(df$dataset), scale = 'row',
                                  # angle_col = 315,
                                  show_colnames = T,
                                  show_rownames = T,
                                  cluster_rows = F,
                                  cluster_cols = F,
                                  breaks = my.breaks,
                                  color = colorRampPalette(rev(brewer.pal(n=length(my.breaks)-2, name = 'RdBu')))(length(my.breaks)),
                                  border_color = 'white',
                                  fontsize = 8))
  p <- ht_median 
  ggsave(p, file = paste0("/bigdata/zlin/PanCancer_ICI/figures/Functional_score_Mac/Macro_FunctionScore_heatmap_", unique(df$dataset), ".pdf"), width = w, height = h)
  return(p)
}
funcscore_ht(list_mac[[11]], dataset, w = , h = 3)


# Dataset names
datasets <- c('SKCM_Becker', 'BCC_Yost', 'BRCA_Bassez1', 'BRCA_Bassez2', 'HNSC_Franken', 'CRC_Li', 'CRC_Chen', 'PCa_Hawley')
# Endothelial
leuko <- msigdbr(species = "Homo sapiens", category = 'C5') |>
  filter(gs_name  == 'GOBP_LEUKOCYTE_ADHESION_TO_ARTERIAL_ENDOTHELIAL_CELL') |>
  pull(gene_symbol) |>
  unique()
collagen <- msigdbr(species = "Homo sapiens", category = 'C2') |>
  filter(gs_name  == 'REACTOME_COLLAGEN_FORMATION') |>
  pull(gene_symbol) |>
  unique()
angio <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_ANGIOGENESIS') |> 
  pull(gene_symbol) |> 
  unique()
gs_list <- list(Angiogenesis = angio,
                `Leukocyte-endothelial Adhesion` = leuko,
                `Collagen Adhesion` = collagen)
list_endo <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qread(paste0('/bigdata/zlin/PanCancer_ICI/data/', dataset, '/seu_r2.qs')) |> 
    subset(subset = celltype_main %in% c('Endo')) |> 
    subset(subset = patient %in% c("CRC_Li_P28", "HNSC_Franken_P16", "HNSC_Franken_P20", "SKCM_Becker_P10", "SKCM_Becker_P6"), invert = T) |> 
    NormalizeData() |> 
    AddModuleScore_UCell(features = gs_list)
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  return(seu@meta.data[,c('celltype_r2', 'dataset', 'sample', names(gs_list))])
})
qsave(list_endo, '/bigdata/zlin/PanCancer_ICI/tables/score_endo.qs')
df_score <- do.call(rbind, list_mac)
df_score |> 
  group_by(celltype_r2) |> 
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) |> 
  column_to_rownames(var = 'celltype_r2') |> 
  apply(2, scales::rescale, to=c(0, 1)) |> t() |>
  Heatmap(column_names_rot = 45, name = 'Mean\n(z-score)', column_title = 'Endothelial Cells',
          heatmap_legend_param = list(
            legend_direction = "horizontal",
            legend_width = unit(2, "cm"), at = c(0, 1),
            title_position = "topcenter"),
          row_names_gp = gpar(fontsize = 10), 
          column_names_gp = gpar(fontsize = 9), 
          col = colorRampPalette(brewer.pal(10, "Blues"))(10),
          width = 9*unit(6, "mm"), 
          height = 12*unit(6, "mm"))

df_score |> 
  group_by(celltype_r2) |> 
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) |> 
  column_to_rownames(var = 'celltype_r2') |> apply(2,scales::rescale, to = c(0,1)) |> 
  t() |> 
  as.data.frame() |> 
  rownames_to_column(var='function') |> 
  ggradar(grid.label.size = 0,  # Affects the grid annotations (0%, 50%, etc.)
          axis.label.size = 7, # Afftects the names of the variables
          group.point.size = 3,
          fill = T,
          fill.alpha = 0.3, legend.position = "bottom") +
  theme(legend.position = c(1, 0),  
        legend.justification = c(1, 0),
        legend.text = element_text(size = 14),
        legend.key = element_rect(fill = NA, color = NA),
        legend.background = element_blank()) + scale_color_d3()
ggsave('/bigdata/zlin/PanCancer_ICI/figures/func_endo.pdf', height = 8, width = 8)



