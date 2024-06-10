#!/usr/bin/env Rscript
rm(list=ls())
pkgs <- c('Seurat','dplyr','tidyr','stringr','ggsci','qs','RColorBrewer','readxl','SingleCellExperiment','scales','pheatmap','ggplotify','ggplot2','patchwork','escape','dittoSeq','ggpubr','rstatix','UCell')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
source('/bigdata/zlin/Melanoma_meta/scripts/ggplot_theme_Publication-2.R')
options(warn = -1)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
     
# T cells
# Dataset names
dataset_names <- c("SKCM_Becker", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "BCC_Yost", "SCC_Yost", "HNSC_IMCISION", "HNSC_Luoma", "NSCLC_Liu", "CRC_Li", "PCa_Hawley")
# CD4T
gene_list <- read_xlsx('/bigdata/zlin/Melanoma_meta/tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 7, skip = 1) |> 
  lapply(function(x){return(x[!is.na(x)])})
gene_list$Exhaustion <- NULL
names(gene_list)[names(gene_list) == 'OXPHOS'] = 'Oxidative phosphorylation'
names(gene_list)[names(gene_list) == 'Lipid metabolism'] = 'Fatty acid metabolism'
gene_list$`Fatty acid metabolism`[gene_list$`Fatty acid metabolism` == 'PLA2G16'] <- 'PLAAT3'
gene_list$`Fatty acid metabolism`[gene_list$`Fatty acid metabolism` == 'ARNTL'] <- 'BMAL1'
gene_list$`Oxidative phosphorylation` <- gsub('\\.', '-', gene_list$`Oxidative phosphorylation`)

# Load and subset datasets
list_seu <- lapply(dataset_names, function(name) {
  print(name)
  seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', name, '/seu_r2.qs')) |> 
    subset(subset = celltype_main == 'CD4+T') |> 
    NormalizeData() |> 
    AddModuleScore_UCell(features = gene_list)
  for(i in 1:length(gene_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gene_list)[i], "_UCell")] <- names(gene_list)[i]
  }
  Idents(seu) <- seu$celltype_r2
  return(seu)
})

Differentiation <- c("Naïve", "Activation/Effector function")
Function <- c("TCR signaling", "Cytotoxicity", "Cytokine/Cytokine receptor",
              "Chemokine/Chemokine receptor", "Stress response", "Adhesion",
              "IFN response", "Treg signature", "Costimulatory molecules")
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
col_order <- c('CD4_Naive','CD4_Tm_CREM-','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM', 'CD4_Tm_CCL5', 
               'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
               'CD4_Treg_Early', 'CD4_Treg_ISG15', 'CD4_Treg_TNFRSF9', 
               'CD4_Th_ISG15', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif')

funcscore_ht <- function(seu, dataset){
  print(unique(seu$dataset))
  FunctionScoreMatrix <- matrix(0,
                                ncol = length(unique(seu$celltype_r2)),
                                nrow = length(MarkerNameVector))
  colnames(FunctionScoreMatrix) <- unique(seu$celltype_r2)
  rownames(FunctionScoreMatrix) <- MarkerNameVector
  for(ci in 1:ncol(FunctionScoreMatrix)) {
    for(ri in 1:nrow(FunctionScoreMatrix)) {
      FunctionVec <- as_tibble(seu@meta.data) |> pull(MarkerNameVector[ri])
      fv <- mean(FunctionVec[seu$celltype_r2 == unique(seu$celltype_r2)[ci]])
      FunctionScoreMatrix[ri, ci] <- fv
    }
  }
  FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
  FunctionScoreMatrix <- FunctionScoreMatrix[, col_order]
  my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
  my.colors <- c(
    colorRampPalette(colors = c("#008DAF", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#E22D4B"))(length(my.breaks)/2))
  signatureType_row <- data.frame(Signature.type = c(
    rep("Differentiation", length(Differentiation)),
    rep("Function", length(Function)),
    rep("Metabolism", length(Metabolism)),
    rep("Apoptosis", length(Apoptosis))))
  rownames(signatureType_row) <- MarkerNameVector
  annoCol <- c('#574141','#638200','#005D89','#B1569A')
  names(annoCol) <- unique(signatureType_row$Signature.type)
  annoCol <- list(Signature.type = annoCol)
  if (missing(dataset) == TRUE){
    dataset <- unique(seu$dataset)
  }
  ht_mean <- as.ggplot(pheatmap(FunctionScoreMatrix, main = 'Mean',scale = 'row',
                                show_colnames = T,
                                show_rownames = T,
                                annotation_row = signatureType_row,
                                annotation_colors = annoCol,
                                gaps_row = c(2, 11, 14),
                                # angle_col = 45,
                                cluster_rows = F,
                                cluster_cols = F,
                                breaks = my.breaks,
                                color = my.colors,
                                border_color = 'black',
                                fontsize = 8,
                                width = 6,
                                height = 4))
  FunctionScoreMatrix <- matrix(0,
                                ncol = length(unique(seu$celltype_r2)),
                                nrow = length(MarkerNameVector))
  colnames(FunctionScoreMatrix) <- unique(seu$celltype_r2)
  rownames(FunctionScoreMatrix) <- MarkerNameVector
  for(ci in 1:ncol(FunctionScoreMatrix)) {
    for(ri in 1:nrow(FunctionScoreMatrix)) {
      FunctionVec <- as_tibble(seu@meta.data) |> pull(MarkerNameVector[ri])
      fv <- median(FunctionVec[seu$celltype_r2 == unique(seu$celltype_r2)[ci]])
      FunctionScoreMatrix[ri, ci] <- fv
    }
  }
  FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
  FunctionScoreMatrix <- FunctionScoreMatrix[, col_order]
  ht_median <- as.ggplot(pheatmap(FunctionScoreMatrix, main = 'Median', scale = 'row',
                                  show_colnames = T,
                                  show_rownames = T,
                                  annotation_row = signatureType_row,
                                  annotation_colors = annoCol,
                                  gaps_row = c(2, 11, 14),
                                  # angle_col = 45,
                                  cluster_rows = F,
                                  cluster_cols = F,
                                  breaks = my.breaks,
                                  color = my.colors,
                                  border_color = 'black',
                                  fontsize = 8,
                                  width = 6,
                                  height = 4))
  p <- (ht_mean + ggtitle(dataset)) / ht_median + plot_layout(guides = 'collect') 
  ggsave(p, file = paste0("/bigdata/zlin/Melanoma_meta/figures/Functional_score_T/CD4_FunctionScore_heatmap_", dataset, ".pdf"), width = 8, height = 10)
  return(p)
}

list_ht_cd4 <- lapply(list_seu, function(seu){funcscore_ht(seu)})
seu <- merge(x=list_seu[[1]], y=list_seu[2:length(list_seu)]) 
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'Naïve'): which(colnames(seu@meta.data) == 'Anti.apoptosis')] <- MarkerNameVector
ht_cd4 <- seu |> funcscore_ht(dataset = 'Pan-cancer')


# CD8T
gene_list <- read_xlsx('/bigdata/zlin/Melanoma_meta/tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 5, skip = 1) |> 
  lapply(function(x){return(x[!is.na(x)])})
names(gene_list)[names(gene_list) == 'Activation:Effector function'] = 'Activation/Effector function'
names(gene_list)[names(gene_list) == 'TCR Signaling'] = 'TCR signaling'
names(gene_list)[names(gene_list) == 'IFN Response'] = 'IFN response'

# Load and subset datasets
list_seu <- lapply(dataset_names, function(name) {
  print(name)
  seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', name, '/seu_r2.qs')) |> 
    subset(subset = celltype_main == 'CD8+T') |> 
    NormalizeData() |> 
    AddModuleScore(features = gene_list, name = "FunctionScore")
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_NME1'] <- 'CD8_STR'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_ZNF683'] <- 'CD8_Trm_ZNF683'
  seu$celltype_r2[seu$celltype_r2 %in% c('CD8_NK-like_EOMES', 'CD8_NK-like_TXK')] <- 'CD8_NK-like'
  for(i in 1:length(gene_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0("FunctionScore", i)] <- names(gene_list)[i]
  }
  Idents(seu) <- seu$celltype_r2
  return(seu)
})
list_seu[[6]]$dataset <- 'BCC_SCC_Yost'
list_seu[[7]]$dataset <- 'SCC_Yost'

Differentiation <- c("Naïve", "Activation/Effector function", "Exhaustion")
Function <- c("TCR signaling", "Cytotoxicity", "Cytokine/Cytokine receptor",
              "Chemokine/Chemokine receptor", "Senescence", "Anergy", "NFKB Signaling", "Stress response", "MAPK Signaling", "Adhesion",
              "IFN response")
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
col_order <- c('CD8_Prolif', 'CD8_Naive', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
               'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13', 'CD8_Tex_OXPHOS-', 
               'CD8_ISG', 'CD8_Temra_CX3CR1', 'CD8_NK-like')

funcscore_ht <- function(seu, dataset){
  print(unique(seu$dataset))
  FunctionScoreMatrix <- matrix(0,
                                ncol = length(unique(seu$celltype_r2)),
                                nrow = length(MarkerNameVector))
  colnames(FunctionScoreMatrix) <- unique(seu$celltype_r2)
  rownames(FunctionScoreMatrix) <- MarkerNameVector
  for(ci in 1:ncol(FunctionScoreMatrix)) {
    for(ri in 1:nrow(FunctionScoreMatrix)) {
      FunctionVec <- as_tibble(seu@meta.data) |> pull(MarkerNameVector[ri])
      fv <- mean(FunctionVec[seu$celltype_r2 == unique(seu$celltype_r2)[ci]])
      FunctionScoreMatrix[ri, ci] <- fv
    }
  }
  FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
  FunctionScoreMatrix <- FunctionScoreMatrix[, col_order]
  my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
  my.colors <- c(
    colorRampPalette(colors = c("#008DAF", "white"))(length(my.breaks)/2),
    colorRampPalette(colors = c("white", "#E22D4B"))(length(my.breaks)/2))
  signatureType_row <- data.frame(Signature.type = c(
    rep("Differentiation", length(Differentiation)),
    rep("Function", length(Function)),
    rep("Metabolism", length(Metabolism)),
    rep("Apoptosis", length(Apoptosis))))
  rownames(signatureType_row) <- MarkerNameVector
  annoCol <- c('#574141','#638200','#005D89','#B1569A')
  names(annoCol) <- unique(signatureType_row$Signature.type)
  annoCol <- list(Signature.type = annoCol)
  if (missing(dataset) == TRUE){
    dataset <- unique(seu$dataset)
  }
  ht_mean <- as.ggplot(pheatmap(FunctionScoreMatrix, main = 'Mean', scale = 'row',
                                show_colnames = T,
                                show_rownames = T,
                                annotation_row = signatureType_row,
                                annotation_colors = annoCol,
                                gaps_row = c(3, 14, 17),
                                # angle_col = 45,
                                cluster_rows = F,
                                cluster_cols = F,
                                breaks = my.breaks,
                                color = my.colors,
                                border_color = 'black',
                                fontsize = 8,
                                width = 6,
                                height = 4))
  FunctionScoreMatrix <- matrix(0,
                                ncol = length(unique(seu$celltype_r2)),
                                nrow = length(MarkerNameVector))
  colnames(FunctionScoreMatrix) <- unique(seu$celltype_r2)
  rownames(FunctionScoreMatrix) <- MarkerNameVector
  for(ci in 1:ncol(FunctionScoreMatrix)) {
    for(ri in 1:nrow(FunctionScoreMatrix)) {
      FunctionVec <- as_tibble(seu@meta.data) |> pull(MarkerNameVector[ri])
      fv <- median(FunctionVec[seu$celltype_r2 == unique(seu$celltype_r2)[ci]])
      FunctionScoreMatrix[ri, ci] <- fv
    }
  }
  FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
  FunctionScoreMatrix <- FunctionScoreMatrix[, col_order]
  ht_median <- as.ggplot(pheatmap(FunctionScoreMatrix, main = 'Median', scale = 'row',
                                  show_colnames = T,
                                  show_rownames = T,
                                  annotation_row = signatureType_row,
                                  annotation_colors = annoCol,
                                  gaps_row = c(3, 14, 17),
                                  # angle_col = 45,
                                  cluster_rows = F,
                                  cluster_cols = F,
                                  breaks = my.breaks,
                                  color = my.colors,
                                  border_color = 'black',
                                  fontsize = 8,
                                  width = 6,
                                  height = 4))
  p <- (ht_mean + ggtitle(dataset)) / ht_median + plot_layout(guides = 'collect') 
  ggsave(p, file = paste0("/bigdata/zlin/Melanoma_meta/figures/Functional_score_T/CD8_FunctionScore_heatmap_", dataset, ".pdf"), width = 8, height = 10)
  return(p)
}

list_ht_cd8 <- lapply(list_seu[-c(11,12)], function(seu){funcscore_ht(seu)})
seu <- merge(x=list_seu[[1]], y=list_seu[2:length(list_seu)]) 
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'Naïve'): which(colnames(seu@meta.data) == 'Anti.apoptosis')] <- MarkerNameVector
ht_cd8 <- seu |> funcscore_ht(dataset = 'Pan-cancer')

# Myeloids
df_mac <- read_xlsx('/bigdata/zlin/Melanoma_meta/tables/1-s2.0-S0092867421000106-mmc5.xlsx', skip = 1)
df_mac$M1[df_mac$M1 == 'IL23'] <- 'IL23A'
df_mac$M2[df_mac$M2 == 'FASL'] <- 'FASLG'
gene_list <- list(M1 = df_mac$M1[!is.na(df_mac$M1)],
                  M2 = df_mac$M2[!is.na(df_mac$M2)],
                  Angiogenesis = df_mac$Angiogenesis[!is.na(df_mac$Angiogenesis)],
                  Phagocytosis = df_mac$Phagocytosis[!is.na(df_mac$Phagocytosis)])
datasets_mye <- c("SKCM_Becker", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "BCC_Yost", "HNSC_IMCISION", "HNSC_Luoma", "CRC_Li", "PCa_Hawley")
list_seu <- lapply(datasets_mye, function(dataset){
  seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
    subset(subset = celltype_main %in%c('Macro')) |> 
    NormalizeData() |> 
    AddModuleScore(features = gene_list, name = "FunctionScore")
  for(i in 1:length(gene_list)){
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0("FunctionScore", i)] <- names(gene_list)[i]
  }
  return(seu)
})
list_seu[[6]]$dataset <- 'BCC_Yost'

funcscore_ht <- function(seu, dataset, w = 4, h = 2){
  print(unique(seu$dataset))
  subtypes <- seu@meta.data |> 
    group_by(sample, cell_type = celltype_r2) |> 
    summarise(count = n(), .groups = 'drop') |> 
    spread(cell_type, count, fill = 0) |> 
    select(-sample) |> 
    apply(2, function(x) sum(x >= 5)) |> 
    sapply(function(x){return(ifelse(x >= 2, 1, ifelse(length(unique(seu$sample)) <= 4 & x>=1, 1, 0)))}) 
  subtypes <- names(subtypes)[subtypes == 1]
  seu <- subset(seu, subset = celltype_r2 %in% subtypes)
  FunctionScoreMatrix <- matrix(0,
                                ncol = length(unique(seu$celltype_r2)),
                                nrow = length(names(gene_list)))
  colnames(FunctionScoreMatrix) <- unique(seu$celltype_r2)
  rownames(FunctionScoreMatrix) <- names(gene_list)
  for(ci in 1:ncol(FunctionScoreMatrix)) {
    for(ri in 1:nrow(FunctionScoreMatrix)) {
      FunctionVec <- as_tibble(seu@meta.data) |> pull(names(gene_list)[ri])
      fv <- mean(FunctionVec[seu$celltype_r2 == unique(seu$celltype_r2)[ci]])
      FunctionScoreMatrix[ri, ci] <- fv
    }
  }
  FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
  my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
  ht_mean <- as.ggplot(pheatmap(FunctionScoreMatrix, main = 'Mean', scale = 'row',
                                # angle_col = 315,
                                show_colnames = T,
                                show_rownames = T,
                                cluster_rows = F,
                                cluster_cols = F,
                                breaks = my.breaks,
                                color = colorRampPalette(rev(brewer.pal(n=length(my.breaks), name = 'RdYlBu')))(length(my.breaks)),
                                border_color = 'white',
                                fontsize = 8))
  
  FunctionScoreMatrix <- matrix(0,
                                ncol = length(unique(seu$celltype_r2)),
                                nrow = length(names(gene_list)))
  colnames(FunctionScoreMatrix) <- unique(seu$celltype_r2)
  rownames(FunctionScoreMatrix) <- names(gene_list)
  for(ci in 1:ncol(FunctionScoreMatrix)) {
    for(ri in 1:nrow(FunctionScoreMatrix)) {
      FunctionVec <- as_tibble(seu@meta.data) |> pull(names(gene_list)[ri])
      fv <- median(FunctionVec[seu$celltype_r2 == unique(seu$celltype_r2)[ci]])
      FunctionScoreMatrix[ri, ci] <- fv
    }
  }
  FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
  my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
  ht_median <- as.ggplot(pheatmap(FunctionScoreMatrix, main = 'Median', scale = 'row',
                                  # angle_col = 315,
                                  show_colnames = T,
                                  show_rownames = T,
                                  cluster_rows = F,
                                  cluster_cols = F,
                                  breaks = my.breaks,
                                  color = colorRampPalette(rev(brewer.pal(n=length(my.breaks), name = 'RdYlBu')))(length(my.breaks)),
                                  border_color = 'white',
                                  fontsize = 8))
  p <- (ht_mean + ggtitle(unique(seu$dataset))) + ht_median
  ggsave(p, file = paste0("/bigdata/zlin/Melanoma_meta/figures/Functional_score_Mac/Macro_FunctionScore_heatmap_", unique(seu$dataset), ".pdf"), width = w, height = h)
  return(p)
}
list_ht_mac <- lapply(list_seu, function(seu){funcscore_ht(seu)})
funcscore_ht(list_seu[[2]], dataset, w = 4.5, h = 2)

seu <- merge(x=list_seu[[1]], y=list_seu[2:length(list_seu)]) 
seu$time_point <- as.character(seu$time_point)
seu$time_point[seu$time_point == 'Post' & seu$interval <= 21] <- 'On'
seu$time_point <- factor(seu$time_point, levels = c('Pre', 'On', 'Post'))
seu$int_cat <- ifelse(seu$interval > 21, 'Post', 'Early on')

# Filter samples with cell count >= 5
df_avg <- seu@meta.data |> 
  filter(response != 'NE') |> 
  group_by(sample) |> 
  mutate(cell_count = n(),
         M1_mean = mean(M1), M1_median = median(M1),
         M2_mean = mean(M2), M2_median = median(M2),
         Angiogenesis_mean = mean(Angiogenesis), Angiogenesis_median = median(Angiogenesis),
         Phagocytosis_mean = mean(Phagocytosis), Phagocytosis_median = median(Phagocytosis)) |> 
  distinct(sample, .keep_all = T) |> 
  filter(cell_count >= 5) |> ungroup()
# Keep paired samples
pt_paired <- df_avg |> 
  distinct(sample, .keep_all = T) |> 
  group_by(patient) |> 
  mutate(sample_count = n())  |> 
  select(sample_count) |> 
  filter(sample_count == 2) |> 
  pull(patient) |> unique()
df_avg <- df_avg |> filter(patient %in% pt_paired)

df_avg |> 
  ggplot(aes(x = time_point, y = Phagocytosis_median)) +
  geom_boxplot(color = 'black', outlier.shape = NA) +
  geom_line(aes(group = patient), color = "gray",linetype = "dashed") +
  geom_point(aes(color = cancertype, shape = modality), size=3, alpha = 0.6) +
  scale_color_manual(values = brewer.pal(length(unique(seu$dataset)),"Set1")) +
  facet_wrap(~response, scales = "free") + 
  theme_Publication() + xlab("") + ylab("Phagocytosis score \n (Median)") +
  # stat_pvalue_manual(stat.test, label = "T-test, p = {p}", y.position = 0.6) +
  labs(color = "Cancer type", 
       shape = "Modality")
ggsave('/bigdata/zlin/Melanoma_meta/figures/Functional_score_Mac/boxplot_phago.pdf', height = 5, width = 6)


df_avg |> 
  ggplot(aes(x = time_point, y = Angiogenesis_median)) +
  geom_boxplot(color = 'black', outlier.shape = NA) +
  geom_line(aes(group = patient), color = "gray",linetype = "dashed") +
  geom_point(aes(color = cancertype, shape = modality), size=3, alpha = 0.6) +
  scale_color_manual(values = brewer.pal(length(unique(seu$dataset)),"Set1")) +
  facet_wrap(~response, scales = "free") + 
  theme_Publication() + xlab("") + ylab("Angiogenesis score \n (Median)") +
  # stat_pvalue_manual(stat.test, label = "T-test, p = {p}", y.position = 0.25) +
  labs(color = "Interval category", 
       shape = "Modality")
ggsave('/bigdata/zlin/Melanoma_meta/figures/Functional_score_Mac/boxplot_angio.pdf', height = 5, width = 6)







