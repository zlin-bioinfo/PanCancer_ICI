rm(list=ls())
pkgs <- c('Seurat', 'qs2', 'dplyr', 'tidyr',  'janitor', 'stringr','BPCells','SeuratExtend','ggplot2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1) 
options(max.print = 1000000)
setwd("/home/zlin/workspace/PanCancer_ICI")
source('scripts/Celltype_classification.R')
# final label adjusting after integration
adj_cd4t <- read.csv('tables/adjusted_cd4t.csv', row.names = 'X')
adj_cd8t <- read.csv('tables/adjusted_cd8tnk.csv', row.names = 'X')
adj_bplasma <- read.csv('tables/adjusted_bplasma.csv', row.names = 'X')
adjusted_myeloids <- read.csv('tables/adjusted_myeloids.csv', row.names = 'X')
adj_nonimmune <- read.csv('tables/adjusted_nonimmune.csv', row.names = 'X')
adj_all <- do.call(rbind, list(adj_cd4t, adj_cd8t, adj_bplasma, adjusted_myeloids, adj_nonimmune))

cohorts <- c('SCC_Yost', 'NSCLC_Liu','SKCM_Becker', 'SKCM_Plozniak', 'BCC_Yost', 
             'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 
             'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
             'NSCLC_Yan', 'NSCLC_Liu',
             'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo','HCC_Ma','RCC_Bi')

lapply(cohorts, function(cohort){
  print(cohort)
  seu <- qs_read(paste0('data/', cohort, '/seu_r2.qs2'))
  seu$celltype_adj <- 'Unchanged'
  common_cells <- intersect(seu$cell.id, adj_all$cell.id)
  seu$celltype_adj[match(common_cells, seu$cell.id)] <- adj_all$celltype_r2[match(common_cells, adj_all$cell.id)] 
  seu$celltype_r2[seu$celltype_r2 == 'PC-early_LTB'] <- "PC_IGHG"
  seu$celltype_r2[seu$celltype_adj == 'PC-trans'] <- 'PC-trans'
  seu$celltype_r2[seu$celltype_adj == 'pDC'] <- 'pDC'
  seu$celltype_main[seu$celltype_adj == 'Cycling'] <- 'Cycling'
  seu$celltype_main[seu$celltype_main == 'Cycling T/NK'] <- 'Cycling'
  seu$celltype_main[seu$celltype_r2 %in% c('GCB-cycling','PC-cycling')] <- 'Cycling'
  seu$celltype_r2[seu$celltype_main == 'Neutrophils'] <- 'Neutrophils'
  seu <- seu |> 
    subset(subset = celltype_adj %in% c("RP-high", "Doublet", "low-quality","Mixed"), invert = T) |> 
    subset(subset = celltype_r2 != 'NK_CD56hiCD16hi') |> 
    subset(subset = celltype_r2 %in% c('Melanocytes(CNA-)','Malignant(CNA+)','Epithelial(CNA-)') & celltype_bped_main %in% c('B-cells', 'CD4+ T-cells', 'CD8+ T-cells', 'DC', 'Monocytes', 'Macrophages', 'Neutrophils', 'NK cells'), invert = T) |> 
    subset(subset = celltype_r2 %in% c('Melanocytes(CNA-)','Malignant(CNA+)','Epithelial(CNA-)') & scGate_multi %in% c('Bcell', 'CD4T', 'CD8T', 'NK', 'Monocyte', 'Macrophage', 'PlasmaCell','panDC','Neutrophils','Mast'), invert = T)
  qs_save(seu, paste0('data/', cohort, '/seu_final.qs2'))
  if (!dir.exists(paste0('data/', cohort, '/seu_final/'))){
    dir.create(paste0('data/', cohort, '/seu_final/'))
  }
  burgertools::Export10X(seu, dir =paste0('data/', cohort, '/seu_final/'),
            append_reductions = NULL, gzip = F)
  print(paste0(cohort, ' done!'))
})

# Load datasets 
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 'SCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'NSCLC_Liu',
              'PCa_Hawley','RCC_Bi','HCC_Guo','HCC_Ma')
list_metadata <- lapply(datasets, function(dataset){ 
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2'))
  metadata <- seu@meta.data
  metadata <- metadata |> filter(celltype_main != 'Neutrophils')
  return(metadata)})
names(list_metadata) <- datasets
qs_save(list_metadata, 'tables/meta_list.qs2') # save as lists
list_metadata <- qs_read('tables/meta_list.qs2')
# Aggregate all metadata
# Frequency by lineage

list_metadata <- list_metadata[c('SKCM_Becker','SKCM_Plozniak', 
                                 'BCC_Yost', 'SCC_Yost',
                                 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
                                 'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
                                 'CRC_Li', 'CRC_Chen', 
                                 'NSCLC_Yan', 'NSCLC_Liu',
                                 'PCa_Hawley','RCC_Bi','HCC_Guo','HCC_Ma')]
putative_malignant <- c('Epithelial(CNA-)','Epithelial(CNA+)','Melanocytes(CNA-)','Melanocytes(CNA+)','Malignant(CNA+)')
common_cols <- Reduce(intersect, lapply(list_metadata, colnames))

meta_list <- lapply(list_metadata, function(metadata){
  print(unique(metadata$cohort))
  if (!('subtype' %in% colnames(metadata))){
    metadata$subtype <- str_split(metadata$cohort, '_', simplify = T)[,1]
  }
  metadata <- metadata[, c(common_cols, 'subtype')]
  # adjusting for cycling cells
  metadata$celltype_r2[metadata$celltype_r2 %in% c('CD4_Tcm','CD4_Treg','CD4_T-ISG','CD4_Tfh','CD4_Tstr','CD4_Tctl','CD4_Th17','CD4_T-naive',
                                                   'CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
                                                   "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
                                                   "CD8_NK-like",'MAIT','gdT') & metadata$celltype_main == 'Cycling'] <- 'Cycling T'
  metadata$celltype_r2[metadata$celltype_r2 %in% c('NK_CD56loCD16hi','NK_CD56hiCD16lo') & metadata$celltype_main == 'Cycling'] <- 'Cycling NK'
  metadata$celltype_r2[metadata$celltype_r2 %in% mye & metadata$celltype_main == 'Cycling'] <- 'Cycling myeloids'
  metadata$celltype_r2[metadata$celltype_r2 %in% nonimmune & metadata$celltype_main == 'Cycling'] <- 'Cycling non-immune'
  metadata$component <- 'celltype'
  metadata$component[metadata$celltype_r2 %in% t_nk] <- 'T_NK'
  metadata$component[metadata$celltype_r2 %in% bplasma] <- 'Bplasma'
  metadata$component[metadata$celltype_r2 %in% mye] <- 'Myeloids'
  metadata$component[metadata$celltype_r2 %in% nonimmune] <- 'Non-immune'
  metadata$component[metadata$celltype_r2 %in% putative_malignant] <- 'Malignant'
  metadata <- metadata|> 
    select(cohort, patient, sample, time_point, celltype_main, celltype_r2, interval, response, res_metric, treatment, component, Phase, subtype) |> 
    group_by(sample) |> 
    dplyr::mutate(count_sample = n()) |> 
    dplyr::mutate(count_cna = sum(celltype_r2 %in% c('Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Malignant(CNA+)')),
                  count_t = sum(celltype_r2 %in% c('CD4_Tcm','CD4_Treg','CD4_T-ISG','CD4_Tfh','CD4_Tstr','CD4_Tctl','CD4_Th17','CD4_T-naive',
                                                   'CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
                                                   "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
                                                   "CD8_NK-like",'MAIT','gdT','Cycling T')),
                  count_tcycling = sum(celltype_r2 == 'Cycling T')) |> 
    dplyr::mutate(freq_cna = count_cna/count_sample,
                  freq_t = count_t/count_sample,
                  freq_tcycling = count_tcycling/count_sample) |> 
    dplyr::mutate(non_malignant_count = sum(!celltype_r2 %in% putative_malignant)) |> 
    dplyr::mutate(freq_non_malignant = non_malignant_count/count_sample,
                  keep = case_when(non_malignant_count >= 200 ~ 'Yes',
                                   non_malignant_count < 200 ~ 'No')) |> 
    # group_by(celltype_main, .add = TRUE) |>
    # dplyr::mutate(count_main = n()) |>
    # dplyr::mutate(freq_main = count_main/count_sample) |>
    group_by(celltype_r2, sample) |>   
    dplyr::mutate(count_r2 = n()) |>  
    dplyr::mutate(freq_r2 = count_r2/count_sample) |> 
    group_by(component, sample) |>  
    dplyr::mutate(count_component = n()) |> 
    dplyr::mutate(freq_r2_comp = count_r2/count_component) |> 
    ungroup() |>
    filter(keep == 'Yes')
  return(metadata)
})

metadata <- do.call(rbind, meta_list) |> as.data.frame() 
metadata$subtype[metadata$subtype %in% c("ER+", "HER2+")] <- 'BRCA(ER/HER+)'
metadata$cohort[metadata$cohort %in% c('BCC_Yost','SCC_Yost')] <- 'BCC&SCC_Yost'
metadata$modality <- ifelse(metadata$treatment %in% c('aPD1+CTLA4', 'aPDL1+CTLA4'), 'Dual', 'Mono')
metadata <- metadata |> 
  mutate(response = case_when((cohort == 'HNSC_Luoma' & is.na(response)) ~ 'NE',
                              (cohort == 'RCC_Bi' & response %in% c('ICB_NE','NoICB')) ~ 'NE',
                              (cohort == 'RCC_Bi' & response %in% c('ICB_PD','ICB_SD')) ~ 'NR',
                              (cohort == 'RCC_Bi' & response == 'ICB_PR') ~ 'R',
                              response %in% c('PR','CR','RE') ~ 'R',
                              response %in% c('SD','PD') ~ 'NR',
                              .default = response))
metadata$subset <- 'All TME'
metadata$subset[metadata$cohort %in% c("HNSC_Luoma", "HNSC_vanderLeun", "TNBC_Zhang","HCC_Guo")] <- "CD45+sorted"
metadata$subset[metadata$cohort == "BCC&SCC_Yost" & metadata$patient %in% c('BCC_Yost_su009','BCC_Yost_su012')] <- "CD45+CD3+sorted"
metadata$subset[metadata$cohort == "BCC&SCC_Yost"] <- "CD45+CD3+sorted"
metadata$subset[metadata$cohort == "NSCLC_Liu"] <- "T cells only"
metadata$tx_status <- 'Baseline'
metadata$tx_status[metadata$time_point %in% c('On','ICI_exposed')] <- 'Treated'
write.csv(metadata, 'tables/meta_all.csv', row.names = F)
metadata <- read.csv('tables/meta_all.csv')
filter_sample <- metadata |> 
  distinct(sample, .keep_all = T) |> 
  pull(sample)
pt <- metadata |>
  distinct(sample, .keep_all = TRUE) |>
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |>
  pull(patient)

metadata <- metadata |> distinct(patient, .keep_all = T)
metadata$paired <- 'No'
metadata$paired[metadata$patient %in% pt] <- 'Yes'
metadata$interval[metadata$paired == 'No'] <- NA

metadata$x <- 1:nrow(metadata)
n <- 8
dfs <- lapply(1:n, function(i) {
  df_rep <- metadata
  df_rep$y <- i
  return(df_rep)
})
final_df <- do.call(rbind, dfs)
final_df$y_cancertype <- NA
final_df$y_cancertype[final_df$y == 8] <- 0.5
final_df$y_cohort <- NA
final_df$y_cohort[final_df$y == 7] <- 1
final_df$y_subset <- NA
final_df$y_subset[final_df$y == 6] <- 1.5
final_df$y_treatment <- NA
final_df$y_treatment[final_df$y == 5] <- 2
final_df$y_res_metric <- NA
final_df$y_res_metric[final_df$y == 4] <- 2.5
final_df$y_res <- NA
final_df$y_res[final_df$y == 3] <- 3
final_df$y_paired <- NA
final_df$y_paired[final_df$y == 2] <- 3.5
final_df$y_int <- NA
final_df$y_int[final_df$y == 1] <- 4
final_df$interval[final_df$interval >= 42] <- 42

library(MetBrewer)
library(ggnewscale)
library(ggsci)
library(ggplot2)
library(MetBrewer)
library(SeuratExtend)
ggplot(final_df) + 
  geom_tile(aes(x,y_int, width=0.4, height=0.4, fill=interval), width = 0.4) +
  scale_fill_gradient(name = 'Interval(days)', low = "yellow", high = "red", na.value = NA, limit = c(0,42), breaks = seq(0, 42, by = 21), labels = c('0','21','42+'),
                      guide = guide_colorbar(frame.colour = "black",
                                             ticks = TRUE, 
                                             title.position = "top",
                                             label.position = "bottom",
                                             barwidth = 4,
                                             barheight = 0.5, 
                                             direction = 'horizontal')) +
  new_scale_fill() +
  geom_tile(aes(x,y_paired, width=0.4, height=0.4, fill = paired)) +
  scale_fill_cosmic() +
  guides(fill = guide_legend(title = "Paired", nrow = 4)) +
  new_scale_fill() +
  geom_tile(aes(x,y_res, width=0.4, height=0.4, fill = response)) +
  scale_fill_manual(name = 'Response',values = c('R' = '#CC0C00FF',
                                                 'NR' = '#5C88DAFF',
                                                 'NE' = '#84BD00FF')) +
  guides(fill = guide_legend(nrow = 4)) +
  new_scale_fill() +
  geom_tile(aes(x,y_res_metric, width=0.4, height=0.4, fill = res_metric)) +
  scale_fill_manual(values = met.brewer('Klimt',length(unique(metadata$res_metric))), name = 'Response Metrics') +
  guides(fill = guide_legend(nrow = 4)) +
  new_scale_fill() +
  geom_tile(aes(x,y_treatment, width=0.4, height=0.4, fill = treatment)) +
  scale_fill_manual(values = met.brewer('Austria',length(unique(metadata$treatment))), name = 'Treatment') +
  guides(fill = guide_legend(nrow = 4)) +
  new_scale_fill() +
  geom_tile(aes(x,y_subset, width=0.4, height=0.4, fill = subset)) +
  scale_fill_jama(name = 'Subset') +
  guides(fill = guide_legend(nrow = 4)) +
  new_scale_fill() + 
  geom_tile(aes(x,y_cohort, width=0.4, height=0.4, fill = cohort)) +
  scale_fill_manual(values = color_pro(length(unique(metadata$cohort)), 1, sort = 'diff'), name = 'Cohort') +
  guides(fill = guide_legend(nrow = 4)) +
  new_scale_fill() +
  geom_tile(aes(x,y_cancertype, width=0.4, height=0.4, fill = subtype)) +
  scale_fill_manual(values = met.brewer('Juarez',length(unique(metadata$subtype))), name = 'Cancer Type') +
  guides(fill = guide_legend(nrow = 4)) +
  theme_void() +
  theme(legend.key.size = unit(0.4, 'cm'),  # Adjust legend key size
        legend.text = element_text(size = 8),
        legend.position = "bottom") +
  annotate("text", y = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4), x = -0.5, 
           label = c("Cancer Type", "Cohort", "Subset", "Treatment", "Response Metric", "Response", "Paired", "Interval"),
           vjust = 0.5, hjust = 1, size = 3) 
ggsave('figures/tile.pdf', height = 2, width = 16)

# Main
datasets <- c('SKCM_Plozniak', 'BCC_Yost', 'SCC_Yost', 
             'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 
             'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
             'NSCLC_Yan', 'NSCLC_Liu',
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
# seu # 1534201 cells
DimPlot(seu, group.by = 'celltype_main', reduction = 'umap.full', alpha = 0.1, 
        cols = rev(color_pro(length(unique(seu$celltype_main)), 1)), 
        label = T, label.size = 4, repel = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
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
ggsave('figures/UMAP/UMAP_main.png', height = 6, width = 8, dpi = 300)
seu$cohort[seu$cohort %in% c("BCC_Yost", "SCC_Yost")] <- "BCC&SCC_Yost"
DimPlot(seu, group.by = 'cohort', reduction = 'umap.full', alpha=0.9, 
        cols = rev(color_pro(length(unique(seu$cohort)), 1))) +
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5),
                                     legend.text=element_text(size=16)) +
  guides(color = guide_legend(ncol=1, override.aes = list(size = 4)))
ggsave('figures/UMAP/UMAP_cohort.png', height = 5, width = 7, dpi = 300)
# UMAP SKCM_Becker
seu <- qs_read(paste0('data/SKCM_Becker/seu_final.qs2'))
metadata <- read.csv('tables/meta_all.csv')
filter_sample <- metadata |> 
  distinct(sample, .keep_all = T) |> 
  pull(sample)
seu <- subset(seu, subset = sample %in% filter_sample)
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
seu <- seu |> 
  NormalizeData() |>
  FindVariableFeatures()  |>
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) |>
  RunPCA(verbose=FALSE) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony') |> 
  RunUMAP(dims = 1:20, reduction = 'harmony') |> 
  JoinLayers()
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c('CD4_T-naive','CD4_Tcm','CD4_Treg','CD4_T-ISG','CD4_Tfh','CD4_Tstr','CD4_Tctl','CD4_Th17',
                                                                          'CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
                                                                          "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
                                                                          "CD8_NK-like",'MAIT','gdT','NK_CD56loCD16hi','NK_CD56hiCD16lo')] <- 'Cycling T/NK'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein")] <- 'Endo'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c("Pericytes","SMC")] <- 'Mural'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c("Myofibroblasts", "CAF_SFRP2",
                                                                          "CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap")] <- 'CAF'
seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c('Mast','pDC','cDC1',
                                                                          'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'mregDC', 'MoDC',
                                                                          'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
                                                                          'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro-ISG',
                                                                          'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2')] <- 'Cycling myeloids'
seu$celltype_main[seu$celltype_r2 %in% mye] <- 'Myeloids'
seu$celltype_main[seu$celltype_r2 == 'pDC'] <- 'pDC'
seu$celltype_main[seu$celltype_main %in% c('Melanocytes(CNA-)','Melanocytes(CNA+)')] <- 'Melanocytes'
seu$celltype_main <- factor(seu$celltype_main, levels = c('CD4+T','CD8+T','NK','Cycling T/NK','B','Plasma','pDC','Myeloids','Endo','Mural','CAF','Melanocytes'))
DimPlot(seu, group.by = 'celltype_main', 
        cols = rev(color_pro(length(unique(seu$celltype_main)), 1)), 
        label = T, label.size = 5, repel = F, pt.size = 0.1) +
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5),
                                     legend.text=element_text(size=16)) +
  guides(color = guide_legend(ncol=1, override.aes = list(size = 4)))
ggsave('figures/UMAP/UMAP_val_cohort.png', height = 5, width = 7, dpi = 300)

# marker genes
datasets <- c('SKCM_Becker','SKCM_Plozniak', 'BCC_Yost', 'SCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'NSCLC_Yan', 'NSCLC_Liu',
              'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo','HCC_Ma','RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2'))
  seu <- seu[,!is.na(seu$celltype_r2)]
  seu <- seu |>
    subset(subset = celltype_main %in% c('Cycling T/NK','Cycling','Epithelial(CNA-)','Melanocytes(CNA-)'), invert = T)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)]) 
# Calculate the number of cells expressing each gene
expr_threshold <- 0.01 * ncol(seu)
# Subset genes expressed in at least 1% of cells
seu <- subset(seu, features = rownames(seu)[rowSums(seu@assays$RNA@counts > 0) >= expr_threshold])
seu <- JoinLayers(seu)
seu <- NormalizeData(seu)
Idents(seu) <- seu$celltype_r2
marker_cosg <- COSG::cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=100, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_all.csv', row.names = F)
gc()

# mean expression of PDCD1,2
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'RCC_Bi','HCC_Guo','HCC_Ma')
list_metadata <- lapply(datasets, function(dataset){ 
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2'))
  seu$celltype_main[seu$celltype_main == 'Cycling' & seu$celltype_r2 %in% c('Mast','pDC','cDC1',
                                                                            'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'mregDC', 'MoDC',
                                                                            'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
                                                                            'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro-ISG',
                                                                            'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2')] <- 'Cycling myeloids'
  seu <- subset(seu, subset = celltype_r2 %in% mye)
  seu <- NormalizeData(seu)
  seu <- seu[c('CD274','PDCD1LG2'),]
  avg_expr_pdl12 <- seu |> 
    AverageExpression(layer = 'data', group.by = 'sample') |> 
    data.frame(check.names = F) |> 
    t() |> data.frame(check.names = F)
  rownames(avg_expr_pdl12) <- str_replace(rownames(avg_expr_pdl12), 'RNA.','')
  rownames(avg_expr_pdl12) <- str_replace_all(rownames(avg_expr_pdl12), '-','_')
  avg_expr_pdl12$mean_pdl12 <- (avg_expr_pdl12[,1]+ avg_expr_pdl12[,2])/2
  write.csv(avg_expr_pdl12, paste0('data/', dataset, '/avg_expr_pdl12.csv'))
  return(avg_expr_pdl12)
  })


# metadata$prior <- 'No'
# metadata$prior[metadata$dataset == 'BRCA_Bassez2'] <- 'Yes'
# metadata$prior[metadata$cancertype == 'BCC' & !metadata$patient %in% c('BCC_Yost_su004')] <- 'Yes'
# metadata$prior[metadata$cancertype == 'SCC'] <- 'Yes'
# metadata$cancertype
# metadata$dataset <- as.character(metadata$dataset)
# metadata$freq_r2 <- as.numeric(metadata$freq_r2)
# metadata$count_r2 <- as.numeric(metadata$count_r2)
# metadata$interval <- as.numeric(metadata$interval)
# metadata <- metadata |> filter(!patient %in% rm_pt)
# metadata$freq_r2_comp <- as.numeric(metadata$freq_r2_comp)
# metadata$interval <- as.numeric(metadata$interval)
# metadata$int_cat <- ifelse(metadata$interval < 21, '< 21d', '>= 21d')
# Set the relative frequency for major groups to 0 for some samples
# meta_filtered <- meta_int |> mutate(freq_r2_comp = case_when((count_component < 40 & count_r2 < 3) ~ 0, .default = freq_r2_comp))
# write.csv(meta_filtered, 'tables/meta_filtered.csv', row.names = F)

# # Ro/e
# roie <- function(metadata, by_component){
#   roie_list <- lapply(unique(metadata$patient), function(p){
#     metadata <- subset(metadata, patient == p)
#     if (by_component == T){
#       roie_list_c <- list()
#       for (i in unique(metadata$component)){
#         meta_sub <- metadata |> filter(component == i)
#         observed_table <- table(meta_sub$time_point, meta_sub$celltype_r2) 
#         expected_table <- chisq.test(observed_table)$expected
#         roie <- data.frame(observed_table/expected_table)
#         observed_df <- data.frame(observed_table)
#         roie$count <- observed_df$Freq
#         roie$component <- i
#         roie_list_c[[i]] <- roie
#       }
#       roie <- do.call(rbind, roie_list_c)
#     } else {
#       observed_table <- table(metadata$time_point, metadata$celltype_r2) 
#       expected_table <- chisq.test(observed_table)$expected
#       roie <- data.frame(observed_table/expected_table)
#       observed_df <- data.frame(observed_table)
#       roie$count <- observed_df$Freq
#     }
#     names(roie)[1] <- 'time'
#     roie$patient <- p
#     roie$time <- factor(roie$time, levels = c('Pre','Post'))
#     names(roie)[2:3] <- c('celltype','ratio')
#     return(roie)
#   })
#   roie_mat <- do.call(rbind, roie_list)
#   meta <- metadata |> 
#     distinct(patient, .keep_all = T) |> 
#     select(interval, cancertype, response, res_metric, treatment, patient, time_point, prior, dataset, modality)
#   roie_mat <- left_join(roie_mat, meta, by = 'patient')
#   roie_mat <- select(roie_mat, !time_point)
#   return(roie_mat)
# }
# mat_roie <- roie(meta_int, by_component = T)
# write.csv(mat_roie, 'tables/roie.csv', row.names = F)
# 
# ROIE <- function(crosstab){
#   rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
#   rowsum.matrix[,1] <- rowSums(crosstab)
#   colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
#   colsum.matrix[1,] <- colSums(crosstab)
#   allsum <- sum(crosstab)
#   roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
#   row.names(roie) <- row.names(crosstab)
#   colnames(roie) <- colnames(crosstab)
#   return(roie)
# }



