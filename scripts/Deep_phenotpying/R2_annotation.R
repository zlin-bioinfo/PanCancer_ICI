rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','SingleR','ggsci','dior','tibble','qs','BiocParallel','scGate','Matrix','SingleCellExperiment','scran','parallel','scGate','burgertools')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Loading dataset
SKCM_Becker <- qread('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/seu_r1.qs')
BRCA_Bassez1 <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r1.qs')
BRCA_Bassez2 <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/seu_r1.qs')
TNBC_Zhang <- qread('/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/seu_r1.qs')
BCC_Yost <- qread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/seu_r1.qs')
BCC_Yost1 <- qread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/seu_r1.qs') %>% 
  subset(subset = patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su012'), invert = T)
BCC_Yost1$dataset <- 'BCC_Yost'
BCC_Yost2 <- qread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/seu_r1.qs') %>% 
  subset(subset = patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su012')) 
BCC_Yost2$dataset <- 'BCC_Yost'
SCC_Yost <- qread('/bigdata/zlin/Melanoma_meta/data/SCC_Yost/seu_r1.qs')
SCC_Yost$dataset <- 'SCC_Yost'
HNSC_IMCISION <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_IMCISION/seu_r1.qs')
HNSC_Luoma <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/seu_r1.qs')
NSCLC_Liu <- qread('/bigdata/zlin/Melanoma_meta/data/NSCLC_Liu/seu_r1.qs')
CRC_Li <- qread('/bigdata/zlin/Melanoma_meta/data/CRC_Li/seu_r1.qs')
PCa_Hawley <- qread('/bigdata/zlin/Melanoma_meta/data/PCa_Hawley/seu_r1.qs')
TNBC_Shiao <- qread('/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/seu_r1.qs')
HNSC_Franken <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_Franken/seu_r1.qs')

# datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, NSCLC_Liu, CRC_Li, PCa_Hawley)

# T_NK 
ref_T <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/T_ref.qs")
ref_CD4 <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/T_CD4_ref.qs")
ref_CD8 <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/T_CD8_ref.qs")
ref_main_NK <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/NK_main_ref.qs")
ref_fine_NK <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/NK_fine_ref.qs")

# Run SingleR with mutiple reference datasets
# arguments:
# seu: seurat object 
# ref_list: reference (list)
# major: major cell types (T/NK, Myeloids)
# type_name: for file name of the predicjtion results
# datset: for file name of the prediction results
SingleRMultiRef <- function(seu, ref_list, if.subset = F, major, type_name, n_cores = 50){
  label_list <- lapply(ref_list, function(x){
    label <- x$label
    return(label)
  })
  # Create the formatted input string
  # Set the number of combinations
  num_combinations <- length(ref_list)
  # Generate random combinations
  random_combinations <- replicate(num_combinations, paste(sample(letters, 4), collapse = ""))
  input_string <- paste0(random_combinations, "=", "ref_list[[", seq_along(ref_list), "]]", collapse = ", ")
  final_input <- paste("ref = list(", input_string, ")", sep = "")
  if (if.subset == T){
    seu <- subset(seu, subset = celltype_major %in% major) 
  }
  pred <- seu %>% 
    as.SingleCellExperiment() %>% 
    logNormCounts() %>% 
    SingleR(ref = eval(parse(text = final_input)), labels = label_list, assay.type.test=1, BPPARAM=MulticoreParam(n_cores))
  print('Predictiion done!')
  qsave(pred, paste0("/bigdata/zlin/Melanoma_meta/data/", unique(seu$dataset), "/", type_name, ".qs"))
  return(pred)
}

hierSingleRMultiRef_T_NK <- function(seu, main_ref_list, ref_list_cd4, ref_list_cd8, ref_list_nk_main, ref_list_nk_fine, major = c("CD8+ T-cells","CD4+ T-cells","NK cells")){
  print(unique(seu$dataset))
  seu <- subset(seu, subset = celltype_major %in% major) 
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures()%>%
    ScaleData() %>%
    RunPCA() %>% 
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) 
  print('Gating NK cells')
  NK_KLRD1 <- gating_model(name = "NK_KLRD1", signature = c("KLRD1","CD3D-"))
  NK_NCAM1 <- gating_model(name = "NK_NCAM1", signature = c("NCAM1","CD3D-"))
  seu_1<- scGate(data = seu, model = NK_KLRD1)
  seu_2<- scGate(data = seu, model = NK_NCAM1)
  barcode_NK <- unique(c(colnames(seu_1)[seu_1$is.pure == 'Pure'],
                         colnames(seu_2)[seu_2$is.pure == 'Pure']))
  print(paste0(length(barcode_NK), ' NK cells were detected.'))
  write.csv(barcode_NK, paste0('/bigdata/zlin/Melanoma_meta/data/', unique(seu$dataset), '/barcode_NK.csv'))
  print('Gating gdT cells')
  gdT <- gating_model(name = "gdT", signature = c("CD3D","TRGC2","TRDV2","TRGV9","TRGV10","TRDC","CD8A-","CD4-"))
  seu_gdT<- scGate(data = seu, model = gdT)
  barcode_gdT <- colnames(seu_gdT)[seu_gdT$is.pure == 'Pure']
  # barcode_gdT <- setdiff(barcode_gdT, barcode_NK)
  print(paste0(length(barcode_gdT), ' gdT cells were detected'))
  write.csv(barcode_gdT, paste0('/bigdata/zlin/Melanoma_meta/data/', unique(seu$dataset), '/barcode_gdT.csv'))
  scGate_models_DB <- get_scGateDB()
  seu_cd4cd8 <- seu[, !colnames(seu) %in% c(barcode_NK, barcode_gdT)]
  print('Gating CD4+T cells')
  seu_cd4 <- scGate(seu_cd4cd8, model = scGate_models_DB$human$generic$CD4T)
  barcode_cd4 <- colnames(seu_cd4)[seu_cd4$is.pure == 'Pure']
  print('Gating CD8+T cells')
  seu_cd8 <- scGate(seu_cd4cd8, model = scGate_models_DB$human$generic$CD8T)
  barcode_cd8 <- colnames(seu_cd8)[seu_cd8$is.pure == 'Pure']
  dual <- intersect(barcode_cd4, barcode_cd8)
  barcode_cd4 <- setdiff(barcode_cd4, dual)
  barcode_cd8 <- setdiff(barcode_cd8, dual)
  print('CD4/CD8')
  # Unresolved CD4/CD8 cells
  pred_T <- SingleRMultiRef(seu[, !colnames(seu) %in% c(barcode_NK, barcode_gdT, barcode_cd4, barcode_cd8)], ref_list = main_ref_list, type_name = 'T')
  barcode_CD4 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD4'], barcode_cd4)
  barcode_CD8 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD8'], barcode_cd8)
  print('CD4')
  pred_CD4 <- SingleRMultiRef(seu = seu[, barcode_CD4], ref_list = ref_list_cd4, type_name = 'CD4')
  print('CD8')
  pred_CD8 <- SingleRMultiRef(seu = seu[, barcode_CD8], ref_list = ref_list_cd8, type_name = 'CD8')
  print('NK main')
  pred_NK_main <- SingleRMultiRef(seu = seu[, barcode_NK], ref_list = ref_list_nk_main, type_name = 'NK_main')
  print('NK fine')
  pred_NK_fine <- SingleRMultiRef(seu = seu[, barcode_NK], ref_list = ref_list_nk_fine, type_name = 'NK_fine')
}
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, NSCLC_Liu, CRC_Li, PCa_Hawley)
mclapply(datasets, function(seu){
  hierSingleRMultiRef_T_NK(seu, main_ref_list = ref_T, 
                          ref_list_cd4 = ref_CD4, ref_list_cd8 = ref_CD8, 
                          ref_list_nk_main = ref_main_NK, ref_list_nk_fine = ref_fine_NK,
                          major = c("CD8+ T-cells","CD4+ T-cells","NK cells"))}, mc.cores = 100)
hierSingleRMultiRef_T_NK(HNSC_Franken, main_ref_list = ref_T, 
                         ref_list_cd4 = ref_CD4, ref_list_cd8 = ref_CD8, 
                         ref_list_nk_main = ref_main_NK, ref_list_nk_fine = ref_fine_NK,
                         major = c("CD8+ T-cells","CD4+ T-cells","NK cells"))

# T cells only (CD45+CD3+ sorted)
hierSingleRMultiRef_T <- function(seu, main_ref_list, ref_list_cd4, ref_list_cd8){
  print(unique(seu$dataset))
  seu <- NormalizeData(seu) %>%
    FindVariableFeatures()%>%
    ScaleData() %>%
    RunPCA() %>% 
    RunUMAP(dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1)
  print('Gating gdT cells')
  gdT <- gating_model(name = "gdT", signature = c("CD3D","TRGC2","TRDV2","TRGV9","TRGV10","TRDC","CD8A-","CD4-"))
  seu_gdT<- scGate(data = seu, model = gdT)
  barcode_gdT <- colnames(seu_gdT)[seu_gdT$is.pure == 'Pure']
  # barcode_gdT <- setdiff(barcode_gdT, barcode_NK)
  print(paste0(length(barcode_gdT), ' gdT cells were detected'))
  write.csv(barcode_gdT, paste0('/bigdata/zlin/Melanoma_meta/data/', unique(seu$dataset), '/barcode_gdT_sorted.csv'))
  scGate_models_DB <- get_scGateDB()
  seu_cd4cd8 <- seu[,!colnames(seu) %in% barcode_gdT]
  print('Gating CD4+T cells')
  seu_cd4 <- scGate(seu_cd4cd8, model = scGate_models_DB$human$generic$CD4T)
  barcode_cd4 <- colnames(seu_cd4)[seu_cd4$is.pure == 'Pure']
  print('Gating CD8+T cells')
  seu_cd8 <- scGate(seu_cd4cd8, model = scGate_models_DB$human$generic$CD8T)
  barcode_cd8 <- colnames(seu_cd8)[seu_cd8$is.pure == 'Pure']
  dual <- intersect(barcode_cd4, barcode_cd8)
  barcode_cd4 <- setdiff(barcode_cd4, dual)
  barcode_cd8 <- setdiff(barcode_cd8, dual)
  print('CD4/CD8')
  # Unresolved CD4/CD8 cells
  pred_T <- SingleRMultiRef(seu[, !colnames(seu) %in% c(barcode_gdT, barcode_cd4, barcode_cd8)], ref_list = main_ref_list, type_name = 'T') # For CD45CD3 sorted (SCC_Yost)
  barcode_CD4 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD4'], barcode_cd4)
  barcode_CD8 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD8'], barcode_cd8)
  print('CD4')
  pred_CD4 <- SingleRMultiRef(seu = seu[, barcode_CD4], ref_list = ref_list_cd4, type_name = 'sorted_CD4')
  print('CD8')
  pred_CD8 <- SingleRMultiRef(seu = seu[, barcode_CD8], ref_list = ref_list_cd8, type_name = 'sorted_CD8')
}
hierSingleRMultiRef_T(BCC_Yost2, main_ref_list = ref_T, ref_list_cd4 = ref_CD4, ref_list_cd8 = ref_CD8)
hierSingleRMultiRef_T(SCC_Yost, main_ref_list = ref_T, ref_list_cd4 = ref_CD4, ref_list_cd8 = ref_CD8)
# Myeloids
ref_Myeloids <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/Myeloids_major_ref.qs")
ref_cdc2 <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/Myeloids_cdc2_ref.qs")
ref_mono <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/Myeloids_mono_ref.qs")
ref_macro <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/Myeloids_macro_ref.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, TNBC_Shiao, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, CRC_Li, PCa_Hawley)

hierSingleRMultiRef_Myeloids <- function(seu, main_ref_list, ref_list_cdc2, ref_list_mono, ref_list_macro){
  pred <- SingleRMultiRef(seu = seu, ref_list = main_ref_list, if.subset = T, major = c('DC', 'Monocytes', 'Macrophages'), type_name = 'Myeloids_major')
  print('Major done!')
  barcode_cdc2 <- rownames(pred)[pred$pruned.labels == 'cDC2']
  barcode_mono <- rownames(pred)[pred$pruned.labels == 'Mono']
  barcode_macro <- rownames(pred)[pred$pruned.labels == 'Macro']
  print('cDC2s')
  pred_cdc2 <- SingleRMultiRef(seu = seu[, barcode_cdc2], ref_list = ref_list_cdc2, type_name = 'cDC2')
  print('Monocytes')
  pred_mono <- SingleRMultiRef(seu = seu[, barcode_mono], ref_list = ref_list_mono, type_name = 'Mono')
  print('Macrophages')
  pred_macro <- SingleRMultiRef(seu = seu[, barcode_macro], ref_list = ref_list_macro, type_name = 'Macro')
}
mclapply(datasets, function(seu){
  hierSingleRMultiRef_Myeloids(seu, main_ref_list = ref_Myeloids, ref_list_cdc2 = ref_cdc2, 
                               ref_list_mono = ref_mono, ref_list_macro = ref_macro)}, mc.cores = 100)
hierSingleRMultiRef_Myeloids(HNSC_Franken, main_ref_list = ref_Myeloids, 
                             ref_list_mono = ref_mono, ref_list_macro = ref_macro, ref_list_cdc2 = ref_cdc2)

# run SingleR with one reference dataset
SingleROneRef <- function(seu, ref, label, major, type_name, n_cores = 50){
  pred <- seu %>% 
    subset(subset = celltype_major %in% major) %>% 
    as.SingleCellExperiment() %>% 
    logNormCounts() %>% 
    SingleR(ref = ref, labels = label, assay.type.test=1, BPPARAM=MulticoreParam(50))
  qsave(pred, paste0("/bigdata/zlin/Melanoma_meta/data/", unique(seu$dataset), "/", type_name, ".qs"))
}

# B cells 
ref_B <- MonacoImmuneData()
# bped <- celldex::BlueprintEncodeData()
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, CRC_Li, PCa_Hawley)
ref_B <- ref_B[, str_detect(ref_B$label.fine, 'B cells') | str_detect(ref_B$label.fine, 'Plasmablasts')]
lapply(datasets, function(seu){
  print(unique(seu$dataset))
  SingleROneRef(seu, ref = ref_B, label = ref_B$label.fine, major = 'B-cells', type_name = 'B')})
SingleROneRef(HNSC_Franken, ref = ref_B, label = ref_B$label.fine, major = 'B-cells', type_name = 'B')
# pan-cancer B
ref_Bplasma <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/Pan_B_major.qs")
ref_b <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/B.qs")
ref_gcb <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/GCB.qs")
ref_plasma <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/Plasma.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, TNBC_Shiao, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, HNSC_Franken, CRC_Li, PCa_Hawley)
hierSingleRMultiRef_Bplasma <- function(seu, main_ref_list, ref_list_b, ref_list_gcb, ref_list_plasma){
  pred <- SingleRMultiRef(seu = seu, ref_list = main_ref_list, if.subset = T, major = c('B-cells', 'Plasma'), type_name = 'Bplasma_major')
  print('Major done!')
  pred_b <- rownames(pred)[pred$pruned.labels == 'B']
  pred_gcb <- rownames(pred)[pred$pruned.labels == 'GCB']
  pred_plasma <- rownames(pred)[pred$pruned.labels == 'Plasma']
  print('B')
  pred_b <- SingleRMultiRef(seu = seu[, pred_b], ref_list = ref_list_b, type_name = 'pan-B')
  if (length(pred_gcb)>0){
    print('GCB')
    pred_gcb <- SingleRMultiRef(seu = seu[, pred_gcb], ref_list = ref_list_gcb, type_name = 'GCB')
  }
  if (length(pred_plasma)>0){
  print('Plasma')
  pred_plasma <- SingleRMultiRef(seu = seu[, pred_plasma], ref_list = ref_list_plasma, type_name = 'Plasma')
  }
}
mclapply(datasets, function(seu){
  hierSingleRMultiRef_Bplasma(seu, main_ref_list = ref_Bplasma, ref_list_b = ref_b, 
                               ref_list_gcb = ref_gcb, ref_list_plasma = ref_plasma)}, mc.cores = 100)
hierSingleRMultiRef_Bplasma(PCa_Hawley, main_ref_list = ref_Bplasma, ref_list_b = ref_b, 
                             ref_list_gcb = ref_gcb, ref_list_plasma = ref_plasma)
# Endothelial cells
ref_Endo <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/Endo_ref.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, BCC_Yost1, CRC_Li, PCa_Hawley)
lapply(datasets, function(seu){SingleROneRef(seu, ref = ref_Endo, label = ref_Endo$label, major = 'Endothelial cells', type_name = 'Endo')})
SingleROneRef(HNSC_Franken, ref = ref_Endo, label = ref_Endo$label, major = 'Endothelial cells', type_name = 'Endo')

# CAF
ref_CAF <- qread("/bigdata/zlin/Melanoma_meta/data/Ref_SingleR/CAF_ref.qs")
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, BCC_Yost1, CRC_Li, PCa_Hawley)
lapply(datasets, function(seu){
  print(unique(seu$dataset))
  SingleRMultiRef(seu, ref_list = ref_CAF, major = c('Fibroblasts','Myocytes'), if.subset = T, type_name = 'CAF')})
SingleRMultiRef(HNSC_Franken, ref_list = ref_CAF, major = c('Fibroblasts','Myocytes'), if.subset = T, type_name = 'CAF')

dataset <- 'SKCM_Becker'
cell_subtype <- c('CD4','CD8','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','Endo','CAF')
SKCM_Becker$celltype_major <- as.character(SKCM_Becker$celltype_major)
SKCM_Becker <- subset(SKCM_Becker, subset = celltype_major %in% c('Melanoma', 'Epithelial cells'), invert = T)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2','Macro', 'Mono'))
SKCM_Becker$celltype_r2 <- pred_concat$celltype_r2[match(colnames(SKCM_Becker), pred_concat$cell)]
# SKCM_Becker$celltype_r2[SKCM_Becker$celltype_major == 'Plasma'] <- 'Plasma'
# SKCM_Becker$celltype_r2[SKCM_Becker$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- SKCM_Becker |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |>
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |>
#   RunPCA(verbose=FALSE) |>
#   RunUMAP(dims = 1:20) |>
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 5)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# SKCM_Becker$celltype_r2[colnames(SKCM_Becker) %in% colnames(seu)[seu$seurat_clusters == 0]] <- 'Plasmablast'
seu <- SKCM_Becker |>
  subset(subset = celltype_major %in% c("Plasma","B-cells")) |>
  NormalizeData() |>
  FindVariableFeatures()
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
SKCM_Becker$celltype_r2[SKCM_Becker$celltype_major == 'pDC'] <- 'pDC'
SKCM_Becker$celltype_major[SKCM_Becker$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
SKCM_Becker$celltype_major[SKCM_Becker$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(SKCM_Becker$celltype_r2, SKCM_Becker$celltype_major, useNA = 'ifany')
SKCM_Becker <- SKCM_Becker[,!is.na(SKCM_Becker$celltype_r2)]
SKCM_Becker$celltype_main <- 'celltype'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_r2, 'CD4')] <- 'CD4+T'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_r2, 'CD8')] <- 'CD8+T'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_r2, 'gdT')] <- 'CD8+T'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_major, 'NK')] <- 'NK'
SKCM_Becker$celltype_main[SKCM_Becker$celltype_r2 %in% bcells] <- 'B'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_r2, 'Plasma')] <- 'Plasma'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_major, 'pDC')] <- 'pDC'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_r2, 'cDC')] <- 'cDC'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_r2, 'Mono')] <- 'Mono'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_r2, 'Macro')] <- 'Macro'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_r2, 'EC_')] <- 'Endo'
SKCM_Becker$celltype_main[str_detect(SKCM_Becker$celltype_r2, 'CAF') | str_detect(SKCM_Becker$celltype_r2, 'Myofibroblast')] <- 'CAF'
table(SKCM_Becker$celltype_main, SKCM_Becker$celltype_major, useNA = 'ifany')
SKCM_Becker <- subset(SKCM_Becker, subset = celltype_r2 == 'CD56highCD16high', invert = T)
# keep only site-matched samples
SKCM_Becker <- subset(SKCM_Becker, subset = patient %in% c('SKCM_Becker_P2', 'SKCM_Becker_P5', 'SKCM_Becker_P6','SKCM_Becker_P7','SKCM_Becker_P9','SKCM_Becker_P10','SKCM_Becker_P11','SKCM_Becker_P12'))
qsave(SKCM_Becker, '/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/seu_r2.qs')
SKCM_Becker <- qread('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/seu_r2.qs')

dataset <- 'BRCA_Bassez1'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','Endo','CAF','pan-B','Plasma','GCB')
BRCA_Bassez1$celltype_major <- as.character(BRCA_Bassez1$celltype_major)
BRCA_Bassez1 <- subset(BRCA_Bassez1, subset = celltype_major %in% c('Epithelial cells'), invert = T)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2','Macro', 'Mono'))
BRCA_Bassez1$celltype_r2 <- pred_concat$celltype_r2[match(colnames(BRCA_Bassez1), pred_concat$cell)]
# BRCA_Bassez1$celltype_r2[BRCA_Bassez1$celltype_major == 'Plasma'] <- 'Plasma'
# BRCA_Bassez1$celltype_r2[BRCA_Bassez1$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- BRCA_Bassez1 |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |> 
#   RunPCA(verbose=FALSE) |> 
#   RunUMAP(dims = 1:20) |> 
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 3)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# BRCA_Bassez1$celltype_r2[colnames(BRCA_Bassez1) %in% colnames(seu)[seu$seurat_clusters == 13]] <- 'Plasmablast'
seu <- BRCA_Bassez1 |>
  subset(subset = celltype_major %in% c("Plasma",'B-cells')) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() 
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
seu <- seu[,!is.na(seu$celltype_r2)]
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
BRCA_Bassez1$celltype_r2[BRCA_Bassez1$celltype_major == 'pDC'] <- 'pDC'
BRCA_Bassez1$celltype_r2[BRCA_Bassez1$celltype_major == 'Mast'] <- 'Mast'
BRCA_Bassez1$celltype_major[BRCA_Bassez1$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
BRCA_Bassez1$celltype_major[BRCA_Bassez1$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(BRCA_Bassez1$celltype_r2, BRCA_Bassez1$celltype_major, useNA = 'ifany')
BRCA_Bassez1 <- BRCA_Bassez1[,!is.na(BRCA_Bassez1$celltype_r2)]
BRCA_Bassez1$celltype_main <- 'celltype'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_r2, 'CD4')] <- 'CD4+T'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_r2, 'CD8')] <- 'CD8+T'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_r2, 'gdT')] <- 'CD8+T'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_major, 'NK')] <- 'NK'
BRCA_Bassez1$celltype_main[BRCA_Bassez1$celltype_r2 %in% bcells] <- 'B'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_r2, 'Plasma')] <- 'Plasma'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_major, 'pDC')] <- 'pDC'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_r2, 'cDC')] <- 'cDC'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_r2, 'Mono')] <- 'Mono'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_r2, 'Macro')] <- 'Macro'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_r2, 'EC_')] <- 'Endo'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_r2, 'CAF') | str_detect(BRCA_Bassez1$celltype_r2, 'Myofibroblast')] <- 'CAF'
BRCA_Bassez1$celltype_main[str_detect(BRCA_Bassez1$celltype_major, 'Mast')] <- 'Mast'
table(BRCA_Bassez1$celltype_main, BRCA_Bassez1$celltype_major, useNA = 'ifany')
BRCA_Bassez1 <- subset(BRCA_Bassez1, subset = celltype_r2 == 'CD56highCD16high', invert = T)
qsave(BRCA_Bassez1, '/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r2.qs')
BRCA_Bassez1 <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r2.qs')

dataset <- 'BRCA_Bassez2'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','Endo','CAF','pan-B','Plasma','GCB')
BRCA_Bassez2$celltype_major <- as.character(BRCA_Bassez2$celltype_major)
BRCA_Bassez2 <- subset(BRCA_Bassez2, subset = celltype_major %in% c('Epithelial cells'), invert = T)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono'))
BRCA_Bassez2$celltype_r2 <- pred_concat$celltype_r2[match(colnames(BRCA_Bassez2), pred_concat$cell)]
# BRCA_Bassez2$celltype_r2[BRCA_Bassez2$celltype_major == 'Plasma'] <- 'Plasma'
# BRCA_Bassez2$celltype_r2[BRCA_Bassez2$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- BRCA_Bassez2 |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |> 
#   RunPCA(verbose=FALSE) |> 
#   RunUMAP(dims = 1:20) |> 
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 3)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# BRCA_Bassez2$celltype_r2[colnames(BRCA_Bassez2) %in% colnames(seu)[seu$seurat_clusters == 6]] <- 'Plasmablast'
seu <- BRCA_Bassez2 |>
  subset(subset = celltype_major %in% c("Plasma",'B-cells')) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() 
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
BRCA_Bassez2$celltype_r2[BRCA_Bassez2$celltype_major == 'pDC'] <- 'pDC'
BRCA_Bassez2$celltype_r2[BRCA_Bassez2$celltype_major == 'Mast'] <- 'Mast'
BRCA_Bassez2$celltype_major[BRCA_Bassez2$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
BRCA_Bassez2$celltype_major[BRCA_Bassez2$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(BRCA_Bassez2$celltype_r2, BRCA_Bassez2$celltype_major, useNA = 'ifany')
BRCA_Bassez2 <- BRCA_Bassez2[,!is.na(BRCA_Bassez2$celltype_r2)]
BRCA_Bassez2$celltype_main <- 'celltype'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_r2, 'CD4')] <- 'CD4+T'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_r2, 'CD8')] <- 'CD8+T'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_r2, 'gdT')] <- 'CD8+T'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_major, 'NK')] <- 'NK'
BRCA_Bassez2$celltype_main[BRCA_Bassez2$celltype_r2 %in% bcells] <- 'B'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_r2, 'Plasma')] <- 'Plasma'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_major, 'pDC')] <- 'pDC'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_r2, 'cDC')] <- 'cDC'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_r2, 'Mono')] <- 'Mono'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_r2, 'Macro')] <- 'Macro'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_r2, 'EC_')] <- 'Endo'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_r2, 'CAF') | str_detect(BRCA_Bassez2$celltype_r2, 'Myofibroblast')] <- 'CAF'
BRCA_Bassez2$celltype_main[str_detect(BRCA_Bassez2$celltype_major, 'Mast')] <- 'Mast'
table(BRCA_Bassez2$celltype_main, BRCA_Bassez2$celltype_major, useNA = 'ifany')
BRCA_Bassez2 <- subset(BRCA_Bassez2, subset = celltype_r2 %in% 'CD56highCD16high', invert = T)
qsave(BRCA_Bassez2, '/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/seu_r2.qs')
BRCA_Bassez2 <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/seu_r2.qs')

dataset <- 'TNBC_Zhang'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-B','Plasma','GCB')
TNBC_Zhang$celltype_major <- as.character(TNBC_Zhang$celltype_major)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
TNBC_Zhang$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono'))
TNBC_Zhang$celltype_r2 <- pred_concat$celltype_r2[match(colnames(TNBC_Zhang), pred_concat$cell)]
# TNBC_Zhang$celltype_r2[TNBC_Zhang$celltype_major == 'Plasma'] <- 'Plasma'
# TNBC_Zhang$celltype_r2[TNBC_Zhang$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- TNBC_Zhang |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() 
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# TNBC_Zhang$celltype_r2[colnames(TNBC_Zhang) %in% colnames(seu)[seu$seurat_clusters == 10]] <- 'Plasmablast'
# TNBC_Zhang$celltype_r2[colnames(TNBC_Zhang) %in% colnames(seu)[seu$seurat_clusters == 20]] <- 'GC-like B'
seu <- TNBC_Zhang |>
  subset(subset = celltype_major %in% c("Plasma",'B-cells')) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() 
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
seu <- seu[,!is.na(seu$celltype_r2)]
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
TNBC_Zhang$celltype_r2[TNBC_Zhang$celltype_major == 'pDC'] <- 'pDC'
TNBC_Zhang$celltype_r2[TNBC_Zhang$celltype_major == 'Mast'] <- 'Mast'
TNBC_Zhang$celltype_major[TNBC_Zhang$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
TNBC_Zhang$celltype_major[TNBC_Zhang$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(TNBC_Zhang$celltype_r2, TNBC_Zhang$celltype_major, useNA = 'ifany')
TNBC_Zhang <- TNBC_Zhang[,!is.na(TNBC_Zhang$celltype_r2)]
TNBC_Zhang$celltype_main <- 'celltype'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_r2, 'CD4')] <- 'CD4+T'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_r2, 'CD8')] <- 'CD8+T'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_r2, 'gdT')] <- 'CD8+T'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_major, 'NK')] <- 'NK'
TNBC_Zhang$celltype_main[TNBC_Zhang$celltype_r2 %in% bcells] <- 'B'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_r2, 'Plasma')] <- 'Plasma'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_major, 'pDC')] <- 'pDC'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_r2, 'cDC')] <- 'cDC'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_r2, 'Mono')] <- 'Mono'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_r2, 'Macro')] <- 'Macro'
TNBC_Zhang$celltype_main[str_detect(TNBC_Zhang$celltype_major, 'Mast')] <- 'Mast'
table(TNBC_Zhang$celltype_main, TNBC_Zhang$celltype_major, useNA = 'ifany')
TNBC_Zhang <- subset(TNBC_Zhang, subset = celltype_r2 == 'CD56highCD16high', invert = T)
qsave(TNBC_Zhang, '/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/seu_r2.qs')
TNBC_Zhang <- qread('/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/seu_r2.qs')

dataset <- 'BCC_Yost'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','Endo','CAF','pan-B','Plasma','GCB')
BCC_Yost$celltype_major <- as.character(BCC_Yost$celltype_major)
BCC_Yost <- subset(BCC_Yost, subset = celltype_major %in% c('Epithelial cells', 'Melanocytes'), invert = T)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono'))
# su009 su012
cell_subtype <- c('sorted_CD4','sorted_CD8')
pred_list2 <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list2) <- cell_subtype
pred_concat2 <- do.call(rbind, lapply(pred_list2, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT2 <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT_sorted.csv')) 
pred_gdT2 <- data.frame(cell = barcode_gdT2[,'x'], celltype_r2 = 'gdT')
pred_concat <- purrr::reduce(list(pred_concat, pred_gdT2, pred_concat2), rbind)
BCC_Yost$celltype_r2 <- pred_concat$celltype_r2[match(colnames(BCC_Yost), pred_concat$cell)]
# BCC_Yost$celltype_r2[BCC_Yost$celltype_major == 'Plasma'] <- 'Plasma'
# BCC_Yost$celltype_r2[BCC_Yost$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- BCC_Yost |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |> 
#   RunPCA(verbose=FALSE) |> 
#   RunUMAP(dims = 1:20) |> 
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 3)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# BCC_Yost$celltype_r2[colnames(BCC_Yost) %in% colnames(seu)[seu$seurat_clusters == 16]] <- 'Plasmablast'
seu <- BCC_Yost |>
  subset(subset = celltype_major %in% c("Plasma","B-cells")) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData()
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
seu <- seu[,!is.na(seu$celltype_r2)]
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
BCC_Yost$celltype_r2[BCC_Yost$celltype_major == 'pDC'] <- 'pDC'
BCC_Yost$celltype_r2[BCC_Yost$celltype_major == 'Mast'] <- 'Mast'
BCC_Yost$celltype_major[BCC_Yost$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), unique(pred_list2[[1]]$pruned.labels), unique(pred_list2[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
BCC_Yost$celltype_major[BCC_Yost$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(BCC_Yost$celltype_r2, BCC_Yost$celltype_major, useNA = 'ifany')
BCC_Yost <- BCC_Yost[,!is.na(BCC_Yost$celltype_r2)]
BCC_Yost$celltype_main <- 'celltype'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'CD4')] <- 'CD4+T'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'CD8')] <- 'CD8+T'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'gdT')] <- 'CD8+T'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_major, 'NK')] <- 'NK'
BCC_Yost$celltype_main[BCC_Yost$celltype_r2 %in% bcells] <- 'B'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'Plasma') & !BCC_Yost$patient %in% c('su009', 'su012')] <- 'Plasma'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_major, 'pDC')] <- 'pDC'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'cDC')] <- 'cDC'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'Mono')] <- 'Mono'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'Macro')] <- 'Macro'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'EC_')] <- 'Endo'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_major, 'Mast')] <- 'Mast'
BCC_Yost$celltype_main[str_detect(BCC_Yost$celltype_r2, 'CAF') | str_detect(BCC_Yost$celltype_r2, 'Myofibroblast')] <- 'CAF'
table(BCC_Yost$patient, BCC_Yost$celltype_main, useNA = 'ifany')
BCC_Yost <- BCC_Yost[,!BCC_Yost$celltype_main == 'celltype']
BCC_Yost <- subset(BCC_Yost, subset = celltype_r2 == 'CD56highCD16high', invert = T)
BCC_Yost$dataset <- 'BCC_Yost'
qsave(BCC_Yost, '/bigdata/zlin/Melanoma_meta/data/BCC_Yost/seu_r2.qs')
BCC_Yost <- qread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/seu_r2.qs')

# SCC_Yost
dataset <- 'SCC_Yost'
cell_subtype <- c('CD4','CD8')
SCC_Yost$celltype_major <- as.character(SCC_Yost$celltype_major)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/sorted_', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
SCC_Yost$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT_sorted.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
mapping <- setNames(pred_concat$celltype_r2, pred_concat$cell)
SCC_Yost$celltype_r2 <- mapping[colnames(SCC_Yost)]
table(SCC_Yost$celltype_r2, SCC_Yost$celltype_major, useNA = 'ifany')
SCC_Yost <- SCC_Yost[,!is.na(SCC_Yost$celltype_r2)]
SCC_Yost$celltype_main <- 'celltype'
SCC_Yost$celltype_main[str_detect(SCC_Yost$celltype_r2, 'CD4')] <- 'CD4+T'
SCC_Yost$celltype_main[str_detect(SCC_Yost$celltype_r2, 'CD8')] <- 'CD8+T'
SCC_Yost$celltype_main[str_detect(SCC_Yost$celltype_r2, 'gdT')] <- 'CD8+T'
table(SCC_Yost$celltype_r2, SCC_Yost$celltype_main, useNA = 'ifany')
qsave(SCC_Yost, '/bigdata/zlin/Melanoma_meta/data/SCC_Yost/seu_r2.qs')

dataset <- 'HNSC_IMCISION'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-B','Plasma','GCB')
HNSC_IMCISION$celltype_major <- as.character(HNSC_IMCISION$celltype_major)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
HNSC_IMCISION$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono'))
HNSC_IMCISION$celltype_r2 <- pred_concat$celltype_r2[match(colnames(HNSC_IMCISION), pred_concat$cell)]
# HNSC_IMCISION$celltype_r2[HNSC_IMCISION$celltype_major == 'Plasma'] <- 'Plasma'
# HNSC_IMCISION$celltype_r2[HNSC_IMCISION$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- HNSC_IMCISION |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |> 
#   RunPCA(verbose=FALSE) |> 
#   RunUMAP(dims = 1:20) |> 
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 3)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# HNSC_IMCISION$celltype_r2[colnames(HNSC_IMCISION) %in% colnames(seu)[seu$seurat_clusters == 11]] <- 'Plasmablast'
seu <- HNSC_IMCISION |>
  subset(subset = celltype_major %in% c("Plasma","B-cells")) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData()
seu <- seu[,!is.na(seu$celltype_r2)]
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
seu <- seu[,!is.na(seu$celltype_r2)]
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
HNSC_IMCISION$celltype_r2[HNSC_IMCISION$celltype_major == 'pDC'] <- 'pDC'
HNSC_IMCISION$celltype_r2[HNSC_IMCISION$celltype_major == 'Mast'] <- 'Mast'
HNSC_IMCISION$celltype_major[HNSC_IMCISION$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
HNSC_IMCISION$celltype_major[HNSC_IMCISION$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(HNSC_IMCISION$celltype_r2, HNSC_IMCISION$celltype_major, useNA = 'ifany')
HNSC_IMCISION <- HNSC_IMCISION[,!is.na(HNSC_IMCISION$celltype_r2)]
HNSC_IMCISION$celltype_main <- 'celltype'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_r2, 'CD4')] <- 'CD4+T'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_r2, 'CD8')] <- 'CD8+T'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_r2, 'gdT')] <- 'CD8+T'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_major, 'NK')] <- 'NK'
HNSC_IMCISION$celltype_main[HNSC_IMCISION$celltype_r2 %in% bcells] <- 'B'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_r2, 'Plasma')] <- 'Plasma'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_major, 'pDC')] <- 'pDC'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_r2, 'cDC')] <- 'cDC'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_r2, 'Mono')] <- 'Mono'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_r2, 'Macro')] <- 'Macro'
HNSC_IMCISION$celltype_main[str_detect(HNSC_IMCISION$celltype_major, 'Mast')] <- 'Mast'
table(HNSC_IMCISION$celltype_main, HNSC_IMCISION$celltype_major, useNA = 'ifany')
HNSC_IMCISION <- subset(HNSC_IMCISION, subset = celltype_r2 == 'CD56highCD16high', invert = T)
qsave(HNSC_IMCISION, '/bigdata/zlin/Melanoma_meta/data/HNSC_IMCISION/seu_r2.qs')
HNSC_IMCISION <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_IMCISION/seu_r2.qs')

dataset <- 'HNSC_Luoma'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-B','Plasma','GCB')
HNSC_Luoma$celltype_major <- as.character(HNSC_Luoma$celltype_major)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
HNSC_Luoma$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono'))
HNSC_Luoma$celltype_r2 <- pred_concat$celltype_r2[match(colnames(HNSC_Luoma), pred_concat$cell)]
# HNSC_Luoma$celltype_r2[HNSC_Luoma$celltype_major == 'Plasma'] <- 'Plasma'
# HNSC_Luoma$celltype_r2[HNSC_Luoma$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- HNSC_Luoma |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |> 
#   RunPCA(verbose=FALSE) |> 
#   RunUMAP(dims = 1:20) |> 
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 3)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# HNSC_Luoma$celltype_r2[colnames(HNSC_Luoma) %in% colnames(seu)[seu$seurat_clusters == 9]] <- 'Plasmablast'
seu <- HNSC_Luoma |>
  subset(subset = celltype_major %in% c("Plasma","B-cells")) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData()
seu <- seu[,!is.na(seu$celltype_r2)]
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
seu <- seu[,!is.na(seu$celltype_r2)]
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
HNSC_Luoma$celltype_r2[HNSC_Luoma$celltype_major == 'pDC'] <- 'pDC'
HNSC_Luoma$celltype_r2[HNSC_Luoma$celltype_major == 'Mast'] <- 'Mast'
HNSC_Luoma$celltype_major[HNSC_Luoma$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
HNSC_Luoma$celltype_major[HNSC_Luoma$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(HNSC_Luoma$celltype_r2, HNSC_Luoma$celltype_major, useNA = 'ifany')
HNSC_Luoma <- HNSC_Luoma[,!is.na(HNSC_Luoma$celltype_r2)]
HNSC_Luoma$celltype_main <- 'celltype'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_r2, 'CD4')] <- 'CD4+T'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_r2, 'CD8')] <- 'CD8+T'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_r2, 'gdT')] <- 'CD8+T'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_major, 'NK')] <- 'NK'
HNSC_Luoma$celltype_main[HNSC_Luoma$celltype_r2 %in% bcells] <- 'B'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_r2, 'Plasma')] <- 'Plasma'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_major, 'pDC')] <- 'pDC'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_r2, 'cDC')] <- 'cDC'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_r2, 'Mono')] <- 'Mono'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_r2, 'Macro')] <- 'Macro'
HNSC_Luoma$celltype_main[str_detect(HNSC_Luoma$celltype_major, 'Mast')] <- 'Mast'
table(HNSC_Luoma$celltype_main, HNSC_Luoma$celltype_major, useNA = 'ifany')
HNSC_Luoma <- subset(HNSC_Luoma, subset = celltype_r2 == 'CD56highCD16high', invert = T)
qsave(HNSC_Luoma, '/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/seu_r2.qs')
HNSC_Luoma <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/seu_r2.qs')

dataset <- 'NSCLC_Liu'
cell_subtype <- c('CD4','CD8','NK_main')
NSCLC_Liu$celltype_major <- as.character(NSCLC_Liu$celltype_major)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
NSCLC_Liu$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
# pred_concat <- rbind(pred_concat, pred_gdT)
mapping <- setNames(pred_concat$celltype_r2, pred_concat$cell)
NSCLC_Liu$celltype_r2 <- mapping[colnames(NSCLC_Liu)]
NSCLC_Liu$celltype_major[NSCLC_Liu$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels))] <- 'T_cells'
NSCLC_Liu$celltype_major[NSCLC_Liu$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(NSCLC_Liu$celltype_r2, NSCLC_Liu$celltype_major, useNA = 'ifany')
NSCLC_Liu <- NSCLC_Liu[,!is.na(NSCLC_Liu$celltype_r2)]
NSCLC_Liu <- subset(NSCLC_Liu, subset = celltype_major == 'NK_cells', invert = T)
NSCLC_Liu$celltype_main <- 'celltype'
NSCLC_Liu$celltype_main[str_detect(NSCLC_Liu$celltype_r2, 'CD4')] <- 'CD4+T'
NSCLC_Liu$celltype_main[str_detect(NSCLC_Liu$celltype_r2, 'CD8')] <- 'CD8+T'
table(NSCLC_Liu$celltype_main, NSCLC_Liu$celltype_major, useNA = 'ifany')
qsave(NSCLC_Liu, '/bigdata/zlin/Melanoma_meta/data/NSCLC_Liu/seu_r2.qs')

dataset <- 'CRC_Li'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','Endo','CAF','pan-B','Plasma','GCB')
CRC_Li$celltype_major <- as.character(CRC_Li$celltype_major)
CRC_Li <- subset(CRC_Li, subset = celltype_major %in% c('Epithelial cells'), invert = T)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
CRC_Li$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono'))
CRC_Li$celltype_r2 <- pred_concat$celltype_r2[match(colnames(CRC_Li), pred_concat$cell)]
# CRC_Li$celltype_r2[CRC_Li$celltype_major == 'Plasma'] <- 'Plasma'
# CRC_Li$celltype_r2[CRC_Li$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- CRC_Li |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |> 
#   RunPCA(verbose=FALSE) |> 
#   RunUMAP(dims = 1:20) |> 
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 3)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# CRC_Li$celltype_r2[colnames(CRC_Li) %in% colnames(seu)[seu$seurat_clusters == 12]] <- 'Plasmablast'
seu <- CRC_Li |>
  subset(subset = celltype_major %in% c("Plasma","B-cells")) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData()
seu <- seu[,!is.na(seu$celltype_r2)]
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
seu <- seu[,!is.na(seu$celltype_r2)]
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
CRC_Li$celltype_r2[CRC_Li$celltype_major == 'pDC'] <- 'pDC'
CRC_Li$celltype_r2[CRC_Li$celltype_major == 'Mast'] <- 'Mast'
CRC_Li$celltype_major[CRC_Li$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
CRC_Li$celltype_major[CRC_Li$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(CRC_Li$celltype_r2, CRC_Li$celltype_major, useNA = 'ifany')
CRC_Li <- CRC_Li[,!is.na(CRC_Li$celltype_r2)]
CRC_Li$celltype_main <- 'celltype'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_r2, 'CD4')] <- 'CD4+T'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_r2, 'CD8')] <- 'CD8+T'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_r2, 'gdT')] <- 'CD8+T'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_major, 'NK')] <- 'NK'
CRC_Li$celltype_main[CRC_Li$celltype_r2 %in% bcells] <- 'B'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_r2, 'Plasma')] <- 'Plasma'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_major, 'pDC')] <- 'pDC'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_r2, 'cDC')] <- 'cDC'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_r2, 'Mono')] <- 'Mono'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_r2, 'Macro')] <- 'Macro'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_r2, 'EC_')] <- 'Endo'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_r2, 'CAF') | str_detect(CRC_Li$celltype_r2, 'Myofibroblast')] <- 'CAF'
CRC_Li$celltype_main[str_detect(CRC_Li$celltype_major, 'Mast')] <- 'Mast'
table(CRC_Li$celltype_main, CRC_Li$celltype_major, useNA = 'ifany')
CRC_Li <- subset(CRC_Li, subset = celltype_r2 =='CD56highCD16high', invert = T)
qsave(CRC_Li, '/bigdata/zlin/Melanoma_meta/data/CRC_Li/seu_r2.qs')
CRC_Li <- qread('/bigdata/zlin/Melanoma_meta/data/CRC_Li/seu_r2.qs')

dataset <- 'PCa_Hawley'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','Endo','CAF','pan-B','Plasma')
PCa_Hawley$celltype_major <- as.character(PCa_Hawley$celltype_major)
PCa_Hawley <- subset(PCa_Hawley, subset = celltype_major %in% c('Epithelial cells'), invert = T)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
PCa_Hawley$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
# barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
# pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono'))
PCa_Hawley$celltype_r2 <- pred_concat$celltype_r2[match(colnames(PCa_Hawley), pred_concat$cell)]
# PCa_Hawley$celltype_r2[PCa_Hawley$celltype_major == 'Plasma'] <- 'Plasma'
# PCa_Hawley$celltype_r2[PCa_Hawley$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- PCa_Hawley |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |> 
#   RunPCA(verbose=FALSE) |> 
#   RunUMAP(dims = 1:20) |> 
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 1)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# PCa_Hawley$celltype_r2[colnames(PCa_Hawley) %in% colnames(seu)[seu$seurat_clusters == 12]] <- 'Plasmablast'
seu <- PCa_Hawley |>
  subset(subset = celltype_major %in% c("Plasma","B-cells")) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData()
seu <- seu[,!is.na(seu$celltype_r2)]
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
seu <- seu[,!is.na(seu$celltype_r2)]
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
PCa_Hawley$celltype_r2[PCa_Hawley$celltype_major == 'pDC'] <- 'pDC'
PCa_Hawley$celltype_r2[PCa_Hawley$celltype_major == 'Mast'] <- 'Mast'
PCa_Hawley$celltype_major[PCa_Hawley$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
PCa_Hawley$celltype_major[PCa_Hawley$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(PCa_Hawley$celltype_r2, PCa_Hawley$celltype_major, useNA = 'ifany')
PCa_Hawley <- PCa_Hawley[,!is.na(PCa_Hawley$celltype_r2)]
PCa_Hawley$celltype_main <- 'celltype'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'CD4')] <- 'CD4+T'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'CD8')] <- 'CD8+T'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'gdT')] <- 'CD8+T'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_major, 'NK')] <- 'NK'
PCa_Hawley$celltype_main[PCa_Hawley$celltype_r2 %in% bcells] <- 'B'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'Plasma')] <- 'Plasma'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_major, 'pDC')] <- 'pDC'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'cDC')] <- 'cDC'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'Mono')] <- 'Mono'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'Macro')] <- 'Macro'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'EC_')] <- 'Endo'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_r2, 'CAF') | str_detect(PCa_Hawley$celltype_r2, 'Myofibroblast')] <- 'CAF'
PCa_Hawley$celltype_main[str_detect(PCa_Hawley$celltype_major, 'Mast')] <- 'Mast'
table(PCa_Hawley$celltype_main, PCa_Hawley$celltype_major, useNA = 'ifany')
PCa_Hawley <- subset(PCa_Hawley, subset = celltype_r2 =='CD56highCD16high', invert = T)
qsave(PCa_Hawley, '/bigdata/zlin/Melanoma_meta/data/PCa_Hawley/seu_r2.qs')
PCa_Hawley <- qread('/bigdata/zlin/Melanoma_meta/data/PCa_Hawley/seu_r2.qs')

dataset <- 'TNBC_Shiao'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-B','Plasma','GCB')
TNBC_Shiao$celltype_major <- as.character(TNBC_Shiao$celltype_major)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
TNBC_Shiao$celltype_r2 <- 'Unresolved'
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2', 'Macro', 'Mono'))
TNBC_Shiao$celltype_r2 <- pred_concat$celltype_r2[match(colnames(TNBC_Shiao), pred_concat$cell)]
# TNBC_Shiao$celltype_r2[TNBC_Shiao$celltype_major == 'Plasma'] <- 'Plasma'
# TNBC_Shiao$celltype_r2[TNBC_Shiao$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- TNBC_Shiao |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |> 
#   RunPCA(verbose=FALSE) |> 
#   RunUMAP(dims = 1:20) |> 
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 3)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# TNBC_Shiao$celltype_r2[colnames(TNBC_Shiao) %in% colnames(seu)[seu$seurat_clusters == 44]] <- 'Plasmablast'
seu <- TNBC_Shiao |>
  subset(subset = celltype_major %in% c("Plasma","B-cells")) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData()
seu <- seu[,!is.na(seu$celltype_r2)]
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
seu <- seu[,!is.na(seu$celltype_r2)]
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
TNBC_Shiao$celltype_r2[TNBC_Shiao$celltype_major == 'pDC'] <- 'pDC'
TNBC_Shiao$celltype_r2[TNBC_Shiao$celltype_major == 'Mast'] <- 'Mast'
TNBC_Shiao$celltype_major[TNBC_Shiao$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
TNBC_Shiao$celltype_major[TNBC_Shiao$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(TNBC_Shiao$celltype_r2, TNBC_Shiao$celltype_major, useNA = 'ifany')
TNBC_Shiao <- TNBC_Shiao[,!is.na(TNBC_Shiao$celltype_r2)]
TNBC_Shiao$celltype_main <- 'celltype'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_r2, 'CD4')] <- 'CD4+T'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_r2, 'CD8')] <- 'CD8+T'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_r2, 'gdT')] <- 'CD8+T'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_major, 'NK')] <- 'NK'
TNBC_Shiao$celltype_main[TNBC_Shiao$celltype_r2 %in% bcells] <- 'B'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_r2, 'Plasma')] <- 'Plasma'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_major, 'pDC')] <- 'pDC'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_r2, 'cDC')] <- 'cDC'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_r2, 'Mono')] <- 'Mono'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_r2, 'Macro')] <- 'Macro'
TNBC_Shiao$celltype_main[str_detect(TNBC_Shiao$celltype_major, 'Mast')] <- 'Mast'
table(TNBC_Shiao$celltype_main, TNBC_Shiao$celltype_major, useNA = 'ifany')
TNBC_Shiao <- subset(TNBC_Shiao, subset = (celltype_r2 == 'CD56highCD16high') | (time_point == 'Post2'), invert = T)
qsave(TNBC_Shiao, '/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/seu_r2.qs')
TNBC_Shiao <- qread('/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/seu_r2.qs')

dataset <- 'HNSC_Franken'
cell_subtype <- c('CD4','CD8','NK_main','Myeloids_major','cDC2','Mono','Macro','Endo','CAF','pan-B','Plasma','GCB')
HNSC_Franken$celltype_major <- as.character(HNSC_Franken$celltype_major)
HNSC_Franken <- subset(HNSC_Franken, subset = celltype_major %in% c('Melanoma', 'Epithelial cells'), invert = T)
pred_list <- lapply(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/', cell_subtype, ".qs"), qread)
names(pred_list) <- cell_subtype
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/barcode_gdT.csv')) 
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2','Macro', 'Mono'))
HNSC_Franken$celltype_r2 <- pred_concat$celltype_r2[match(colnames(HNSC_Franken), pred_concat$cell)]
# HNSC_Franken$celltype_r2[HNSC_Franken$celltype_major == 'Plasma'] <- 'Plasma'
# HNSC_Franken$celltype_r2[HNSC_Franken$celltype_r2 == 'Plasmablasts'] <- 'Plasma'
# seu <- HNSC_Franken |>
#   subset(subset = celltype_r2 %in% c("Plasma")) |> 
#   NormalizeData() |>
#   FindVariableFeatures() |>
#   ScaleData() |> 
#   RunPCA(verbose=FALSE) |> 
#   RunUMAP(dims = 1:20) |> 
#   FindNeighbors(dims = 1:20) |> FindClusters(resolution = 3)
# DotPlot(seu, features = c('MS4A1','CD79A','BTG1','HLA-A',
#                           'IGHM','IGHD','CD86','CD27','CD40','BCL6','AICDA',
#                           'SUGCT','RGS13','NEIL1','CDCA7','POU2AF1','IGHG1','IGHG2','IGHA1','IGHA2',
#                           'MKI67','TOP2A','JCHAIN','MZB1','XBP1','SDC1','CD38')) + RotatedAxis()
# HNSC_Franken$celltype_r2[colnames(HNSC_Franken) %in% colnames(seu)[seu$seurat_clusters == 12]] <- 'Plasmablast'
seu <- HNSC_Franken |>
  subset(subset = celltype_major %in% c("Plasma","B-cells")) |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData()
seu <- seu[,!is.na(seu$celltype_r2)]
genes_to_check = c('TCL1A','YBX3','FCER2',
                   'CD27','CD19','IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                   'HSPA1A','DNAJB1','HSPA1B',
                   'MT1X','MT2A','SLC30A1',
                   'TNF','EGR1',
                   'ZNF331','VPS37B','NR4A2','CD69','CD86','IGHM','IGHD',
                   'CD83','CCR7','NFKBID',
                   'CRIP1','S100A10','S100A4', 'ITGB1',
                   'DUSP4','FCRL5','ZEB2','ITGAX',
                   'ENO1','PSME2','NME1',
                   'SUGCT','NEIL1','MEF2B',
                   'LMO2','GMDS','MARCKSL1',
                   'STMN1','HMGB2','MKI67','TOP2A',
                   'TXNDC5','MYDGF',
                   'XBP1','MZB1','JCHAIN','MS4A1')
seu <- seu[,!is.na(seu$celltype_r2)]
DotPlot(seu, group.by = 'celltype_r2', features = genes_to_check) + RotatedAxis() + scale_colour_gradient2(low = "#FFFFD9", mid = "#41B6C4", high = "#081D58")
bcells <- unique(seu$celltype_r2)
HNSC_Franken$celltype_r2[HNSC_Franken$celltype_major == 'pDC'] <- 'pDC'
HNSC_Franken$celltype_r2[HNSC_Franken$celltype_major == 'Mast'] <- 'Mast'
HNSC_Franken$celltype_major[HNSC_Franken$celltype_r2 %in% c(unique(pred_list[[1]]$pruned.labels), unique(pred_list[[2]]$pruned.labels), 'gdT')] <- 'T_cells'
HNSC_Franken$celltype_major[HNSC_Franken$celltype_r2 %in% c(unique(pred_list[[3]]$pruned.labels))] <- 'NK_cells'
table(HNSC_Franken$celltype_r2, HNSC_Franken$celltype_major, useNA = 'ifany')
HNSC_Franken <- HNSC_Franken[,!is.na(HNSC_Franken$celltype_r2)]
HNSC_Franken$celltype_main <- 'celltype'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_r2, 'CD4')] <- 'CD4+T'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_r2, 'CD8')] <- 'CD8+T'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_r2, 'gdT')] <- 'CD8+T'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_major, 'NK')] <- 'NK'
HNSC_Franken$celltype_main[HNSC_Franken$celltype_r2 %in% bcells] <- 'B'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_r2, 'Plasma')] <- 'Plasma'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_major, 'pDC')] <- 'pDC'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_major, 'Mast')] <- 'Mast'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_r2, 'cDC')] <- 'cDC'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_r2, 'Mono')] <- 'Mono'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_r2, 'Macro')] <- 'Macro'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_r2, 'EC_')] <- 'Endo'
HNSC_Franken$celltype_main[str_detect(HNSC_Franken$celltype_r2, 'CAF') | str_detect(HNSC_Franken$celltype_r2, 'Myofibroblast')] <- 'CAF'
table(HNSC_Franken$celltype_main, HNSC_Franken$celltype_major, useNA = 'ifany')
HNSC_Franken <- subset(HNSC_Franken, subset = celltype_r2 == 'CD56highCD16high', invert = T)
qsave(HNSC_Franken, '/bigdata/zlin/Melanoma_meta/data/HNSC_Franken/seu_r2.qs')

# 
datasets <- c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'BCC_Yost', 'SCC_Yost', 'HNSC_IMCISION', 'HNSC_Luoma', 'NSCLC_Liu', 'CRC_Li', 'PCa_Hawley', 'TNBC_Shiao', 'HNSC_Franken')
lapply(datasets, function(dataset){
  print(dataset)
  seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) 
  # Modification
  # CD4+T
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Tn_IL7R-'] <- 'CD4_Prolif'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Tn_ADSL'] <- 'CD4_Naive'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_pre-Tfh_CXCR5+'] <- 'CD4_pre-Tfh_CXCR5'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Tm_CAPG+CREM-'] <- 'CD4_Tm_CAPG'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Tm_TNF'] <- 'CD4_Tm_CREM-'
  seu$celltype_r2[seu$celltype_r2 %in% c('CD4_Treg_TNFRSF9-', 'CD4_Treg_S1PR1')] <- 'CD4_Treg_Early'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Th_ISG'] <- 'CD4_Th_ISG15'
  seu$celltype_r2[seu$celltype_r2 == 'CD4_Treg_ISG'] <- 'CD4_Treg_ISG15'
  # CD8
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_IL7R'] <- 'CD8_Tcm_IL7R'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tex_TCF7'] <- 'CD8_Tpex_TCF7'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_NME1'] <- 'CD8_Prolif'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_ZNF683'] <- 'CD8_Trm_ZNF683'
  seu$celltype_r2[seu$celltype_r2 %in% c('CD8_NK-like_EOMES', 'CD8_NK-like_TXK')] <- 'CD8_NK-like'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tm_ZNF683'] <- 'CD8_Trm_ZNF683'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_MAIT_SLC4A10'] <- 'MAIT'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_Tex_OXPHOS-'] <- 'CD8_Tex_CXCL13'
  seu$celltype_r2[seu$celltype_r2 == 'CD8_ISG'] <- 'CD8_Tex_ISG15'
  # NK
  seu$celltype_r2[seu$celltype_r2 == 'CD56highCD16low'] <- 'NK_CD56hiCD16lo'
  seu$celltype_r2[seu$celltype_r2 == 'CD56lowCD16high'] <- 'NK_CD56loCD16hi'
  # DC
  seu$celltype_r2[seu$celltype_r2 == 'cDC2_CD1A'] <- 'DC_LC-like'
  seu$celltype_r2[seu$celltype_r2 == 'cDC2_FCN1'] <- 'MoDC'
  seu$celltype_r2[seu$celltype_r2 == 'cDC2_CXCR4hi'] <- 'cDC2_CD1C'
  seu$celltype_r2[seu$celltype_r2 == 'cDC3'] <- 'MigrDC'
  # Macro
  seu$celltype_r2[seu$celltype_r2 == 'Macro_GPNMB'] <- 'Macro_TREM2'
  seu$celltype_r2[seu$celltype_r2 == 'Macro_IL1B'] <- 'Macro_TNF'
  seu$celltype_r2[seu$celltype_r2 == 'Macro_NLRP3'] <- 'Macro_IL1B'
  # # B/Plasma
  # seu$celltype_r2 <- str_replace(seu$celltype_r2,' cells','')
  # seu$celltype_r2[seu$celltype_r2 == 'Exhausted B'] <- 'GC-like B'
  # seu$celltype_r2[seu$celltype_r2 == 'Non-switched memory B'] <- 'Memory IgM+ B'
  # seu$celltype_r2[seu$celltype_r2 == 'Switched memory B'] <- 'Memory IgM- B'
  # B/Plasma
  seu$celltype_r2[seu$celltype_r2 == 'B.01.TCL1A+naiveB'] <- 'B_Naive'
  seu$celltype_r2[seu$celltype_r2 == 'B.02.IFIT3+B'] <- 'B_ISG15'
  seu$celltype_r2[seu$celltype_r2 == 'B.03.HSP+B'] <- 'B_HSP'
  seu$celltype_r2[seu$celltype_r2 == 'B.04.MT1X+B'] <- 'B_MT2A'
  seu$celltype_r2[seu$celltype_r2 == 'B.05.EGR1+ACB'] <- 'ACB_EGR1'
  seu$celltype_r2[seu$celltype_r2 == 'B.06.NR4A2+ACB2'] <- 'ACB_NR4A2'
  seu$celltype_r2[seu$celltype_r2 == 'B.07.CCR7+ACB3'] <- 'ACB_CCR7'
  seu$celltype_r2[seu$celltype_r2 == 'B.08.ITGB1+SwBm'] <- 'B_Memory'
  seu$celltype_r2[seu$celltype_r2 == 'B.09.DUSP4+AtM'] <- 'B_AtM'
  seu$celltype_r2[seu$celltype_r2 == 'B.10.ENO1+Pre_GCB'] <- 'GCB_Pre'
  seu$celltype_r2[seu$celltype_r2 == 'B.11.SUGCT+DZ_GCB'] <- 'GCB_SUGCT'
  seu$celltype_r2[seu$celltype_r2 == 'B.12.LMO2+LZ_GCB'] <- 'GCB_LMO2'
  seu$celltype_r2[seu$celltype_r2 == 'B.13.Cycling_GCB'] <- 'GCB_Prolif'
  seu$celltype_r2[seu$celltype_r2 == 'B.14.Plasmablast'] <- 'Plasmablast'
  seu$celltype_r2[seu$celltype_r2 == 'B.15.Plasma cell'] <- 'Plasma_cell'
  # non-immune
  seu$celltype_r2[str_detect(seu$celltype_r2, 'EndMT')] <- 'EndMT'
  
  qsave(seu, paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) 
  if (!dir.exists(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2/'))){
    dir.create(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2/'))
  }
  Export10X(seu, dir =paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2/'), 
            append_reductions = NULL, gzip = F)
})





