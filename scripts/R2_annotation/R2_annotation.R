rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','SingleR','ggsci','tibble','qs','qs2','BiocParallel','scGate','Matrix','SingleCellExperiment','scran','parallel','scGate','janitor','RColorBrewer')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
options(max.print = 10000)
# setwd("/home/zlin/workspace/PanCancer_ICI")
# Loading dataset
SKCM_Becker <- qs_read('data/SKCM_Becker/seu_r1.qs2')
SKCM_Becker$cohort <- 'SKCM_Becker'
SKCM_Plozniak <- qs_read('data/SKCM_Plozniak/seu_r1.qs2')
BRCA_Bassez1 <- qs_read('data/BRCA_Bassez1/seu_r1.qs2')
BRCA_Bassez2 <- qs_read('data/BRCA_Bassez2/seu_r1.qs2')
TNBC_Zhang <- qs_read('data/TNBC_Zhang/seu_r1.qs2')
TNBC_Shiao <- qs_read('data/TNBC_Shiao/seu_r1.qs2')
HNSC_Franken <- qs_read('data/HNSC_Franken/seu_r1.qs2')
HNSC_vanderLeun <- qs_read('data/HNSC_vanderLeun/seu_r1.qs2')
HNSC_vanderLeun$cohort <- 'HNSC_vanderLeun'
HNSC_Luoma <- qs_read('data/HNSC_Luoma/seu_r1.qs2')
CRC_Li <- qs_read('data/CRC_Li/seu_r1.qs2')
CRC_Chen <- qs_read('data/CRC_Chen/seu_r1.qs2')
BCC_Yost <- qs_read('data/BCC_Yost/seu_r1.qs2')
BCC_Yost1 <- BCC_Yost |> 
  subset(subset = patient %in% c("BCC_Yost_su009", "BCC_Yost_su012"), invert = T)
BCC_Yost2 <- BCC_Yost |> 
  subset(subset = patient %in% c("BCC_Yost_su009", "BCC_Yost_su012"))
SCC_Yost <- qs_read('data/SCC_Yost/seu_r1.qs2')
NSCLC_Liu <- qs_read('data/NSCLC_Liu/seu_r1.qs2')
PCa_Hawley <- qs_read('data/PCa_Hawley/seu_r1.qs2')

# Run SingleR with mutiple reference datasets
# arguments:
# seu: seurat object 
# ref_list: reference (list)
# major: major cell types (T/NK, Myeloids)
# type_name: for file name of the predicjtion results
# datset: for file name of the prediction results
SingleRMultiRef <- function(seu, ref_list, if.subset = F, major, type_name, n_cores = 10){
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
  qs_save(pred, paste0("data/", unique(seu$cohort), "/", type_name, ".qs2"))
  return(pred)
}
# run SingleR with one reference dataset
SingleROneRef <- function(seu, ref, label, major, type_name, n_cores = 10){
  pred <- seu |> 
    subset(subset = celltype_major %in% major) |> 
    as.SingleCellExperiment() |> 
    logNormCounts() |> 
    SingleR(ref = ref, labels = label, assay.type.test=1, BPPARAM=MulticoreParam(n_cores))
  qs_save(pred, paste0("data/", unique(seu$cohort), "/", type_name, ".qs2"))
}

# T_NK
ref_T <- qs_read("data/Ref_SingleR/T_ref.qs2")
ref_CD4 <- qs_read("data/Ref_SingleR/T_CD4_ref.qs2")
ref_CD8 <- qs_read("data/Ref_SingleR/T_CD8_ref.qs2")
ref_main_NK <- qs_read("data/Ref_SingleR/NK_main_ref.qs2")
ref_fine_NK <- qs_read("data/Ref_SingleR/NK_fine_ref.qs2")

hierSingleRMultiRef_T_NK <- function(seu, main_ref_list, ref_cd4, ref_list_cd8, ref_list_nk_main, ref_list_nk_fine, major = c("CD8+ T-cells","CD4+ T-cells","Cycling T/NK","NK cells")){
  print(unique(seu$cohort))
  seu <- subset(seu, subset = celltype_major %in% major) 
  seu <- NormalizeData(seu) 
  print('Gating NK cells')
  NK_KLRD1 <- gating_model(name = "NK_KLRD1", signature = c("KLRD1","CD3D-"))
  NK_NCAM1 <- gating_model(name = "NK_NCAM1", signature = c("NCAM1","CD3D-"))
  seu_1<- scGate(data = seu, model = NK_KLRD1)
  seu_2<- scGate(data = seu, model = NK_NCAM1)
  barcode_NK <- unique(c(colnames(seu_1)[seu_1$is.pure == 'Pure'],
                         colnames(seu_2)[seu_2$is.pure == 'Pure']))
  print(paste0(length(barcode_NK), ' NK cells were detected.'))
  write.csv(barcode_NK, paste0('data/', unique(seu$cohort), '/barcode_NK.csv'))
  print('Gating gdT cells')
  gdT <- gating_model(name = "gdT", signature = c("CD3D","TRDV1","TRGC2","TRDV2","TRGV5","TRGV9","TRGV10","TRDC","CD8A-","CD4-"))
  seu_gdT<- scGate(data = seu, model = gdT)
  barcode_gdT <- colnames(seu_gdT)[seu_gdT$is.pure == 'Pure']
  # barcode_gdT <- setdiff(barcode_gdT, barcode_NK)
  print(paste0(length(barcode_gdT), ' gdT cells were detected'))
  write.csv(barcode_gdT, paste0('data/', unique(seu$cohort), '/barcode_gdT.csv'))
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
  pred_CD4 <- SingleROneRef(seu = seu[, barcode_CD4], ref = ref_cd4, label = ref_cd4$label, type_name = 'CD4_nm_2023', major = major)
  print('CD8')
  pred_CD8 <- SingleROneRef(seu = seu[, barcode_CD8], ref = ref_list_cd8[[1]], label = ref_list_cd8[[1]]$label, type_name = 'CD8_science_2021', major = major)
  pred_CD8 <- SingleROneRef(seu = seu[, barcode_CD8], ref = ref_list_cd8[[2]], label = ref_list_cd8[[2]]$label, type_name = 'CD8_nm_2023', major = major)
  print('NK main')
  pred_NK_main <- SingleRMultiRef(seu = seu[, barcode_NK], ref_list = ref_list_nk_main, type_name = 'NK_main')
  print('NK fine')
  pred_NK_fine <- SingleRMultiRef(seu = seu[, barcode_NK], ref_list = ref_list_nk_fine, type_name = 'NK_fine')
}

mclapply(list(SKCM_Becker, SKCM_Plozniak, BCC_Yost1,
              BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, TNBC_Shiao,
              HNSC_Franken, HNSC_vanderLeun, HNSC_Luoma, 
              CRC_Li, CRC_Chen, PCa_Hawley), function(seu){
                hierSingleRMultiRef_T_NK(seu, main_ref_list = ref_T, 
                                         ref_cd4 = ref_CD4, ref_list_cd8 = ref_CD8, 
                                         ref_list_nk_main = ref_main_NK, ref_list_nk_fine = ref_fine_NK,
                                         major = c("CD8+ T-cells","CD4+ T-cells","Cycling T/NK","NK cells"))}, mc.cores = 15)

hierSingleRMultiRef_T_NK(BCC_Yost1, main_ref_list = ref_T,
                         ref_cd4 = ref_CD4, ref_list_cd8 = ref_CD8,
                         ref_list_nk_main = ref_main_NK, ref_list_nk_fine = ref_fine_NK,
                         major = c("CD8+ T-cells","CD4+ T-cells","Cycling T/NK","NK cells"))

# T cells only (CD45+CD3+ sorted)
hierSingleRMultiRef_T <- function(seu, main_ref_list, ref_cd4, ref_list_cd8, major){
  print(unique(seu$cohort))
  seu <- NormalizeData(seu) 
  print('Gating gdT cells')
  gdT <- gating_model(name = "gdT", signature = c("CD3D","TRDV1","TRGC2","TRDV2","TRGV5","TRGV9","TRGV10","TRDC","CD8A-","CD4-"))
  seu_gdT<- scGate(data = seu, model = gdT)
  barcode_gdT <- colnames(seu_gdT)[seu_gdT$is.pure == 'Pure']
  # barcode_gdT <- setdiff(barcode_gdT, barcode_NK)
  print(paste0(length(barcode_gdT), ' gdT cells were detected'))
  write.csv(barcode_gdT, paste0('data/', unique(seu$cohort), '/barcode_gdT_sorted.csv'))
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
  pred_T <- SingleRMultiRef(seu[, !colnames(seu) %in% c(barcode_gdT, barcode_cd4, barcode_cd8)], ref_list = main_ref_list, type_name = 'sorted_T') # For CD45CD3 sorted (SCC_Yost)
  barcode_CD4 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD4'], barcode_cd4)
  barcode_CD8 <- c(rownames(pred_T)[pred_T$pruned.labels == 'CD8'], barcode_cd8)
  print('CD4')
  pred_CD4 <- SingleROneRef(seu = seu[, barcode_CD4], , ref = ref_cd4, label = ref_cd4$label, type_name = 'sorted_CD4_nm_2023', major = major)
  print('CD8')
  pred_CD8 <- SingleROneRef(seu = seu[, barcode_CD8], ref = ref_list_cd8[[1]], label = ref_list_cd8[[1]]$label, type_name = 'sorted_CD8_science_2021', major = major)
  pred_CD8 <- SingleROneRef(seu = seu[, barcode_CD8], ref = ref_list_cd8[[2]], label = ref_list_cd8[[2]]$label, type_name = 'sorted_CD8_nm_2023', major = major)
}
hierSingleRMultiRef_T(BCC_Yost2, main_ref_list = ref_T, ref_cd4 = ref_CD4, ref_list_cd8 = ref_CD8, major = c("CD8+ T-cells","CD4+ T-cells","Cycling T/NK"))
hierSingleRMultiRef_T(SCC_Yost, main_ref_list = ref_T, ref_cd4 = ref_CD4, ref_list_cd8 = ref_CD8, major = c("T-cells"))

# Myeloids (doi.org/10.1016/j.cell.2021.01.010)
ref_Myeloids <- qread("data/Ref_SingleR/Myeloids_major_ref.qs")
ref_cdc2 <- qread("data/Ref_SingleR/Myeloids_cdc2_ref.qs")
ref_mono <- qread("data/Ref_SingleR/Myeloids_mono_ref.qs")
ref_macro <- qread("data/Ref_SingleR/Myeloids_macro_ref.qs")
cohorts <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, TNBC_Shiao, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, CRC_Li, PCa_Hawley)

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
mclapply(list(SKCM_Becker, SKCM_Plozniak, BCC_Yost1,
              BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, TNBC_Shiao,
              HNSC_Franken, HNSC_vanderLeun, HNSC_Luoma, 
              CRC_Li, CRC_Chen, PCa_Hawley), function(seu){
                hierSingleRMultiRef_Myeloids(seu, main_ref_list = ref_Myeloids, 
                                             ref_list_cdc2 = ref_cdc2, 
                                             ref_list_mono = ref_mono, 
                                             ref_list_macro = ref_macro)}, mc.cores = 5)
hierSingleRMultiRef_Myeloids(HNSC_vanderLeun, main_ref_list = ref_Myeloids, 
                             ref_list_mono = ref_mono, ref_list_macro = ref_macro, ref_list_cdc2 = ref_cdc2)

# # B cells 
# ref_B <- MonacoImmuneData()
# bped <- celldex::BlueprintEncodeData()
# cohorts <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, CRC_Li, PCa_Hawley)
# ref_B <- ref_B[, str_detect(ref_B$label.fine, 'B cells') | str_detect(ref_B$label.fine, 'Plasmablasts')]
# lapply(cohorts, function(seu){
#   print(unique(seu$cohort))
#   SingleROneRef(seu, ref = ref_B, label = ref_B$label.fine, major = 'B-cells', type_name = 'B')})
# SingleROneRef(HNSC_Franken, ref = ref_B, label = ref_B$label.fine, major = 'B-cells', type_name = 'B')

# pan-cancer B (doi.org/10.1126/science.adj4857)
ref_Bplasma <- qread("data/Ref_SingleR/Pan_B_major.qs")
ref_b <- qread("data/Ref_SingleR/B.qs")
ref_gcb <- qread("data/Ref_SingleR/GCB.qs")
ref_plasma <- qread("data/Ref_SingleR/Plasma.qs")
ref_asc <- qread("data/Ref_SingleR/B_asc.qs")
hierSingleRMultiRef_Bplasma <- function(seu, main_ref_list, ref_list_b, ref_list_gcb, ref_list_plasma){
  pred <- SingleRMultiRef(seu = seu, ref_list = main_ref_list, if.subset = T, major = c('B-cells', 'Plasma cells'), type_name = 'Bplasma_major')
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
  pred_asc <- SingleROneRef(seu = seu[, pred_plasma], ref = ref_list_plasma, label = ref_list_plasma$label, major = c('B-cells', 'Plasma cells'), type_name = 'Plasma')
  # pred_plasma <- SingleRMultiRef(seu = seu[, pred_plasma], ref_list = ref_list_plasma, type_name = 'Plasma')
  }
}

# pred <- seu[, pred_plasma] |> 
#   as.SingleCellExperiment() |> 
#   logNormCounts() |> 
#   SingleR(ref = ref_asc, labels = ref_asc$label, assay.type.test=1, BPPARAM=MulticoreParam(50))

mclapply(list(SKCM_Becker, SKCM_Plozniak, BCC_Yost,
              BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, TNBC_Shiao,
              HNSC_Franken, HNSC_vanderLeun, HNSC_Luoma, 
              CRC_Li, CRC_Chen, PCa_Hawley), function(seu){
                hierSingleRMultiRef_Bplasma(seu, main_ref_list = ref_Bplasma, 
                                            ref_list_b = ref_b, 
                                            ref_list_gcb = ref_gcb, 
                                            ref_list_plasma = ref_asc)}, mc.cores = 10)

hierSingleRMultiRef_Bplasma(HNSC_vanderLeun, main_ref_list = ref_Bplasma, ref_list_b = ref_b, 
                             ref_list_gcb = ref_gcb, ref_list_plasma = ref_asc)

seu_sub_b <- subset(BRCA_Bassez1, subset = celltype_major %in% c("B-cells", "Plasma")) |> as.SingleCellExperiment() |> logNormCounts()
pred <- SingleR(seu_sub_b, ref = sce, labels = sce$Annotation, BPPARAM=MulticoreParam(30))
seu_sub_b$singler <- pred$pruned.labels
seu_sub_b$singler[is.na(seu_sub_b$singler)] <- 'unknown'
seu_sub_b <- seu_sub_b |> as.Seurat()
genes_to_check <- c('FCER2', 'TCL1A', 'IL4R', 'CD72', 'BACH2', 'IGHD', 'IGHM',
                   'NR4A1', 'NR4A2', 'CREM', 'CD83', 
                   'ISG15', 'IFI44L', 'IFI6', 'IFIT3', 
                   'CD27', 'TNFRSF13B', 'TXNIP', 'GPR183', 
                   'HSPA1A', 'HSPA1B', 'DNAJB1', 'EGR1', 
                   'FCRL4', 'FCRL5', 'ITGAX', 'TBX21', 'CR2', 
                   'CCR1', 'CXCR3', 'PDCD1', 'HOK', 'FCRL3', 'FGR', 
                   'NME1', 'APEX1', 'POLD2', 'POLE3', 'MYC', 
                   'BCL6', 'RGS13', 'ACIDA', 'IL21R', 
                   'MKI67', 'STMN1', 'HMGB2', 'TOP2A', 
                   'CD38', 'MZB1', 'PRDM1', 'IRF4', 'XBP1', 
                   'MS4A1', 'LTB', 'HLA-DRA', 'HLA-DRB1', 'HLA-DPA1', 'HLA-DQA1', 
                   'IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHA1', 'IGHA2', 
                   'IL10', 'IL12A', 'EBI3', 'TGFB1', 'IL35')
DotPlot(seu_sub_b, group.by = 'singler', features = genes_to_check) + RotatedAxis()
df_marker <- readxl::read_xlsx('data/GSE233236(ref_B)/1-s2.0-S0092867424007128-mmc3.xlsx')

# # pan-cancer B (doi.org/10.1016/j.cell.2024.06.038)
# ref_BPlasma <- qread("data/Ref_SingleR/Pan_B-major.qs2")
# ref_b <- qread("data/Ref_SingleR/B_naive-memory.qs2")
# ref_cg <- qread("data/Ref_SingleR/B_cycling-gc.qs2")
# ref_asc <- qread("data/Ref_SingleR/B_asc.qs2")
# hierSingleRMultiRef_BPlasma <- function(seu, main_ref_list, ref_1, ref_2, ref_3){
#   pred <- SingleRMultiRef(seu = seu, ref_list = main_ref_list, if.subset = T, major = c('B-cells', 'Plasma cells'), type_name = 'BPlasma_major')
#   print('Major done!')
#   pred_b <- rownames(pred)[pred$pruned.labels == 'naive/memory']
#   pred_cg <- rownames(pred)[pred$pruned.labels == 'cycling/gc']
#   pred_asc <- rownames(pred)[pred$pruned.labels == 'asc']
#   print('B')
#   pred_b <- SingleROneRef(seu = seu[, pred_b], ref = ref_1, label = ref_1$label, major = c('B-cells', 'Plasma cells'), type_name = 'B_naivememory')
#   if (length(pred_cg)>0){
#     print('GCB')
#     pred_cg <- SingleROneRef(seu = seu[, pred_cg], ref = ref_2, label = ref_2$label, major = c('B-cells', 'Plasma cells'), type_name = 'B_cylinggc')
#   }
#   if (length(pred_asc)>0){
#     print('Plasma')
#     pred_asc <- SingleROneRef(seu = seu[, pred_asc], ref = ref_3, label = ref_3$label, major = c('B-cells', 'Plasma cells'), type_name = 'ASC')
#   }
# }
# cohorts <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, TNBC_Shiao, BCC_Yost1, HNSC_IMCISION, HNSC_Luoma, HNSC_Franken, CRC_Li, CRC_Chen, PCa_Hawley)
# mclapply(cohorts, function(seu){
#   hierSingleRMultiRef_BPlasma(seu, main_ref_list = ref_BPlasma, ref_1 = ref_b, ref_2 = ref_cg, ref_3 = ref_asc)}, mc.cores = 100)
# 
# Pan-cancer Endothelial cells (doi.org/10.1093/nsr/nwae231)
ref_endo_major <- qread('data/Ref_SingleR/Endo_major.qs')
ref_endo_vascular <- qread('data/Ref_SingleR/Endo_vascular.qs')
hierSingleRMultiRef_Endo <- function(seu, main_ref_list, ref_list_vas){
  pred <- SingleRMultiRef(seu = seu, ref_list = main_ref_list, if.subset = T, major = c('Endothelial cells'), type_name = 'pan-Endo')
  print('Major done!')
  pred_vas <- rownames(pred)[pred$pruned.labels == 'vascular']
  print('Vascular')
  pred_vas <- SingleRMultiRef(seu = seu[, pred_vas], ref_list = ref_list_vas, type_name = 'pan-Endo_vas')
}
cohorts <- list(SKCM_Becker, SKCM_Plozniak, BRCA_Bassez1, BRCA_Bassez2, TNBC_Shiao, HNSC_Franken, BCC_Yost1, CRC_Li, CRC_Chen, PCa_Hawley)
mclapply(cohorts, function(seu){
  hierSingleRMultiRef_Endo(seu,
                           main_ref_list = ref_endo_major,
                           ref_list_vas = ref_endo_vascular)}, mc.cores = 5)
hierSingleRMultiRef_Endo(TNBC_Shiao, main_ref_list = ref_endo_major, ref_list_vas = ref_endo_vascular)
# # Markers to validate
# genes_to_check = c('GJA4','GJA5','FBLN5',
#                    'CA4','CD36','RGCC',
#                    'PROX1','LYVE1','CCL2',
#                    'COL4A1','KDR','ESM1',
#                    'ACKR1','SELP','CLU')

# # Endothelial cells
# ref_Endo <- qread("data/Ref_SingleR/Endo_ref.qs2")
# cohorts <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, BCC_Yost1, CRC_Li, CRC_Chen, PCa_Hawley)
# lapply(cohorts, function(seu){SingleROneRef(seu, ref = ref_Endo, label = ref_Endo$label, major = 'Endothelial cells', type_name = 'Endo')})
# SingleROneRef(HNSC_Franken, ref = ref_Endo, label = ref_Endo$label, major = 'Endothelial cells', type_name = 'Endo')

# CAF (doi.org/10.1038/s41467-024-48310-4)
ref_CAF <- qread("data/Ref_SingleR/CAF_ref.qs")
cohorts <- list(SKCM_Becker, SKCM_Plozniak, BRCA_Bassez1, BRCA_Bassez2, TNBC_Shiao, HNSC_Franken, BCC_Yost, CRC_Li, CRC_Chen, PCa_Hawley)
lapply(cohorts, function(seu){
  print(unique(seu$cohort))
  SingleRMultiRef(seu, ref_list = ref_CAF, major = c('Fibroblasts','Myocytes','Mural cells','Pericytes'), if.subset = T, type_name = 'CAF')})
SingleRMultiRef(TNBC_Shiao, ref_list = ref_CAF, major = c('Fibroblasts','Myocytes','Mural cells','Pericytes'), if.subset = T, type_name = 'CAF')

genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38','JCHAIN'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('CXCR1', 'CXCR2', 'PTGS2','OLR1', 'VEGFA'),
                      c('COL3A1', 'FAP', 'COL1A1'), 
                      c('ACTA2', "RGS5", "COX4I2","DCN"),
                      c("DES", "TNNT3", "COX6A2", "ACTC1",  "MYL1"),
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('KRT15','KRT17','KRT19','EPCAM'),
                      c('MKI67','TOP2A')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Neu','Fibro','PC','SMC','Endo','Mela','Epi','Proliferating')
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

add_to_seu_r2 <- function(cd45_sorted = F, .cohort, pred_res, celltype_major_remove, label_overlapped, celltype_r2_remove){
  seu <- qs_read(paste0('data/', .cohort,'/seu_r1.qs2'))
  seu$celltype_major <- as.character(seu$celltype_major)
  pred_list <- lapply(paste0('data/', .cohort, '/', pred_res, ".qs2"), qs_read)
  names(pred_list) <- pred_res
  # adding CD8+ stressed T cell
  pred_list$CD8_science_2021$pruned.labels[pred_list$CD8_nm_2023$pruned.labels == 'CD8_c4_Tstr'] <- 'CD8_stress'
  pred_list$CD8_nm_2023 <- NULL
  pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
    data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
  }))
  barcode_gdT <- read.csv(paste0('data/', cohort, '/barcode_gdT.csv')) 
  if (nrow(barcode_gdT) > 0) {
    pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
    pred_concat <- rbind(pred_concat, pred_gdT)
  } else {
    message("barcode_gdT is empty. Skipping processing for gdT.")
  }
  pred_concat <- filter(pred_concat, !celltype_r2 %in% label_overlapped)
  seu$celltype_r2 <- pred_concat$celltype_r2[match(colnames(seu), pred_concat$cell)]
  seu$celltype_r2[seu$celltype_major == 'pDC'] <- 'pDC'
  seu$celltype_r2[seu$celltype_major == 'Mast'] <- 'Mast'
  # seu$celltype_major[seu$celltype_r2 %in% na.omit(c(unique(pred_list$CD4_nm_2023$labels)), na.omit(unique(pred_list$CD8_science_2021$labels)), 'gdT')] <- 'T_cells'
  # seu$celltype_major[seu$celltype_r2 %in% na.omit(c(unique(pred_list$NK_main$labels)))] <- 'NK_cells'
  print(table(seu$celltype_r2, seu$celltype_major, useNA = 'ifany'))
  # seu <- seu[,!is.na(seu$celltype_r2)]
  seu$celltype_main <- 'celltype'
  seu$celltype_main[str_detect(seu$celltype_r2, 'CD4')] <- 'CD4+T'
  seu$celltype_main[str_detect(seu$celltype_r2, 'CD8')] <- 'CD8+T'
  seu$celltype_main[str_detect(seu$celltype_r2, 'gdT')] <- 'CD8+T'
  seu$celltype_main[seu$celltype_r2 %in% na.omit(c(unique(pred_list$NK_main$pruned.labels)))] <- 'NK'
  seu$celltype_main[seu$celltype_r2 %in% na.omit(c(unique(pred_list$`pan-B`$pruned.labels), unique(pred_list$GCB$pruned.labels), unique(pred_list$Plasma$pruned.labels)))] <- 'B'
  seu$celltype_main[str_detect(seu$celltype_r2, 'PC')] <- 'Plasma'
  seu$celltype_main[seu$celltype_r2 == 'cycling_ASC'] <- 'Plasma'
  seu$celltype_main[str_detect(seu$celltype_major, 'pDC')] <- 'pDC'
  seu$celltype_main[str_detect(seu$celltype_r2, 'cDC')] <- 'cDC'
  seu$celltype_main[str_detect(seu$celltype_r2, 'Mono')] <- 'Mono'
  seu$celltype_main[str_detect(seu$celltype_r2, 'Macro')] <- 'Macro'
  seu$celltype_main[seu$celltype_r2 == 'Mast'] <- 'Mast'
  if (cd45_sorted == F){
    seu$celltype_main[seu$celltype_r2 %in% na.omit(c(unique(pred_list$`pan-Endo`$pruned.labels), unique(pred_list$`pan-Endo_vas`$pruned.labels)))] <- 'Endo'
    seu$celltype_main[str_detect(seu$celltype_r2, 'CAF') | str_detect(seu$celltype_r2, 'Myofibroblast')] <- 'CAF'
    malignant_df <- read.csv(paste0('data/', .cohort, '/infercnv/infercnv_output.csv'), row.names = 1)
    seu$malignant <- malignant_df[colnames(seu),'Malignant']
    seu$celltype_r2[seu$malignant == 'yes'] <- 'Malignant'
    seu$celltype_main[seu$celltype_r2 == 'Malignant'] <- 'Malignant'
  }
  print(table(seu$celltype_major,seu$celltype_main, useNA = 'ifany'))
  return(seu)
}


# TME
## SKCM_Becker BRCA_Bassez1 BRCA_Bassez2 HNSC_Franken CRC_Li CRC_Chen PCa_Hawley
cohort <- 'SKCM_Becker'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2[seu$celltype_major == 'Melanocytes'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Melanocytes'] <- 'Malignant'
DimPlot(seu, group.by = 'celltype_main', cols = getPalette(length(unique(seu$celltype_main))), label = T) /
  DotPlot(seu, group.by = 'celltype_main', features = genes_to_check) + RotatedAxis()
seu@meta.data |> tabyl(celltype_major,celltype_main)
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- subset(seu, subset = celltype_major != 'Neutrophils')
# seu <- subset(seu, subset = celltype_r2 != 'CD56highCD16high')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'SKCM_Plozniak'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2[seu$celltype_major == 'Melanocytes'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Melanocytes'] <- 'Malignant'
seu@meta.data |> tabyl(celltype_major,celltype_main)
DimPlot(seu, group.by = 'celltype_main', cols = getPalette(length(unique(seu$celltype_main))), label = T) /
  DotPlot(seu, group.by = 'celltype_main', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major == 'Epithelial cells', invert = T)
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- subset(seu, subset = celltype_major != 'Neutrophils')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'BRCA_Bassez1'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
DimPlot(seu, group.by = 'celltype_main', cols = getPalette(length(unique(seu$celltype_main))), label = T) /
  DotPlot(seu, group.by = 'celltype_main', features = genes_to_check) + RotatedAxis()
seu@meta.data |> tabyl(celltype_major,celltype_main)
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'BRCA_Bassez2'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
DimPlot(seu, group.by = 'celltype_main', cols = getPalette(length(unique(seu$celltype_main))), label = T) /
  DotPlot(seu, group.by = 'celltype_main', features = genes_to_check) + RotatedAxis()
seu@meta.data |> tabyl(celltype_major,celltype_main)
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'TNBC_Shiao'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
DimPlot(seu, group.by = 'celltype_main', cols = getPalette(length(unique(seu$celltype_main))), label = T) /
  DotPlot(seu, group.by = 'celltype_main', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- subset(seu, subset = celltype_major != 'Neutrophils')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'HNSC_Franken'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
DimPlot(seu, group.by = 'celltype_main', cols = getPalette(length(unique(seu$celltype_main))), label = T) /
  DotPlot(seu, group.by = 'celltype_main', features = genes_to_check) + RotatedAxis()
seu@meta.data |> tabyl(celltype_major,celltype_main)
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- subset(seu, subset = celltype_major != 'Neutrophils')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'CRC_Li'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
DimPlot(seu, group.by = 'celltype_main', cols = getPalette(length(unique(seu$celltype_main))), label = T) /
  DotPlot(seu, group.by = 'celltype_main', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'CRC_Chen'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF'),
                     label_overlapped = c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- subset(seu, subset = celltype_major != 'Neutrophils')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'PCa_Hawley'
pred_res <- c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF','pan-B','Plasma')
seu <- qs_read('data/PCa_Hawley/seu_r1.qs2')
seu$celltype_major <- as.character(seu$celltype_major)
# PCa_Hawley <- subset(PCa_Hawley, subset = celltype_major %in% c('Epithelial cells'), invert = T)
pred_list <- lapply(paste0('data/PCa_Hawley/', pred_res, ".qs2"), qs_read)
names(pred_list) <- pred_res
# adding CD8+ stressed T cell
pred_list$CD8_science_2021$pruned.labels[pred_list$CD8_nm_2023$pruned.labels == 'CD8_c4_Tstr'] <- 'CD8_stress'
pred_list$CD8_nm_2023 <- NULL
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('data/', cohort, '/barcode_gdT.csv')) 
if (nrow(barcode_gdT) > 0) {
  pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
  pred_concat <- rbind(pred_concat, pred_gdT)
} else {
  message("barcode_gdT is empty. Skipping processing for gdT.")
}
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2 <- pred_concat$celltype_r2[match(colnames(seu), pred_concat$cell)]
seu$celltype_r2[seu$celltype_major == 'pDC'] <- 'pDC'
seu$celltype_r2[seu$celltype_major == 'Mast'] <- 'Mast'
seu$celltype_major[seu$celltype_r2 %in% na.omit(c(unique(pred_list$CD4_nm_2023$labels)), unique(pred_list$CD8_science_2021$labels), 'gdT')] <- 'T_cells'
seu$celltype_major[seu$celltype_r2 %in% na.omit(c(unique(pred_list$NK_main$labels)))] <- 'NK_cells'
print(table(seu$celltype_r2, seu$celltype_major, useNA = 'ifany'))
# seu <- seu[,!is.na(seu$celltype_r2)]
seu$celltype_main <- 'celltype'
seu$celltype_main[str_detect(seu$celltype_r2, 'CD4')] <- 'CD4+T'
seu$celltype_main[str_detect(seu$celltype_r2, 'CD8')] <- 'CD8+T'
seu$celltype_main[str_detect(seu$celltype_r2, 'gdT')] <- 'CD8+T'
seu$celltype_main[seu$celltype_r2 %in% na.omit(c(unique(pred_list$NK_main$pruned.labels)))] <- 'NK'
seu$celltype_main[seu$celltype_r2 %in% na.omit(c(unique(pred_list$`pan-B`$pruned.labels), unique(pred_list$GCB$pruned.labels), unique(pred_list$Plasma$pruned.labels)))] <- 'B'
seu$celltype_main[str_detect(seu$celltype_r2, 'PC')] <- 'Plasma'
seu$celltype_main[seu$celltype_r2 == 'cycling_ASC'] <- 'Plasma'
seu$celltype_main[str_detect(seu$celltype_major, 'pDC')] <- 'pDC'
seu$celltype_main[str_detect(seu$celltype_r2, 'cDC')] <- 'cDC'
seu$celltype_main[str_detect(seu$celltype_r2, 'Mono')] <- 'Mono'
seu$celltype_main[str_detect(seu$celltype_r2, 'Macro')] <- 'Macro'
seu$celltype_main[seu$celltype_r2 == 'Mast'] <- 'Mast'
seu$celltype_main[seu$celltype_r2 %in% na.omit(c(unique(pred_list$`pan-Endo`$pruned.labels), unique(pred_list$`pan-Endo_vas`$pruned.labels)))] <- 'Endo'
seu$celltype_main[str_detect(seu$celltype_r2, 'CAF') | str_detect(seu$celltype_r2, 'Myofibroblast')] <- 'CAF'
malignant_df <- read.csv(paste0('data/PCa_Hawley/infercnv/infercnv_output.csv'), row.names = 1)
seu$malignant <- malignant_df[colnames(seu),'Malignant']
seu$celltype_r2[seu$malignant == 'yes'] <- 'Malignant'
seu$celltype_main[seu$celltype_r2 == 'Malignant'] <- 'Malignant'
seu$celltype_r2[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
table(seu$celltype_major, seu$celltype_main, useNA = 'ifany')
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

# CD45+Sorted 
# TNBC_Zhang TNBC_Shiao HNSC_IMCISION HNSC_Luoma 
cohort <- 'TNBC_Zhang'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro'),
                     label_overlapped = c('cDC2','Macro', 'Mono'),
                     cd45_sorted = T)
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- subset(seu, subset = celltype_major != 'Neutrophils')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'HNSC_vanderLeun'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro'),
                     label_overlapped = c('cDC2','Macro', 'Mono'),
                     cd45_sorted = T)
DimPlot(seu, group.by = 'celltype_main', cols = getPalette(length(unique(seu$celltype_main))), label = T) /
  DotPlot(seu, group.by = 'celltype_main', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- subset(seu, subset = celltype_major != 'Neutrophils')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'HNSC_Luoma'
seu <- add_to_seu_r2(.cohort = cohort, 
                     pred_res = c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro'),
                     label_overlapped = c('cDC2','Macro', 'Mono'),
                     cd45_sorted = T)
DimPlot(seu, group.by = 'celltype_main', cols = getPalette(length(unique(seu$celltype_main))), label = T) /
  DotPlot(seu, group.by = 'celltype_main', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- subset(seu, subset = celltype_major != 'Neutrophils')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

# Mixed BCC_Yost
cohort <- 'BCC_Yost'
pred_res <- c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023','NK_main','pan-B','GCB','Plasma','Myeloids_major','cDC2','Mono','Macro','pan-Endo','pan-Endo_vas','CAF','sorted_CD4_nm_2023','sorted_CD8_science_2021','sorted_CD8_nm_2023')
seu <- qs_read('data/BCC_Yost/seu_r1.qs2')
seu$celltype_major <- as.character(seu$celltype_major)
# PCa_Hawley <- subset(PCa_Hawley, subset = celltype_major %in% c('Epithelial cells'), invert = T)
pred_list <- lapply(paste0('data/BCC_Yost/', pred_res, ".qs2"), qs_read)
names(pred_list) <- pred_res
# adding CD8+ stressed T cell
pred_list$CD8_science_2021$pruned.labels[pred_list$CD8_nm_2023$pruned.labels == 'CD8_c4_Tstr'] <- 'CD8_stress'
pred_list$CD8_nm_2023 <- NULL
pred_list$sorted_CD8_science_2021$pruned.labels[pred_list$sorted_CD8_nm_2023$pruned.labels == 'CD8_c4_Tstr'] <- 'CD8_stress'
pred_list$sorted_CD8_nm_2023 <- NULL
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('data/BCC_Yost/barcode_gdT.csv'))
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
barcode_gdT_sorted <- read.csv(paste0('data/BCC_Yost/barcode_gdT_sorted.csv'))
pred_gdT_sorted <- data.frame(cell = barcode_gdT_sorted[,'x'], celltype_r2 = 'gdT')
# pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
# pred_concat <- rbind(pred_concat, pred_gdT)
# # su009 su012
# pred_res <- c('sorted_CD4','sorted_CD8')
# pred_list2 <- lapply(paste0('data/', cohort, '/', cell_subtype, ".qs2"), qread)
# names(pred_list2) <- cell_subtype
# pred_concat2 <- do.call(rbind, lapply(pred_list2, function(df) {
#   data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
# }))
# barcode_gdT2 <- read.csv(paste0('data/', cohort, '/barcode_gdT_sorted.csv')) 
# pred_gdT2 <- data.frame(cell = barcode_gdT2[,'x'], celltype_r2 = 'gdT')
pred_concat <- purrr::reduce(list(pred_concat, pred_gdT, pred_gdT_sorted), rbind)
pred_concat <- filter(pred_concat, !celltype_r2 %in% c('cDC2','Macro', 'Mono','vascular'))
seu$celltype_r2 <- pred_concat$celltype_r2[match(colnames(seu), pred_concat$cell)]
seu$celltype_r2[seu$celltype_major == 'pDC'] <- 'pDC'
seu$celltype_r2[seu$celltype_major == 'Mast'] <- 'Mast'
seu@meta.data |> tabyl(celltype_r2, celltype_major)
seu$celltype_main <- 'celltype'
seu$celltype_main[str_detect(seu$celltype_r2, 'CD4')] <- 'CD4+T'
seu$celltype_main[str_detect(seu$celltype_r2, 'CD8')] <- 'CD8+T'
seu$celltype_main[str_detect(seu$celltype_r2, 'gdT')] <- 'CD8+T'
seu$celltype_main[str_detect(seu$celltype_major, 'NK')] <- 'NK'
seu$celltype_main[seu$celltype_r2 %in% na.omit(c(unique(pred_list$`pan-B`$pruned.labels), unique(pred_list$GCB$pruned.labels), unique(pred_list$Plasma$pruned.labels)))] <- 'B'
seu$celltype_main[str_detect(seu$celltype_r2, 'Plasma') & !seu$patient %in% c('su009', 'su012')] <- 'Plasma'
seu$celltype_main[str_detect(seu$celltype_major, 'pDC')] <- 'pDC'
seu$celltype_main[str_detect(seu$celltype_r2, 'cDC')] <- 'cDC'
seu$celltype_main[str_detect(seu$celltype_r2, 'Mono')] <- 'Mono'
seu$celltype_main[str_detect(seu$celltype_r2, 'Macro')] <- 'Macro'
seu$celltype_main[seu$celltype_r2 %in% na.omit(c(unique(pred_list$`pan-Endo`$pruned.labels), unique(pred_list$`pan-Endo_vas`$pruned.labels)))] <- 'Endo'
seu$celltype_main[str_detect(seu$celltype_major, 'Mast')] <- 'Mast'
seu$celltype_main[str_detect(seu$celltype_r2, 'CAF') | str_detect(seu$celltype_r2, 'Myofibroblast')] <- 'CAF'
# malignant_df <- read.csv(paste0('data/BCC_Yost/infercnv/infercnv_output.csv'), row.names = 1)
# seu$malignant <- malignant_df[colnames(seu),'Malignant']
# seu$celltype_r2[seu$malignant == 'yes'] <- 'Malignant'
# seu$celltype_main[seu$celltype_r2 == 'Malignant'] <- 'Malignant'
seu$celltype_r2[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu$celltype_main[seu$celltype_major == 'Epithelial cells'] <- 'Malignant'
seu@meta.data |> tabyl(celltype_main, celltype_major)
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

# SCC_Yost
cohort <- 'SCC_Yost'
pred_res <- c('sorted_CD4_nm_2023','sorted_CD8_science_2021','sorted_CD8_nm_2023')
seu <- qs_read('data/SCC_Yost/seu_r1.qs2')
seu$celltype_major <- as.character(seu$celltype_major)
pred_list <- lapply(paste0('data/SCC_Yost/', pred_res, ".qs2"), qs_read)
names(pred_list) <- pred_res
# adding CD8+ stressed T cell
pred_list$CD8_science_2021$pruned.labels[pred_list$CD8_nm_2023$pruned.labels == 'CD8_c4_Tstr'] <- 'CD8_stress'
pred_list$CD8_nm_2023 <- NULL
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
barcode_gdT <- read.csv(paste0('data/SCC_Yost/barcode_gdT_sorted.csv'))
pred_gdT <- data.frame(cell = barcode_gdT[,'x'], celltype_r2 = 'gdT')
pred_concat <- rbind(pred_concat, pred_gdT)
seu$celltype_r2 <- pred_concat$celltype_r2[match(colnames(seu), pred_concat$cell)]
seu$celltype_main <- 'celltype'
seu$celltype_main[str_detect(seu$celltype_r2, 'CD4')] <- 'CD4+T'
seu$celltype_main[str_detect(seu$celltype_r2, 'CD8')] <- 'CD8+T'
seu$celltype_main[seu$celltype_r2 == 'gdT'] <- 'CD8+T'
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

cohort <- 'NSCLC_Liu'
pred_res <- c('CD4_nm_2023','CD8_science_2021','CD8_nm_2023')
seu <- qs_read('data/NSCLC_Liu/seu_r1.qs2')
seu$celltype_major <- as.character(seu$celltype_major)
pred_list <- lapply(paste0('data/NSCLC_Liu/', pred_res, ".qs2"), qs_read)
names(pred_list) <- pred_res
# adding CD8+ stressed T cell
pred_list$CD8_science_2021$pruned.labels[pred_list$CD8_nm_2023$pruned.labels == 'CD8_c4_Tstr'] <- 'CD8_stress'
pred_list$CD8_nm_2023 <- NULL
pred_concat <- do.call(rbind, lapply(pred_list, function(df) {
  data.frame(cell = rownames(df), celltype_r2 = df$pruned.labels)
}))
seu$celltype_r2 <- pred_concat$celltype_r2[match(colnames(seu), pred_concat$cell)]
seu$celltype_main <- 'celltype'
seu$celltype_main[str_detect(seu$celltype_r2, 'CD4')] <- 'CD4+T'
seu$celltype_main[str_detect(seu$celltype_r2, 'CD8')] <- 'CD8+T'
seu <- subset(seu, subset = celltype_main != 'celltype')
seu <- seu |> 
  GetAssayData(assay = 'RNA', layer = 'counts') |> 
  CreateSeuratObject(meta.data = seu@meta.data)
qs_save(seu, paste0('data/', cohort, '/seu_r2.qs2'))

# Adjust labels
datasets <- c('SKCM_Becker', 'SKCM_Plozniak', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'BCC_Yost', 'SCC_Yost', 'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 'NSCLC_Liu', 'CRC_Li', 'CRC_Chen', 'PCa_Hawley', 'TNBC_Shiao')
lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) 
  # Modification
  # CD4+T
  # seu$celltype_r2[seu$celltype_r2 == 'CD4_Tn_IL7R-'] <- 'CD4_Prolif'
  # seu$celltype_r2[seu$celltype_r2 == 'CD4_Tn_ADSL'] <- 'CD4_Naive'
  # seu$celltype_r2[seu$celltype_r2 == 'CD4_pre-Tfh_CXCR5+'] <- 'CD4_pre-Tfh_CXCR5'
  # seu$celltype_r2[seu$celltype_r2 == 'CD4_Tm_CAPG+CREM-'] <- 'CD4_Tm_CAPG'
  # seu$celltype_r2[seu$celltype_r2 == 'CD4_Tm_TNF'] <- 'CD4_Tm_CREM-'
  # seu$celltype_r2[seu$celltype_r2 %in% c('CD4_Treg_TNFRSF9-', 'CD4_Treg_S1PR1')] <- 'CD4_Treg_Early'
  # seu$celltype_r2[seu$celltype_r2 == 'CD4_Th_ISG'] <- 'CD4_Th_ISG15'
  # seu$celltype_r2[seu$celltype_r2 == 'CD4_Treg_ISG'] <- 'CD4_Treg_ISG15'
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
  # CAF
  seu$celltype_r2[str_detect(seu$celltype_r2, 'EndMT')] <- 'EndMT'
  # Endo
  seu$celltype_r2[seu$celltype_r2 == 'lymphatics'] <- 'Endo_lymphatic'
  seu$celltype_r2[seu$celltype_r2 == 'arteries'] <- 'Endo_artery'
  seu$celltype_r2[seu$celltype_r2 == 'capillaries'] <- 'Endo_capillary'
  seu$celltype_r2[seu$celltype_r2 == 'tip cell'] <- 'Endo_tip'
  seu$celltype_r2[seu$celltype_r2 == 'veins'] <- 'Endo_vein'
  
  seu$time_point <- ifelse(seu$time_point == 'Post', 'On', seu$time_point)
  
  qsave(seu, paste0('data/', dataset, '/seu_r2.qs2')) 
  if (!dir.exists(paste0('data/', dataset, '/seu_r2/'))){
    dir.create(paste0('data/', dataset, '/seu_r2/'))
  }
  Export10X(seu, dir =paste0('data/', dataset, '/seu_r2/'), 
            append_reductions = NULL, gzip = F)
})





