pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','SingleR','ggsci','tibble','qs','qs2','BiocParallel','Matrix','SingleCellExperiment','scran','parallel','scGate','ggplot2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Reference Preparation
# Myeloids (normalized)
directory <- 'data/GSE154763/'
file_list <- list.files(directory, pattern="expression.csv$")[-1]
expression_list <- lapply(file_list, function(x){t(read.csv(paste0(directory, x), row.names = 1))})
meta_list <- list.files(directory, pattern="metadata.csv$")[-1]
metadata_list <- lapply(meta_list, function(x){read.csv(paste0(directory, x), row.names = 1)})

ref_list <- c()
for (i in 1:length(file_list)){
  print(i)
  sce <- SingleCellExperiment(assays = list(logcounts = expression_list[[i]]), colData = metadata_list[[i]])
  sce$subtype <- sapply(strsplit(sce$MajorCluster,'_'), function(x){paste0(x[2],'_',x[3])})
  sce$majortype <- sapply(strsplit(sce$MajorCluster,'_'), function(x){x[2]})
  # sce <- sce[, !sce$subtype %in% c('Mast_KIT','pDC_LILRA4')]
  ref_list <- c(ref_list, sce)
}
ref_major <- lapply(ref_list, function(x){aggregateReference(x, x$majortype)})
qs_save(ref_major, "data/Ref_SingleR/Myeloids_major_ref.qs2")
ref_mono <- list()
for (i in 1:length(ref_list)){
  sce <- ref_list[[i]]
  sce <- sce[, str_split(sce$subtype, '_', simplify = T)[, 1] == 'Mono']
  sce <- aggregateReference(sce, sce$subtype)
  ref_mono <- c(ref_mono, sce)
}
qsave(ref_mono, "data/Ref_SingleR/Myeloids_mono_ref.qs")
ref_macro <- list()
for (i in 1:length(ref_list)){
  sce <- ref_list[[i]]
  sce <- sce[, str_split(sce$subtype, '_', simplify = T)[, 1] == 'Macro']
  sce <- aggregateReference(sce, sce$subtype)
  ref_macro <- c(ref_macro, sce)
}
qsave(ref_macro, "data/Ref_SingleR/Myeloids_macro_ref.qs")
# cDC2
expr <- list.files(directory, pattern="expression.csv$")[1]
expr <- t(read.csv(paste0(directory, expr), row.names = 1))
meta <- list.files(directory, pattern="metadata.csv$")[1]
meta <- read.csv(paste0(directory, meta),row.names = 1)
sce <- SingleCellExperiment(assays = list(logcounts = expr), colData = meta) 
sce$subtype <- sapply(strsplit(sce$MajorCluster,'_'), function(x){paste0(x[2],'_',x[3])})
ref_cdc2 <- list()
for (i in 1:length(unique(sce$cancer))){
  sce_sub <- sce[, sce$cancer == unique(sce$cancer)[i]]
  sce_sub <- aggregateReference(sce_sub, sce_sub$subtype)
  ref_cdc2 <- c(ref_cdc2, sce_sub)
}
qsave(ref_cdc2, "data/Ref_SingleR/Myeloids_cdc2_ref.qs")

#T/NK (counts)
file_list <- list.files("data/GSE156728(ref_Tcell)", pattern = 'counts', full.names = T)
cancertype <- unique(str_split(file_list,'_',simplify = T)[,3])
df_metadata <- data.table::fread('data/GSE156728(ref_Tcell)/GSE156728_metadata.txt.gz') %>%
  filter(cancerType %in% cancertype, loc=='T') %>% 
  column_to_rownames(var = 'cellID')
df_metadata$meta.cluster <- str_replace_all(df_metadata$meta.cluster, '\\.', '_')
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c01_Tn_TCF7'] <- 'CD4_c01_Naive'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c02_Tn_PASK'] <- 'CD4_c02_pre-Tfh_CXCR5+'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c03_Tn_ADSL'] <- 'CD4_c03_Tn_ADSL'
# df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c04_Tn_il7r'] <- 'CD4_c04_Tn_IL7R-'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c06_Tm_ANXA1'] <- 'CD4_c06_Tm_AREG'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c07_Tm_ANXA2'] <- 'CD4_c07_Tm_TIMP1'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c11_Tm_GZMA'] <- 'CD4_c11_Tm_CAPG+CREM-'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c14_Th17_SLC4A10'] <- 'CD4_c14_Th17_CCR6'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c15_Th17_IL23R'] <- 'CD4_c15_Th17_IL26'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c17_TfhTh1_CXCL13'] <- 'CD4_c17_TfhTh1_IFNG'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c18_Treg_RTKN2'] <- 'CD4_c18_Treg_TNFRSF9-'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c21_Treg_OAS1'] <- 'CD4_c21_Treg_ISG'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD4') & df_metadata$meta.cluster == 'CD4_c22_ISG_IFIT1'] <- 'CD4_c22_Th_ISG'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c01_Tn_MAL'] <- 'CD8_c01_Naive'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c03_Tm_RPS12'] <- 'uncharacterized'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c04_Tm_CD52'] <- 'CD8_c04_Tm_ZNF683'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c05_Tem_CXCR5'] <- 'CD8_c05_Tem_Early'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c08_Tk_TYROBP'] <- 'CD8_c08_NK-like_EOMES'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c09_Tk_KIR2DL4'] <- 'CD8_c09_NK-like_TXK'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c11_Tex_PDCD1'] <- 'CD8_c11_Tex_GZMK'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c13_Tex_myl12a'] <- 'CD8_c13_Tex_OXPHOS-'
df_metadata$meta.cluster[str_detect(df_metadata$meta.cluster, 'CD8') & df_metadata$meta.cluster == 'CD8_c15_ISG_IFIT1'] <- 'CD8_c15_ISG'
df_metadata$meta.cluster[df_metadata$meta.cluster != 'uncharacterized'] <- str_replace(df_metadata$meta.cluster[df_metadata$meta.cluster != 'uncharacterized'], 'c\\d{2}\\_', '')
df_metadata <- filter(df_metadata, !meta.cluster %in% c('uncharacterized', 'CD4_Mix_NME1', 'CD4_Mix_NME2','CD4_c04_Tn_il7r','CD8_Tm_NME1','CD8_Tex_OXPHOS-','CD8_Tm_ZNF683'))

ref_CD4 <- lapply(file_list[str_detect(file_list, 'CD4')], function(x){
  cancertype <- unique(str_split(x,'_',simplify = T)[,3])
  matrix_meta <- filter(df_metadata, cancerType == cancertype, str_detect(meta.cluster, 'CD4')) 
  matrix_count <- data.table::fread(x) %>% tibble::column_to_rownames(var = 'V1') %>% as.sparse()
  matrix_count <- matrix_count[, which(colnames(matrix_count) %in% rownames(matrix_meta))]
  sce <- SingleCellExperiment(assay = list(counts = matrix_count))
  sce$major <- 'CD4'
  return(sce)
})

ref_CD8 <- lapply(file_list[str_detect(file_list, 'CD8')], function(x){
  cancertype <- unique(str_split(x,'_',simplify = T)[,3])
  print(cancertype)
  matrix_meta <- filter(df_metadata, cancerType == cancertype, str_detect(meta.cluster, 'CD8')) 
  matrix_count <- data.table::fread(x) %>% tibble::column_to_rownames(var = 'V1') %>% as.sparse()
  matrix_count <- matrix_count[, which(colnames(matrix_count) %in% rownames(matrix_meta))]
  sce <- SingleCellExperiment(assay = list(counts = matrix_count)) 
  sce$major <- 'CD8'
  return(sce)
})
ref_list_T <- list()
for (i in 1:5){
  print(i)
  common_genes <- intersect(rownames(ref_CD4[[i]]), rownames(ref_CD8[[i]]))
  ref_T <- cbind(ref_CD4[[i]][common_genes,], ref_CD8[[i]][common_genes,])
  ref_list_T[[i]] <- ref_T
}

ref_T <- lapply(ref_list_T, function(sce){
  sce <- logNormCounts(sce) 
  sce <-aggregateReference(sce, sce$major)
  return(sce)
})
qs_save(ref_T, "data/Ref_SingleR/T_ref.qs2")

# Fine level
seu_CD4 <- readRDS('data/CD4_CD8_NM/CD4.rds')
Idents(seu_CD4) <- seu_CD4$cell.type
seu_CD4 <- subset(seu_CD4, subset = TissueType %in% c('Metastatic tumor tissue', 'Primary tumor tissue'))
seu_CD4$cell.type <- as.character(seu_CD4$cell.type)
seu_CD4$cell.type[seu_CD4$cell.type == 'CD4_c2_Tn'] <- 'CD4_RP'
seu_CD4$cell.type[str_detect(seu_CD4$cell.type, '_Tn')] <- 'CD4_naive'
sce <- seu_CD4 |> as.SingleCellExperiment() 
sce <- aggregateReference(sce, sce$cell.type)
qs_save(sce, "data/Ref_SingleR/T_CD4_ref.qs2")

# CD8
# Science 2021
ref_CD8 <- lapply(file_list[str_detect(file_list, 'CD8')], function(x){
  cancertype <- unique(str_split(x,'_',simplify = T)[,3])
  matrix_meta <- filter(df_metadata, cancerType == cancertype, str_detect(meta.cluster, 'CD8')) 
  matrix_count <- data.table::fread(x) %>% tibble::column_to_rownames(var = 'V1') %>% as.sparse()
  matrix_count <- matrix_count[, which(colnames(matrix_count) %in% rownames(matrix_meta))]
  seu <- CreateSeuratObject(counts = matrix_count, meta.data = matrix_meta, min.cells = 5, min.features = 400)
  return(seu)
})
ref_CD8 <- merge(x=ref_CD8[[1]], y=ref_CD8[2:length(ref_CD8)]) |> JoinLayers() |> NormalizeData()
ref_CD8$meta.cluster[ref_CD8$meta.cluster %in% c('CD8_NK-like_TXK','CD8_NK-like_EOMES')] <- 'CD8_NK-like'
sce_sci <- ref_CD8 |> as.SingleCellExperiment()
sce_sci <- sce_sci |> aggregateReference(sce_sci$meta.cluster)

# Nature Medicine 2023 (T stress)
seu_CD8 <- readRDS('data/CD4_CD8_NM/CD8.rds')
Idents(seu_CD8) <- seu_CD8$cell.type
seu_CD8 <- subset(seu_CD8, subset = TissueType %in% c('Metastatic tumor tissue', 'Primary tumor tissue'))
sce_nm <- seu_CD8 |> as.SingleCellExperiment()
sce_nm <- sce_nm |> aggregateReference(sce_nm$cell.type)

sce_list <- list(sce_sci, sce_nm)
qs_save(sce_list, "data/Ref_SingleR/T_CD8_ref.qs2")

# NK
data_dir <- 'data/GSE212890(ref_NK)/'
list.files(data_dir) 
matrix_count <- readMM(paste0(data_dir, list.files(data_dir)[str_detect(list.files(data_dir) , 'counts')])) 
features <- read.csv(paste0(data_dir, list.files(data_dir)[str_detect(list.files(data_dir) , 'genes')]))
barcodes <- read.csv(paste0(data_dir, list.files(data_dir)[str_detect(list.files(data_dir) , 'barcodes')]))
metadata <- read.csv(paste0(data_dir, list.files(data_dir)[str_detect(list.files(data_dir) , 'metadata')]))
colnames(matrix_count) <- features$X0
rownames(matrix_count) <- barcodes$X0
matrix_count <- as(matrix_count, "CsparseMatrix")

ref_NK <- SingleCellExperiment(assay = list(counts = t(matrix_count)), colData = metadata) 
ref_NK <- ref_NK[, ref_NK$meta_tissue == 'Tumor']
cancertype <- unique(ref_NK$meta_histology)
ref_NK <- lapply(cancertype, function(x) {
  return(ref_NK[, ref_NK$meta_histology == x])
})
ref_main_NK <- lapply(ref_NK, function(x) {
  sce <- logNormCounts(x) %>% aggregateReference(.$Majortype)
  return(sce)
})

qs_save(ref_main_NK, "data/Ref_SingleR/NK_main_ref.qs2")
ref_fine_NK <- lapply(ref_NK, function(x) {
  sce <- logNormCounts(x) %>% aggregateReference(.$celltype)
  return(sce)
})
qs_save(ref_fine_NK, "data/Ref_SingleR/NK_fine_ref.qs2")

# CAFs
# https://doi.org/10.1016/j.ccell.2024.08.020
matrix_count <- data.table::fread('data/GSE246215/GSE246215_Fibroblast_counts.csv.gz') |> tibble::column_to_rownames(var = 'GeneName') |> as.sparse()
metadata <- read.csv('data/GSE246215/GSE246215_Fibroblast_metadata.csv.gz') |> tibble::column_to_rownames(var = 'CellName')
seu <- CreateSeuratObject(counts = matrix_count, meta.data = metadata)
seu$Cluster[seu$Cluster == 'c13'] <- 'CAF-ap'
seu$Cluster[seu$Cluster == 'c08'] <- 'Pericytes'
seu$Cluster[seu$Cluster == 'c04'] <- 'CAF-desmo'
seu$Cluster[seu$Cluster %in% c('c01','c18')] <- 'SMC'
seu$Cluster[seu$Cluster == 'c07'] <- 'iCAF_IL6'
seu$Cluster[seu$Cluster == 'c19'] <- 'iCAF_MMP1'
seu$Cluster[seu$Cluster == 'c11'] <- 'Myofibroblasts'
seu$Cluster[seu$Cluster == 'c05'] <- 'CAF-prog'
seu$Cluster[seu$Cluster == 'c16'] <- 'CAF_SFRP2'
seu <- subset(seu, subset = Cluster %in% c('c09','c17','c02','c10','c15','c12','c06','c03','c20','c14'), invert=T)
Idents(seu) <- seu$Cluster
seu <- NormalizeData(seu)
sce <- seu |> as.SingleCellExperiment() |> logNormCounts() 
sce <- sce |> aggregateReference(sce$Cluster)
qsave(sce, "data/Ref_SingleR/CAF_ref_CC.qs")

# Endothelial cells
seu <- schard::h5ad2seurat('data/pan-cancer_Endo/NSR.panE.exprs.all.h5ad.gz')
metadata <- read.csv('data/pan-cancer_Endo/NSR.panE.obs.meta.csv', row.names = 1)
seu <- CreateSeuratObject(counts = seu@assays$RNA@counts, meta.data = metadata)
seu <- subset(seu, subset = tissue == 'T')
seu$major <- ifelse(seu$cellType.major %in% c('veins', 'hypoxia', 'tip cell', 'arteries', 'capillaries'), 'vascular', seu$cellType.major)
ref_list <- lapply(unique(seu$Dataset), function(dataset){
  print(dataset)
  sce <- subset(seu, subset = Dataset == dataset) |> 
    as.SingleCellExperiment() |> 
    logNormCounts() 
  sce <- sce |> aggregateReference(labels = sce$major)
  return(sce)
})
qsave(ref_list, file = 'data/Ref_SingleR/Endo_major.qs')

# Endo_vascular
ref_list <- lapply(unique(seu$Dataset), function(dataset){
  print(dataset)
  sce <- subset(seu, subset = seu$cellType.major %in% c('veins', 'tip cell', 'arteries', 'capillaries') & Dataset == dataset) |> 
    as.SingleCellExperiment() |> 
    logNormCounts() 
  sce <- sce |> aggregateReference(labels = sce$cellType.major)
  return(sce)
})
qsave(ref_list, file = 'data/Ref_SingleR/Endo_vascular.qs')
# Mural cells
endo_mc <- readRDS('data/ref_Endo/panMC.rds')
endo_mc$dataset <- str_split(endo_mc$DonorID, '_', simplify = T)[,1]
sce <- endo_mc |> 
  as.SingleCellExperiment() |> 
  logNormCounts() 
sce <- sce |> aggregateReference(labels = sce$Anno_Tier1)
qsave(sce, file = 'data/Ref_SingleR/Mural_ref.qs')

# B cells
seu <- readRDS("data/blueprint_B/scRNA_data/panB_scRNA_processed_data.rds")
seu <- CreateSeuratObject(counts = seu@assays$RNA@counts, meta.data = seu@meta.data)
seu$majortype <- ifelse(seu$celltype %in% c("B.14.Plasmablast", "B.15.Plasma cell"), 'Plasma',
                        ifelse(seu$celltype %in% c("B.12.LMO2+LZ_GCB", "B.11.SUGCT+DZ_GCB", "B.13.Cycling_GCB"), 'GCB', 'B'))
ref_B_major <- lapply(unique(seu$dataid), function(x) {
  print(x)
  seu_sub <- subset(seu, subset = dataid == x)
  sce <- SingleCellExperiment(assay = list(counts = seu_sub@assays$RNA$counts), colData = seu_sub@meta.data) |> 
    logNormCounts() 
  sce <- aggregateReference(sce,sce$majortype)
  return(sce)
})
qsave(ref_B_major, "data/Ref_SingleR/Pan_B_major.qs")

ref_plasma <- lapply(unique(seu$dataid), function(x) {
  print(x)
  seu_sub <- subset(seu, subset = dataid == x & majortype == 'Plasma')
  sce <- SingleCellExperiment(assay = list(counts = seu_sub@assays$RNA$counts), colData = seu_sub@meta.data) |> 
    logNormCounts() 
  sce <- aggregateReference(sce,sce$celltype)
  return(sce)
})
qsave(ref_plasma, "data/Ref_SingleR/Plasma.qs")

seu_gcb <- subset(seu, subset = majortype == 'GCB')
ref_GCB <- lapply(unique(seu_gcb$dataid)[table(seu_gcb$dataid)>50], function(x) {
  print(x)
  seu_sub <- subset(seu, subset = dataid == x & majortype == 'GCB')
  sce <- SingleCellExperiment(assay = list(counts = seu_sub@assays$RNA$counts), colData = seu_sub@meta.data) |> 
    logNormCounts() 
  sce <- aggregateReference(sce,sce$celltype)
  return(sce)
})
qsave(ref_GCB, "data/Ref_SingleR/GCB.qs")

seu_b <- subset(seu, subset = majortype == 'B')
ref_B <- lapply(names(table(seu_b$dataid))[table(seu_b$dataid)>50], function(x) {
  print(x)
  seu_sub <- subset(seu, subset = dataid == x & majortype == 'B')
  sce <- SingleCellExperiment(assay = list(counts = seu_sub@assays$RNA$counts), colData = seu_sub@meta.data) |> 
    logNormCounts() 
  sce <- aggregateReference(sce,sce$celltype)
  return(sce)
})
qsave(ref_B, "data/Ref_SingleR/B.qs")

# B cells (GSE233236)
count_matrix <- Read10X('data/GSE233236(ref_B)/', gene.column=1)
metadata <- data.table::fread('data/GSE233236(ref_B)/metadata.tsv.gz')
metadata <- tibble::column_to_rownames(metadata, var = 'V1')
seu <- CreateSeuratObject(counts = count_matrix, meta.data = metadata)
seu <- subset(seu, subset = Tissue == 'Tumor')  
seu$cluster <- str_split(seu$Annotation, '_', simplify = T)[,1]
seu$major <- ifelse(seu$cluster %in% c(paste('c', seq(1,9,1), sep = '0'), c('c10','c11')), 'naive/memory',
                    ifelse(seu$cluster %in% paste0('c', seq(12,15,1)), 'cycling/gc', 'asc'))
seu$dataset <- str_split(seu$SampleID, '\\.', simplify = T)[,1]
seu$Annotation <- str_remove(seu$Annotation, "^(c\\d{2}_)")
ref_list <- lapply(unique(seu$dataset), function(ds){
  print(ds)
  sce <- subset(seu, subset = dataset == ds) |> 
    as.SingleCellExperiment() |> 
    logNormCounts() 
  sce <- sce |> aggregateReference(labels = sce$major)
  return(sce)
})
qsave(ref_list, "data/Ref_SingleR/Pan_B-major.qs")

sce <- subset(seu, subset = major %in% 'naive/memory') |> 
  as.SingleCellExperiment() |> 
  logNormCounts() 
sce <- sce |> aggregateReference(sce$Annotation)
qsave(sce, "data/Ref_SingleR/B_naive-memory.qs")

sce <- subset(seu, subset = major %in% 'cycling/gc') |> 
  as.SingleCellExperiment() |> 
  logNormCounts() 
sce <- sce |> aggregateReference(sce$Annotation)
qsave(sce, "data/Ref_SingleR/B_cycling-gc.qs")

asc <- c("PC_IGHG", "PC_IGHA", "early-PC_RGS13", "early-PC_LTB", "cycling_ASC")
sce <- subset(seu, subset = Annotation %in% asc) |> 
  as.SingleCellExperiment() |> 
  logNormCounts() 
sce <- sce |> aggregateReference(sce$Annotation)
qsave(sce, "data/Ref_SingleR/B_asc.qs")




