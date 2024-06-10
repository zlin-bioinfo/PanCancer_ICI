#!/usr/bin/env Rscript
rm(list=ls())
pkgs <- c('Seurat','DoubletFinder','scDblFinder','tidyr','plyr','dplyr','stringr','SingleR','scGate','ggsci','qs','BiocParallel','scRepertoire','harmony','RColorBrewer','scCustomize','SingleCellExperiment','STACAS')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

# scGate model
scGate_models_DB <- get_scGateDB()
# Loading reference
bped <- celldex::BlueprintEncodeData()
bped <- bped[, bped$label.main %in% c("CD4+ T-cells", "CD8+ T-cells", "NK cells", "B-cells", "DC", "Monocytes", "Macrophages", "Neutrophils", "Endothelial cells", "Epithelial cells", "Myocytes", "Fibroblasts", "Melanocytes")]
# cellcycle genelist
exp.mat <- read.table(file = "/bigdata/zlin/Melanoma_meta/tables/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, 
                      as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes 
data(EnsemblGeneTable.Hs)
# # import cell type signature
# sign_celltype <- read.csv('/bigdata/zlin/Melanoma_meta/tables/celltype_signature.csv') |> as.list() |> lapply(function(x) x[x != ''])
# # noise genes
# gene_rm <- readxl::read_excel('/bigdata/zlin/Melanoma_meta/tables/41586_2022_5400_MOESM3_ESM.xlsx', sheet = '1d_Genes excluded', skip = 2)[,-1]
# gene_rm <- as.vector(t(gene_rm)) |> .[!is.na(.)] |> append('MALAT1')
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
paramSweep_v5 <- function (seu, PCs = 1:10, sct = FALSE, num.cores = 1) 
{
  require(Seurat)
  require(fields)
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                 10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts
    n.real.cells <- ncol(data)
  }
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct, mc.cores = num.cores)
    stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3, 
                      n.real.cells, real.cells, pK, pN, data, orig.commands, 
                      PCs, sct)
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                  sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}
doubletFinder_v5 <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                              sct = FALSE, annotations = NULL) 
{
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets)
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets)
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                             ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbor_types[i, ] <- table(doublet_types1[neighbors - 
                                                        n_real.cells]) + table(doublet_types2[neighbors - 
                                                                                                n_real.cells])
          neighbor_types[i, ] <- neighbor_types[i, ]/sum(neighbor_types[i, 
          ])
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                    1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) = levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", 
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, 
                                                                              ct]
      }
    }
    return(seu)
  }
}
# SKCM_Becker    
sample_list = c()
folder <- c('Combicohort','Monocohort')
for (cohort in folder) {
  input_dir = paste0("/bigdata/zlin/Melanoma/data/Melanoma_Becker/",cohort,"/")
  samples = list.files(paste0("/bigdata/zlin/Melanoma/data/Melanoma_Becker/",cohort,"/"))
  sample_list = c(sample_list, samples)
}
sample_list = sample_list[! sample_list =='E19326']
input_dirdte0(rep(c("/bigdata/zlin/Melanoma/data/Melanoma_Becker/Combicohort/", "/bigdata/zlin/Melanoma/data/Melanoma_Becker/Monocohort/"), each=12),
                  sample_list, "/raw_feature_bc_matrix")
seu_list <- lapply(sample_list, function(sample) {
  count_matrix <- Read10X(input_dir[str_detect(input_dir, sample)])
  seu <- CreateSeuratObject(counts = count_matrix, min.cells=3, min.features=200)
  seu$SampleID <- sample
  return(seu)
})
names(seu_list) <- sample_list
saveRDS(seu_list, file = '/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/raw.rds')

# Preprocessing
seu_list <- readRDS('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/raw.rds')
seu_list <- lapply(seu_list, function(seu){
  print(unique(seu$SampleID))
  seu$percent_mito <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu$percent_ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu <- seu |> NormalizeData() |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress=c('nCount_RNA')) |>
    RunPCA() |>
    RunUMAP(dims = 1:20)
  # identify doublets
  # scDblFinder
  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce, dbr.sd=1)
  seu[['scDblFinder.class']] <- sce$scDblFinder.class
  # DoubletFinder
  sweep.list <- paramSweep_v5(seu, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn |>
    arrange(desc(BCmetric))
  pK <- pK[1, 'pK']
  pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
  nExp <- round(ncol(seu) * 0.08)
  seu <- doubletFinder_v5(seu, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:10)
  seu[['DoubletFinder.class']] <- seu@meta.data[, colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]
  seu[[colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]] <- NULL
  seu[[colnames(seu@meta.data)[grepl("pANN", colnames(seu@meta.data))]]] <- NULL
  seu[['singlet']] <- ifelse(seu$scDblFinder.class=='doublet' & seu$DoubletFinder.class=='Doublet', 'no', 'yes')
  # filter low-quality cells and doublets
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mito < 20 & singlet =='yes')
  seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
  seu <- seu%>%
    NormalizeData() |>
    FindVariableFeatures() |> 
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  # Preliminary annotation
  pred_bped_main <- SingleR(test = as.SingleCellExperiment(seu), ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
  seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
  pred_bped_fine <- SingleR(test = as.SingleCellExperiment(seu), ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
  seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
  # scGate
  seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
  return(seu)
})
qsave(seu_list, file = '/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/list.qs')
seu_list <- qread('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/list.qs')
# check sample quality
df <- data.frame(
  nfeature = sapply(seu_list, function(seu_obj){nfeature <- sum(Matrix::rowSums(seu_obj@assays$RNA$counts>0)>3)
  return(nfeature)}),
  ncell = sapply(seu_list, function(seu_obj){ncell <- ncol(seu_obj)
  return(ncell)})
)
sample_rm <- rownames(filter(df, nfeature < 9000 | ncell < 200)); sample_rm
seu_list <- seu_list[-which(names(seu_list) %in% sample_rm)]
# merge samples
seu <- merge(x=seu_list[[1]], y=seu_list[2:length(seu_list)], add.cell.ids = names(seu_list))
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
# adding clinical data
sample_metadata <- read.csv(file = "/bigdata/melanoma/Metadata.csv", sep = '\t', row.names = 2)
seu$treatment <- sample_metadata[seu$SampleID,]$Treatment
seu$treatment <- ifelse(seu$treatment == 'Combi', 'aPD1+CTLA4', 'aPD1')
seu$time_point <- sample_metadata[seu$SampleID,]$Time.point
seu$time_point[seu$time_point=='after'] <- 'Post'
seu$time_point[seu$time_point=='before'] <- 'Pre'
seu$patient <- sample_metadata[seu$SampleID,]$Patient
seu$patient <- str_replace(seu$patient, 'atient ','')
clinical_metadata <- readxl::read_xlsx('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/clinical_info.xlsx')
seu$response <- clinical_metadata$`ICI response`[match(seu$patient, clinical_metadata$Patient)]
seu$site <- clinical_metadata$`Leision type(primary site/metastasis)`[match(seu$patient, clinical_metadata$Patient)]
seu$dataset <- 'SKCM_Becker'
seu$sample <- paste0(seu$dataset, '_', seu$patient, '_', seu$time_point)
seu$patient <- paste0(seu$dataset, '_', seu$patient)
clinical_metadata <- drop_na(clinical_metadata, Age) 
seu$age <- clinical_metadata$Age[match(seu$patient, clinical_metadata$Patient)]
seu$gender <- clinical_metadata$Gender[match(seu$patient, clinical_metadata$Patient)]
seu$prior <- 'No'
# seu$immune <- 'non-immune'
# seu$immune[seu$celltype_bped_fine %in% immune] <- 'immune'
seu$modality <- ifelse(seu$treatment == 'aPD1+CTLA4', 'Dual', 'Mono')
seu$cancertype <- 'SKCM'
seu$interval <- 7
seu$interval[seu$patient == 'SKCM_Becker_P1'] <- 12
seu$interval[seu$patient %in% c('SKCM_Becker_P2', 'SKCM_Becker_P5', 'SKCM_Becker_P6')] <- 9
seu$response <- ifelse(seu$response %in% c('PR', 'CR'), 'RE', ifelse(seu$response %in% c('SD', 'PD'), 'NR', 'NE'))
seu$response[seu$patient == 'SKCM_Becker_P7'] <- 'NR'
seu$res_metric <- 'RECIST'
seu$res_metric[seu$patient == 'SKCM_Becker_P7'] <- 'Survival'
# remove unpaired samples
table(seu$patient, seu$time_point)
seu <- subset(seu, subset = patient %in% c('SKCM_Becker_P4', 'SKCM_Becker_P8'), invert = T)
my_scGate_model <- gating_model(name = "pDC", signature = c("LILRA4")) 
seu <- scGate(data = seu, model = my_scGate_model)
seu <- seu |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(seu, group.by = 'celltype_bped_main', cols = getPalette(length(unique(seu$celltype_bped_main))))
DimPlot(seu, group.by = 'scGate_multi', cols = getPalette(length(unique(seu$scGate_multi))))

# 'CD3D','CD4','CD8A','KLRD1','CD79A','MS4A1','LILRA4','CD68','COL3A1','ACTA2','PECAM1','MLANA','KRT15','MKI67'
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38','CD24'), # Plasma cells 'CD38','CD24'
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('COL3A1','FAP', 'COL1A1','ATCA2'), 
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('KRT15','KRT17','EPCAM'),
                      c('MKI67')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Fibro','Endo','Mela','Epi','Proliferating')
FeaturePlot_scCustom(seurat_object = seu, colors_use = viridis_magma_dark_high, 
                     features = c('MLANA','MITF')) / 
  DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters)))) / 
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'celltype_bped_main', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'scGate_multi', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$celltype_major == 'Melanocytes'] <- 'Melanoma'
seu$celltype_major[seu$seurat_clusters %in% c('0','1','3','7','9','13','16','22','25')] <- 'Melanoma'
seu$celltype_major[seu$seurat_clusters == '21'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' & seu$celltype_major == 'B-cells'] <- 'Plasma'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'CD4T'] <- 'CD4+ T-cells'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'CD8T'] <- 'CD8+ T-cells'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('unknown', 'Neutrophils') | (scGate_multi == 'Multi'), invert = T)
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/seu_r1.qs')

# Franken
dataset = 'HNSC_Franken'
# dir.create(file.path('output', 'QC', dataset))
matrix_count <- Read10X('/bigdata/zlin/Melanoma/data/additional_datasets/Franken/HNSCC/')
matrix_meta <- read.csv('/bigdata/zlin/Melanoma/data/additional_datasets/Franken/HNSCC_metadata.csv',row.names = 1)
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta)
seu <- subset(seu, subset = SampleType %in% c('On-treatment', 'Pre-treatment'))
seu$Treatment <- mapvalues(seu$Treatment,from = c('Durvalumab','Durvalumab-Tremelimumab'), to = c('aPDL1','aPDL1+CTLA4'))
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='Treatment')] <- "treatment"
seu$SampleType <- mapvalues(seu$SampleType, from = c('Pre-treatment','On-treatment'), to = c('Pre','Post'))
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='SampleType')] <- "time_point"
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='Patient')] <- "patient"
seu$dataset <- dataset
seu$patient <- paste0(seu$dataset, '_', 'P', seu$patient)
seu$sample <- paste0(seu$patient,'_', seu$time_point)
seu$dataset <- dataset
seu$cancertype <- 'HNSC'
seu$res_metric <- 'T-cell expansion&MRI'
seu$subtype <- 'HNSC'
seu$patient <- str_replace(seu$patient, 'BIOKEY', 'BRCA_Bassez1')
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='timepoint')] <- "time_point"
e <- paste0('HNSC_Franken_P', c(1,7,8,9,13,14,16,17,19,20))
ne <- paste0('HNSC_Franken_P', c(2,3,4,10,11,12,15,18))
seu$expansion <- 'NA'
seu$expansion[seu$patient %in% e] <- 'Yes'
seu$expansion[seu$patient %in% ne] <- 'No'
df_mri <- data.frame(pt = paste0('HNSC_Franken_P', 1:20),
                     mri_pre = c(43,37,36,34,28,14,43,43,30,33,44,25,64,41,47,46,25,37,32,49),
                     mri_post = c(42,45,52,43,26,13,39,45,27,35,43,27,55,35,NA,27,25,49,34,47))
seu$mri_pre <- df_mri$mri_pre[match(seu$patient, df_mri$pt)]
seu$mri_post <- df_mri$mri_post[match(seu$patient, df_mri$pt)]
seu$mri_change <- seu$mri_post - seu$mri_pre
seu$mri_fc <- seu$mri_post/seu$mri_pre
seu$response <- 'NR'
seu$response[seu$expansion == 'Yes'] <- 'RE'
seu$res_metric <- 'T-cell expansion'
seu$res_metric[seu$expansion == 'Yes' & seu$mri_change <= 0] <- 'T-cell expansion+MRI'
seu$response[seu$expansion == 'Yes' & seu$mri_change > 0] <- 'NR'
seu$res_metric[seu$expansion == 'Yes' & seu$mri_change > 0] <- 'MRI'
seu$interval <- 14

# seu_list <- SplitObject(seu, split.by = 'sample')
# # check sample quality
# df <- data.frame(
#   nfeature = sapply(seu_list, function(seu){nfeature <- sum(Matrix::rowSums(seu@assays$RNA$counts>0)>3)
#   return(nfeature)}),
#   ncell = sapply(seu_list, function(seu){ncell <- ncol(seu)
#   return(ncell)})
# )
# # remove samples with <200 cells or <9000 genes 
# sample_rm <- rownames(filter(df, nfeature < 9000 | ncell < 200)); sample_rm
# NA

seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) |>
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
sce <- as.SingleCellExperiment(seu)
pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/HNSC_Franken/processing.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_Franken/processing.qs')
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('CXCR1', 'CXCR2', 'PTGS2','OLR1', 'VEGFA'),
                      c('COL3A1','FAP', 'COL1A1','ATCA2'), 
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('KRT15','KRT17','EPCAM'),
                      c('MKI67')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Neu','Fibro','Endo','Mela','Epi','Proliferating')
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '20'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' & seu$celltype_major == 'B-cells'] <- 'Plasma'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'CD4T'] <- 'CD4+ T-cells'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'CD8T'] <- 'CD8+ T-cells'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'Fibroblast'] <- 'Fibroblasts'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'Epithelial'] <- 'Epithelial cells'
FeaturePlot_scCustom(seurat_object = seu, colors_use = viridis_magma_dark_high, 
                     features = c('CD3D','CD4','CD8A','KLRD1','CD79A','MS4A1','LILRA4','CD68','COL3A1','ACTA2','PECAM1','MLANA','KRT15','MKI67'))
FeaturePlot_scCustom(seurat_object = seu, colors_use = viridis_magma_dark_high, 
                     features = c('CD3D'))
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
table(seu$celltype_major)
seu <- subset(seu, subset = (celltype_major %in% c('Neutrophils','unknown','Melanocytes')) | (scGate_multi == 'Multi'), invert = T)
seu$prior <- 'No'
seu$modality <- ifelse(seu$treatment == 'aPDL1+CTLA4', 'Dual', 'Mono')
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/HNSC_Franken/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_Franken/seu_r1.qs')

# BRCA_Bassez anti-PD1 dataset 
dataset = 'BRCA_Bassez1'
matrix_count <- readRDS('/bigdata/zlin/Melanoma/data/additional_datasets/Bassez/1863-counts_cells_cohort1.rds')
matrix_meta <- read.csv('/bigdata/zlin/Melanoma/data/additional_datasets/Bassez/1872-BIOKEY_metaData_cohort1_web.csv',row.names = 1)
matrix_tcr <- read.csv('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/1879-BIOKEY_barcodes_vdj_combined_cohort1.csv')
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta)
seu <- subset(seu, subset = expansion == "n/a", invert = TRUE)
seu[['cdr3s_nt']] <- matrix_tcr$cdr3_nt[match(str_split(colnames(seu),'-', simplify = T)[,1], matrix_tcr$barcode)]
extract_unique_tcr <- function(cdr3s_nt) {
  if (is.na(cdr3s_nt)) return(NA)
  tra <- str_extract_all(cdr3s_nt, "TRA:[^;]*")
  trb <- str_extract_all(cdr3s_nt, "TRB:[^;]*")
  if (length(tra[[1]]) == 1 & length(trb[[1]]) == 1) {
    return(cdr3s_nt)
  } else {
    return(NA)
  }
}
seu[['cdr3s_nt_unique']] <-  sapply(seu$cdr3s_nt, extract_unique_tcr)
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='patient_id')] <- "patient"
seu$patient <- str_replace(seu$patient, 'BIOKEY', 'BRCA_Bassez1')
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='timepoint')] <- "time_point"
seu$time_point <- ifelse(seu$time_point == 'On', 'Post', 'Pre')
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='expansion')] <- "response"
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='BC_type')] <- "subtype"
seu$sample <- paste0(seu$patient,'_', seu$time_point)
seu$treatment <- 'aPD1'
seu$dataset <- 'BRCA_Bassez1'
seu$cancertype <- 'BRCA'
seu$response <- ifelse(seu$response == 'E', 'RE', 'NR')
seu$res_metric <- 'T-cell expansion'
seu$site <- 'n/a'
seu$orig.ident <- NULL
seu$nCount_RNA <- NULL
seu$nFeature_RNA <- NULL
seu_list <- SplitObject(seu, split.by = 'sample')
# check sample quality
df <- data.frame(
  nfeature = sapply(seu_list, function(seu){nfeature <- sum(Matrix::rowSums(seu@assays$RNA$counts>0)>3)
  return(nfeature)}),
  ncell = sapply(seu_list, function(seu){ncell <- ncol(seu)
  return(ncell)})
)
# remove samples with <200 cells or <9000 genes 
sample_rm <- rownames(filter(df, nfeature < 9000 | ncell < 200)); sample_rm
# BRCA_Bassez1_23_Pre
remove(seu_list)
seu <- subset(seu, subset = patient == 'BRCA_Bassez1_23', invert=T)
seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  # CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) |> 
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
sce <- as.SingleCellExperiment(seu)
pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/processing.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/processing.qs')
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('COL3A1','FAP', 'COL1A1','ATCA2'), 
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('KRT15','KRT17','EPCAM'),
                      c('MKI67')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Fibro','Endo','Mela','Epi','Proliferating')
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '25'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' & seu$celltype_major == 'B-cells'] <- 'Plasma'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'CD4T'] <- 'CD4+ T-cells'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'CD8T'] <- 'CD8+ T-cells'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'Fibroblast'] <- 'Fibroblasts'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'Epithelial'] <- 'Epithelial cells'
FeaturePlot_scCustom(seurat_object = seu, colors_use = viridis_magma_dark_high, 
                     features = c('CD3D','CD4','CD8A','KLRD1','CD79A','MS4A1','LILRA4','CD68','COL3A1','ACTA2','PECAM1','MLANA','KRT15','MKI67'))
FeaturePlot_scCustom(seurat_object = seu, colors_use = viridis_magma_dark_high, 
                     features = c('CD3D'))
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
table(seu$celltype_major)
seu <- subset(seu, subset = (celltype_major %in% c('Neutrophils','unknown','Melanocytes')) | (scGate_multi == 'Multi'), invert = T)
seu$age <- df$Age[match(seu$patient, df$ID)]
seu$prior <- 'No'
seu$modality <- ifelse(seu$treatment == 'aPD1+CTLA4', 'Dual', 'Mono')
seu$cancertype <- ifelse(seu$subtype == 'HER2+', 'HER2+BC', ifelse(seu$subtype == 'ER+', 'ER+BC', 'TNBC'))
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r1.qs')

# BRCA_Bassez Chemo+anti-PD1 therapy
matrix_count <- readRDS('/bigdata/zlin/Melanoma/data/additional_datasets/Bassez/1867-counts_cells_cohort2.rds')
matrix_meta <- read.csv('/bigdata/zlin/Melanoma/data/additional_datasets/Bassez/1871-BIOKEY_metaData_cohort2_web.csv',row.names = 1)
matrix_tcr <- read.csv('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/1880-BIOKEY_barcodes_vdj_combined_cohort2.csv')
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta)
seu[['cdr3s_nt']] <- matrix_tcr$cdr3_nt[match(str_split(colnames(seu),'-', simplify = T)[,1], matrix_tcr$barcode)]
extract_unique_tcr <- function(cdr3s_nt) {
  if (is.na(cdr3s_nt)) return(NA)
  tra <- str_extract_all(cdr3s_nt, "TRA:[^;]*")
  trb <- str_extract_all(cdr3s_nt, "TRB:[^;]*")
  if (length(tra[[1]]) == 1 & length(trb[[1]]) == 1) {
    return(cdr3s_nt)
  } else {
    return(NA)
  }
}
seu[['cdr3s_nt_unique']] <-  sapply(seu$cdr3s_nt, extract_unique_tcr)
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='patient_id')] <- "patient"
seu$patient <- str_replace(seu$patient, 'BIOKEY', 'BRCA_Bassez2')
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='timepoint')] <- "time_point"
seu$time_point <- ifelse(seu$time_point == 'On', 'Post', 'Pre')
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='expansion')] <- "response"
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='BC_type')] <- "subtype"
seu$sample <- paste0(seu$patient,'_', seu$time_point)
seu$treatment <- 'aPD1(pre-Chemo)'
seu$dataset <- 'BRCA_Bassez2'
seu$cancertype <- 'BRCA'
seu$response <- ifelse(seu$response == 'E', 'RE', 'NR')
seu$res_metric <- 'T-cell expansion'
seu$site <- 'n/a'
seu$orig.ident <- NULL
seu$nCount_RNA <- NULL
seu$nFeature_RNA <- NULL
seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) |> 
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
sce <- as.SingleCellExperiment(seu)
pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/processing.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/processing.qs')
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('COL3A1','FAP', 'COL1A1','ATCA2'), 
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('KRT15','KRT17','EPCAM'),
                      c('MKI67')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Fibro','Endo','Mela','Epi','Proliferating')
DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'celltype_bped_main', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'scGate_multi', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '22'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' & seu$celltype_major == 'B-cells'] <- 'Plasma'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'Fibroblast'] <- 'Fibroblasts'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
table(seu$celltype_major)
seu <- subset(seu, subset = (celltype_major %in% c('unknown','Melanocytes')) | (scGate_multi == 'Multi'), invert = T)

# clinical
df <- data.frame(
  ID = c(32:42),
  `E/NE` = c("E", "NE", "E", "E", "NE", "NE", "NE", "NE", "NE", "NE", "NE"),
  Age = c("30-40", "51-60", "30-40", "71-80", "41-50", "61-70", "51-60", "51-60", "51-60", "30-40", "61-70"),
  `Histological type` = c("IBC-NST with medullary pattern", "IBC-NST", "IBC-NST", "IBC-NST", "IBC-NST*", "IBC-NST", "IBC-NST", "IBC-NST", "IBC-NST", "IBC-NST", "IBC-NST"),
  Type = c("TNBC", "ER-/PR-/HER2+", "TNBC", "TNBC", "TNBC", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "ER+/PR+/HER2-", "TNBC", "TNBC", "ER+/PR+/HER2-"),
  `Time between biopsy and surgery` = c(12, 12, 12, 8, 9, 7, 9, 9, 9, 8, 8),
  pTNM = c("ypT1N0", "ypT0N0", "ypT2N0", "ypT1N2", "ypT1N1", "ypT2N3", "ypT2N1", "ypT3N3", "ypT1bN0(sn)", "ypT2N1a", "ypT2(m)N1a"),
  cTNM = c("ycT1N0", "ycT2N1", "ycT2N0", "ycT2N1", "ycT1N0", "ycT2N0", "ycT1N0", "ycT2N1", "ycT2(m)N0M0", "ycT2N1M0", "ycT4bN0M0"),
  sTILs = c("ND", "No tumor", "34/35", "No tumor", "15/6", "1/2", "No tumor/3", "2/1", "8/ND", "1.2/0.2", "0.4/0.4"),
  Ki67 = c("Intermediate/Intermediate", "ND/ND", "High/Intermediate", "Low/ND", "High/High", "ND/ND", "Low/Low", "High/Low", "Low/Low", "High/High", "Low/Low"),
  Menopausal = c("Pre", "Post", "Pre", "Post", "Post", "Post", "Pre", "Post", "Post", "Pre", "Post")
)
df$ID <- paste0('BIOKEY_', df$ID)
df <- df[df$ID %in% unique(seu$patient),]
seu$age <- df$Age[match(seu$patient, df$ID)]
seu$interval <- 9
seu$prior <- 'Yes'
seu$modality <- ifelse(seu$treatment == 'aPD1+CTLA4', 'Dual', 'Mono')
seu$cancertype <- ifelse(seu$subtype == 'HER2+', 'HER2+BC', ifelse(seu$subtype == 'ER+', 'ER+BC', 'TNBC'))
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/seu_r1.qs')

# TNBC
input_dir <- "/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/"
count_matrix <- Read10X(input_dir,gene.column=1)
clinical_data <- readxl::read_xlsx("/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/1-s2.0-S1535610821004992-mmc2.xlsx", skip = 1) |>
  janitor::clean_names() |>
  slice(2:23) |>
  tidyr::fill(treatment) 
patient_id <- filter(clinical_data, tumor == 'Y', x11 == 'Y')$patient_id
count_matrix <- count_matrix[,str_split(colnames(count_matrix),'\\.', simplify = T)[,2] %in% paste0(rep(c('Pre_','Post_'),each=12),rep(paste0(patient_id, '_t'),2))]
seu <- CreateSeuratObject(counts = count_matrix, min.cells=3, min.features=200)
seu$sample <- paste0(str_split(str_split(rownames(seu@meta.data),'\\.', simplify = T)[,2],'_',simplify=T)[,1], '_',
                     str_split(str_split(rownames(seu@meta.data),'\\.', simplify = T)[,2],'_',simplify=T)[,2])
seu$time_point <- str_split(seu$sample, '_', simplify = T)[,1]
seu$patient <- str_split(seu$sample, '_', simplify = T)[,2]
seu$response <- clinical_data$clinical_efficacy_number[match( seu$patient,clinical_data$patient_id)]
seu$treatment <- clinical_data$treatment[match( seu$patient,clinical_data$patient_id)]
seu$treatment[seu$treatment == 'Anti-PD-L1+ Chemo'] <- 'aPDL1+Chemo'
seu$dataset <- 'TNBC_Zhang'
seu$patient <- paste0(seu$dataset, '_', seu$patient)
seu$sample <- paste0(seu$dataset, '_', seu$sample)
seu$cancertype <- 'TNBC'
seu$treatment <- as.character(seu$treatment)
seu$interval <- 28
seu$response <- ifelse(seu$response %in% c('PR', 'CR'), 'RE', 'NR')
seu$res_metric <- 'RECIST'
seu$modality <- 'Mono'
seu$prior <- 'No'
seu <- subset(seu, subset = treatment == 'aPDL1+Chemo')
seu_list <- SplitObject(seu, split.by = "sample")
# df <- data.frame(
#   nfeature = sapply(seu_list, function(seu_obj){nfeature <- sum(Matrix::rowSums(seu_obj@assays$RNA$counts>0)>3)
#   return(nfeature)}),
#   ncell = sapply(seu_list, function(seu_obj){ncell <- ncol(seu_obj)
#   return(ncell)})
# )
# sample_rm <- rownames(filter(df, nfeature < 9000 | ncell < 200)); sample_rm # "Post_P016" 
seu <- subset(seu, subset = patient == "TNBC_Zhang_P016", invert=TRUE)
seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) |> 
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
sce <- as.SingleCellExperiment(seu)
pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/processing.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/processing.qs')
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac')
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '17'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' & seu$celltype_major == 'B-cells'] <- 'Plasma'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
table(seu$celltype_major)
seu <- subset(seu, subset = (celltype_major %in% c('Epithelial cells','unknown','Endothelial cells')) | (scGate_multi == 'Multi'), invert = T)
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/seu_r1.qs')

# Yost
dataset = 'BCC_Yost'
matrix_count <- data.table::fread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/GSE123813_bcc_scRNA_counts.txt') |> 
  tibble::column_to_rownames(var = 'V1') |> as.sparse()
matrix_meta <- data.table::fread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/GSE123813_bcc_all_metadata.txt') |> 
  tibble::column_to_rownames(var = 'cell.id') 
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta) |> 
  subset(subset = orig.ident %in% c('bcc.su010.pre.tcell','bcc.su010.post.tcell'), invert = T)
seu$cancertype <- 'BCC'
seu$res_metric <- 'RECIST'
matrix_tcr <- read.table('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/GSE123813_bcc_tcr.txt')
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'treatment')] <- "time_point"
seu$time_point <- mapvalues(seu$time_point, from = c('pre','post'), to = c('Pre','Post'))
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'cluster')] <- "celltype_orig"
seu$dataset <- 'BCC/SCC_Yost'
seu$patient <- paste0(seu$dataset, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$prior <- 'Yes'
seu$prior[seu$patient == 'BCC/SCC_Yost_su004'] <- 'No'
seu$site <- 'n/a'
seu$treatment <- 'aPD1'
seu$modality <- 'Mono'
seu_list <- SplitObject(seu, split.by = 'sample')
# check sample quality
(df <- data.frame(
  nfeature = sapply(seu_list, function(seu_obj){nfeature <- sum(Matrix::rowSums(seu_obj@assays$RNA$counts>0)>3)
  return(nfeature)}),
  ncell = sapply(seu_list, function(seu_obj){ncell <- ncol(seu_obj)
  return(ncell)})
))
sample_rm <- rownames(filter(df, nfeature < 9000 | ncell < 200)); sample_rm
# remove(seu_list)
seu <- subset(seu, subset = patient == 'BCC/SCC_Yost_su002', invert=T)
# # Find those that have both "_Pre" and "_Post" counterparts
# with_both <- unique(seu$patient)[
#   sapply(unique(seu$patient), function(x) {
#     sum(unique(seu$sample) %in% c(paste0(x, "_Pre"), paste0(x, "_Post"))) == 2
#   })
# ]
# # keep paired samples
# seu <- subset(seu, subset = patient %in% with_both)
seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) |> 
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
sce <- as.SingleCellExperiment(seu)
pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/BCC_Yost/processing.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/processing.qs')
tcr <- read.table('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/GSE123813_bcc_tcr.txt')
seu[['cdr3s_nt']] <- tcr$cdr3s_nt[match(rownames(seu@meta.data), rownames(tcr))]
# clinical
clinical <- readxl::read_xlsx('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/41591_2019_522_MOESM2_ESM.xlsx', skip = 2, n_max = 16)
clinical <- clinical[clinical$Patient %in% unique(str_replace(seu$patient, 'BCC/SCC_Yost_','')),]
clinical$`scRNA days pre treatment`[clinical$Patient == 'su001'] <- -78
clinical$`scRNA days post treatment`[clinical$Patient == 'su003'] <- 121
clinical$interval <- as.numeric(clinical$`scRNA days post treatment`) 
seu$interval <- clinical$interval[match(str_replace(seu$patient, 'BCC/SCC_Yost_',''), clinical$Patient)]
seu$response <- clinical$Response[match(str_replace(seu$patient, 'BCC/SCC_Yost_',''), clinical$Patient)]
seu$response <- ifelse(seu$response == 'Yes', 'RE', 'NR')

genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38','CD24'), # Plasma cells 'CD38','CD24'
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('COL3A1','FAP', 'COL1A1','ATCA2'), 
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('KRT15','KRT17','EPCAM'),
                      c('MKI67')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Fibro','Endo','Mela','Epi','Proliferating')
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters)))) / DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'scGate_multi', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'panDC'] <- 'DC'
seu$celltype_major[seu$seurat_clusters == '15'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' & seu$celltype_major == 'B-cells'] <- 'Plasma'
seu$celltype_major[seu$celltype_major == 'unknown' & seu$scGate_multi == 'CD4T'] <- 'Fibroblasts'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
table(seu$celltype_major)
seu <- subset(seu, subset = (celltype_major %in% c('unknown','Melanocytes')) | (scGate_multi == 'Multi'), invert = T)
seu <- subset(seu, subset = (patient %in% c('BCC/SCC_Yost_su009','BCC/SCC_Yost_su012') & (celltype_major %in% c('B-cells', 'Endothelial cells', 'Fibroblasts','Monocytes' ,'Plasma'))), invert = T)
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/BCC_Yost/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/seu_r1.qs')


# SCC_Yost
matrix_count <- data.table::fread('/bigdata/zlin/Melanoma_meta/data/SCC_Yost/GSE123813_scc_scRNA_counts.txt') |>
  tibble::column_to_rownames(var = 'V1') |> as.sparse()
matrix_meta <- data.table::fread('/bigdata/zlin/Melanoma_meta/data/SCC_Yost/GSE123813_scc_metadata.txt') |>
  tibble::column_to_rownames(var = 'cell.id')
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta)
seu$subtype <- 'SCC'
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'treatment')] <- "time_point"
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'timepoint')] <- "time_point"
seu$time_point <- mapvalues(seu$time_point, from = c('pre','post'), to = c('Pre','Post'))
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'cluster')] <- "celltype_orig"
seu$sample <- paste0(seu$patient,'_', seu$time_point)
seu <- subset(seu, subset = patient %in% c('su010', 'su013'), invert = T)
seu$dataset <- 'BCC/SCC_Yost'
seu$cancertype <- 'SCC'
seu$res_metric <- 'RECIST'
seu$site <- 'n/a'
seu$prior <- 'Yes'
seu$treatment <- 'aPD1'
# seu_list <- SplitObject(seu, split.by = 'sample')
# # check sample quality
# (df <- data.frame(
#   nfeature = sapply(seu_list, function(seu_obj){nfeature <- sum(Matrix::rowSums(seu_obj@assays$RNA$counts>0)>3)
#   return(nfeature)}),
#   ncell = sapply(seu_list, function(seu_obj){ncell <- ncol(seu_obj)
#   return(ncell)})
# ))
seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu[['celltype_major']] <- 'T_cells'
tcr <- read.table('/bigdata/zlin/Melanoma_meta/data/SCC_Yost/GSE123813_scc_tcr.txt')
seu[['cdr3s_nt']] <- tcr$cdr3s_nt[match(rownames(seu@meta.data), rownames(tcr))]

#clinical
clinical <- readxl::read_xlsx('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/41591_2019_522_MOESM2_ESM.xlsx', skip = 2, n_max = 16)
clinical <- clinical[clinical$Patient %in% unique(seu$patient),]
clinical$`scRNA days pre treatment`[clinical$Patient == 'su001'] <- -78
clinical$interval <- as.numeric(clinical$`scRNA days post treatment`) 
seu$interval <- clinical$interval[match(seu$patient, clinical$Patient)]
seu$response <- clinical$Response[match(seu$patient, clinical$Patient)]
seu$response <- ifelse(seu$response == 'Yes (CR)', 'RE', 'NR')
seu$patient <- paste0(seu$dataset, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$modality <- 'Mono'
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/SCC_Yost/seu_r1.qs')

# Liu NSCLC
matrix_count <- readRDS('/bigdata/zlin/Melanoma_meta/data/NSCLC_Liu/GSE179994_all.Tcell.rawCounts.rds')
matrix_meta <- data.table::fread('/bigdata/zlin/Melanoma_meta/data/NSCLC_Liu/GSE179994_Tcell.metadata.tsv') |> 
  tibble::column_to_rownames('cellid') 
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta)
matrix_tcr <- read.delim('/bigdata/zlin/Melanoma_meta/data/NSCLC_Liu/GSE179994_all.scTCR.tsv', row.names = 1)
seu$percent_mito <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mito < 20)
seu[['patient']] <- str_split(seu$sample, '\\.', simplify = T)[,1]
seu[['time_point']] <- str_split(seu$sample, '\\.', simplify = T)[,2]
seu_list <- SplitObject(seu, split.by = "sample")
(df <- data.frame(
  nfeature = sapply(seu_list, function(seu_obj){nfeature <- sum(Matrix::rowSums(seu_obj@assays$RNA$counts>0)>3)
  return(nfeature)}),
  ncell = sapply(seu_list, function(seu_obj){ncell <- ncol(seu_obj)
  return(ncell)})
))
# sample_rm <- rownames(filter(df, nfeature < 10000 | ncell < 100)); sample_rm
# rm(seu_list)
pt_keep <- c('P1', 'P10', 'P13', 'P19', 'P29', 'P30', 'P33', 'P35')
seu <- subset(seu, subset = patient %in% pt_keep)
seu <- subset(seu, subset = sample %in% c('P1.post.1', 'P1.post.2', 'P13.post.2'), invert = T)
seu[['dataset']] <-  'NSCLC_Liu'
seu$patient <- paste0(seu$dataset, '_', seu$patient)
seu[['sample']] <- paste0(seu$patient, '_', seu$time_point)
seu[['response']] <- 'RE'
seu[['interval']] <-  30
seu$interval[seu$patient == 'P1'] <- 213-14
seu$interval[seu$patient == 'P10'] <- 61-7
seu$interval[seu$patient == 'P13'] <- 94-17
seu$interval[seu$patient == 'P19'] <- 62-14
seu$interval[seu$patient == 'P29'] <- 55-8
seu$interval[seu$patient == 'P30'] <- 57-8
seu$interval[seu$patient == 'P33'] <- 81-14
seu$interval[seu$patient == 'P35'] <- 151-7
seu[['treatment']] <- 'aPD1+Chemo'
seu$modality <- 'Mono'
seu$prior <- 'No'
seu$cancertype <- 'NSCLC'
seu$celltype_major <- 'T_cells'
seu$time_point <- ifelse(seu$time_point == 'pre', 'Pre', 'Post')
seu$res_metric <- 'RECIST'
seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
matrix_tcr <- matrix_tcr[intersect(colnames(seu), rownames(matrix_tcr)),]
seu$clone <- matrix_tcr$clone.id[match(colnames(seu), rownames(matrix_tcr))]
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/NSCLC_Liu/seu_r1.qs')

# CRC_Li
# in terminal
# ls GSM* | awk -F '_' '{print $1}' | uniq| while read i; do mkdir $i; mv *$i*gz $i; done
# find -name "*matrix.mtx.gz" | while read i; do mv $i $(dirname $i)/matrix.mtx.gz; done
# find -name "*features.tsv.gz" | while read i; do mv $i $(dirname $i)/features.tsv.gz; done
# find -name "*barcodes.tsv.gz" | while read i; do mv $i $(dirname $i)/barcodes.tsv.gz; done
sample_list = list.files("/bigdata/zlin/Melanoma_meta/data/CRC_Li/", pattern = 'GSM')
input_dir_list = paste0(rep("/bigdata/zlin/Melanoma_meta/data/CRC_Li/", 14), sample_list)
seu_list <- lapply(sample_list, function(sample) {
  input_dir <- input_dir_list[str_detect(input_dir_list, sample)]
  count_matrix <- Read10X(input_dir)
  seu <- CreateSeuratObject(counts = count_matrix, min.cells=3, min.features=200)
  seu$SampleID <- sample
  return(seu)
})
names(seu_list) <- sample_list
qsave(seu_list, file = '/bigdata/zlin/Melanoma_meta/data/CRC_Li/raw.qs')

seu_list <- qread('/bigdata/zlin/Melanoma_meta/data/CRC_Li/raw.qs')
sample_list = list.files("/bigdata/zlin/Melanoma_meta/data/CRC_Li/", pattern = 'GSM')
seu_list <- lapply(seu_list, function(seu){
  print(unique(seu$SampleID))
  seu$percent_mito <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu$percent_ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu <- seu |> NormalizeData() |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress=c('nCount_RNA')) |>
    RunPCA() |>
    RunUMAP(dims = 1:20)
  # identify doublets
  # scDblFinder
  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce, dbr.sd=1)
  seu[['scDblFinder.class']] <- sce$scDblFinder.class
  # DoubletFinder
  sweep.list <- paramSweep_v3(seu, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn |>
    arrange(desc(BCmetric))
  pK <- pK[1, 'pK']
  pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
  nExp <- round(ncol(seu) * 0.08)
  seu <- doubletFinder_v3(seu, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:10)
  seu[['DoubletFinder.class']] <- seu@meta.data[, colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]
  seu[[colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]] <- NULL
  seu[[colnames(seu@meta.data)[grepl("pANN", colnames(seu@meta.data))]]] <- NULL
  seu[['singlet']] <- ifelse(seu$scDblFinder.class=='doublet' & seu$DoubletFinder.class=='Doublet', 'no', 'yes')
  # filter low-quality cells and doublets
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mito < 20 & singlet =='yes')
  seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
  seu <- seu%>%
    NormalizeData() |>
    FindVariableFeatures() |> 
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  # Preliminary annotation
  pred_bped_main <- SingleR(test = as.SingleCellExperiment(seu), ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
  seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
  pred_bped_fine <- SingleR(test = as.SingleCellExperiment(seu), ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
  seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
  # scGate
  seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
  return(seu)
})
qsave(seu_list, file = '/bigdata/zlin/Melanoma_meta/data/CRC_Li/list.qs')
seu_list <- qread('/bigdata/zlin/Melanoma_meta/data/CRC_Li/list.qs')
# check sample quality
(df <- data.frame(
  nfeature = sapply(seu_list, function(seu){nfeature <- sum(Matrix::rowSums(seu@assays$RNA$counts>0)>3)
  return(nfeature)}),
  ncell = sapply(seu_list, function(seu){ncell <- ncol(seu)
  return(ncell)})
))
meta_df <- data.frame(row.names = sample_list, 
                      patient = rep(c('P21','P24','P25','P27','P28','P30','P31'), each=2), 
                      time_point = rep(c('Post','Pre'),7),
                      treatment = rep('aPD1',14),
                      response = rep('pCR', 14))
meta_df$response[meta_df$patient == 'P31'] <- 'non-pCR'
meta_df$sample <- paste0(meta_df$patient, '_', meta_df$time_point)
# merge samples
for (i in 1:14){
  seu_list[[i]]$sample <- meta_df[i,'sample']
  seu_list[[i]]$response <- meta_df[i,'response']
  seu_list[[i]]$treatment <- meta_df[i,'treatment']
  seu_list[[i]]$patient <- meta_df[i,'patient']
  seu_list[[i]]$time_point <- meta_df[i,'time_point']
}
seu <- merge(x=seu_list[[1]], y=seu_list[2:length(seu_list)], add.cell.ids = meta_df$sample)
seu$dataset <- 'CRC_Li'
seu$patient <- paste0(seu$dataset,'_', seu$patient)
seu$sample <- paste0(seu$dataset,'_', seu$sample)
seu$prior <- 'No'
seu$modality <- 'Mono'
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('COL3A1','FAP', 'COL1A1','ATCA2'), 
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('KRT15','KRT17','EPCAM'),
                      c('MKI67')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Fibro','Endo','Mela','Epi','Proliferating')
DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'celltype_bped_main', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'scGate_multi', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' & seu$celltype_major == 'B-cells'] <- 'Plasma'
seu$celltype_major[seu$scGate_multi == 'Endothelial' & seu$celltype_major == 'unknown'] <- 'Endothelial cells'
seu$celltype_major[seu$scGate_multi == 'Epithelial' & seu$celltype_major == 'unknown'] <- 'Endothelial cells'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
table(seu$celltype_major)
seu <- subset(seu, subset = (celltype_major %in% c('unknown','Melanocytes','Neutrophils')) | (scGate_multi == 'Multi'), invert = T)
table(seu$celltype_major)
seu$cancertype <- 'CRC'
seu$interval <- 84
seu$response <- ifelse(seu$response == 'pCR', 'RE', 'NR')
seu$res_metric <- 'Pathology'
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/CRC_Li/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/CRC_Li/seu_r1.qs')

# IMCISION
matrix_count <- read.table('/bigdata/zlin/Melanoma_meta/data/HNSC_IMCISION/GSM7324294_Count_data_IMCISION.txt') |> as.sparse()
matrix_meta <- read.table('/bigdata/zlin/Melanoma_meta/data/HNSC_IMCISION/GSM7324295_Meta_data_IMCISION.txt', header = TRUE, sep = '\t')
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta)
seu$sample <- paste0(seu$patient, '_', seu$timepoint)
seu[["percent_mito"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent_ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
# Visualize QC metrics as a violin plot
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mito < 20 )
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3)
seu_list <- SplitObject(seu, split.by = "sample")
(df <- data.frame(
  nfeature = sapply(seu_list, function(seu_obj){nfeature <- sum(Matrix::rowSums(seu_obj@assays$RNA$counts>0)>3)
  return(nfeature)}),
  ncell = sapply(seu_list, function(seu_obj){ncell <- ncol(seu_obj)
  return(ncell)})
))
sample_rm <- rownames(filter(df, nfeature < 9000 | ncell < 200)); sample_rm
# [1] "Pat38_pre"  "Pat12_pre"  "Pat38_post" "Pat12_post" "Pat04_post"
rm(seu_list)
pt_rm <- str_split(sample_rm,'_', simplify = T)[,1]
seu <- subset(seu, subset = patient %in% pt_rm, invert=TRUE)
seu$cancertype <- 'HNSC'
names(seu@meta.data)[names(seu@meta.data) == 'timepoint'] <- 'time_point'
seu$time_point <- ifelse(seu$time_point == 'pre', 'Pre', 'Post')
seu$treatment <- 'aPD1+CTLA4'
seu$interval <- 28
seu$dataset <- 'HNSC_IMCISION'
seu$res_metric <- 'RECIST'
seu$prior <- 'No'
seu$modality <- 'Dual'
seu$patient <- paste0(seu$dataset, '_', seu$patient)
seu$sample <- paste0(seu$dataset, '_', seu$sample)
seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) |> 
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
sce <- as.SingleCellExperiment(seu)
pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac')
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$scGate_multi == 'panDC'] <- 'DC'
seu$celltype_major[seu$seurat_clusters == '15'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' & seu$celltype_major == 'B-cells'] <- 'Plasma'
seu$celltype_major[seu$scGate_multi == 'CD4T' & seu$celltype_major == 'unknown'] <- 'CD4+ T-cells'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
table(seu$celltype_major)
seu <- subset(seu, subset = (celltype_major %in% c('Epithelial cells','unknown','Endothelial cells','Neutrophils')) | (scGate_multi == 'Multi'), invert = T)
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/HNSC_IMCISION/seu_r1.qs')

#HNSC Luoma
files <- list.files('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/')
files_scrna_pre <- files[str_detect(files, 'pre-Tx_GEX_sc_tumor')]
pattern <- "P\\d+"
pt_rna <- str_extract(files_scrna_pre, pattern)
combined_pattern <- paste(paste0(pt_rna, '_post-Tx_GEX_sc_tumor'), collapse = "|")
files_scrna_post <- files[str_detect(files, combined_pattern)]
count_list <- purrr::map(paste0('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/' ,c(files_scrna_pre, files_scrna_post)), Read10X_h5)
clin_df <- readxl::read_xlsx('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/1-s2.0-S0092867422007231-mmc1.xlsx') |> 
  filter(`Pat. ID` %in% pt_rna) |> 
  tibble::column_to_rownames(var = 'Pat. ID')
clin_df_dup <- clin_df
rownames(clin_df_dup) <- paste0(rownames(clin_df), '_')
clin_df <- rbind(clin_df, clin_df_dup)
pt_rna_dup <- rep(pt_rna, times=2)
seu_list <- list()
timepoint <- rep(c('Pre','Post'), each=6)
treatment <- rep(c('Nivo+Ipi', 'Nivo', 'Nivo+Ipi', 'Nivo', 'Nivo', 'Nivo+Ipi'), times = 2)
for (i in 1:length(count_list)){
  seu <- CreateSeuratObject(counts = count_list[[i]], min.cells = 3, min.features = 200)
  seu[['time_point']] <- rep(c('Pre','Post'), each=6)[[i]]
  seu[['patient']] <- rep(pt_rna_dup, times=2)[[i]]
  seu[['sample']] <- paste0(seu$patient, '_', seu$time_point)
  seu[['treatment']] <- treatment[[i]]
  seu[['modality']] <- clin_df[i, 'Cohort']
  seu[['interval']] <- 28
  seu[['response']] <- clin_df[i, 'RECIST response excluding non measurable']
  print(unique(seu$sample))
  seu$percent_mito <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu$percent_ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu <- seu |> NormalizeData() |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress=c('nCount_RNA')) |>
    RunPCA() |>
    RunUMAP(dims = 1:20)
  # identify doublets
  # scDblFinder
  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce, dbr.sd=1)
  seu[['scDblFinder.class']] <- sce$scDblFinder.class
  # DoubletFinder
  sweep.list <- paramSweep_v5(seu, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn |>
    arrange(desc(BCmetric))
  pK <- pK[1, 'pK']
  pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
  nExp <- round(ncol(seu) * 0.08)
  seu <- doubletFinder_v5(seu, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:10)
  seu[['DoubletFinder.class']] <- seu@meta.data[, colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]
  seu[[colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]] <- NULL
  seu[[colnames(seu@meta.data)[grepl("pANN", colnames(seu@meta.data))]]] <- NULL
  seu[['singlet']] <- ifelse(seu$scDblFinder.class=='doublet' & seu$DoubletFinder.class=='Doublet', 'no', 'yes')
  # filter low-quality cells and doublets
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent_mito < 20 & singlet =='yes')
  seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
  seu <- seu%>%
    NormalizeData() |>
    FindVariableFeatures() |> 
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  # Preliminary annotation
  pred_bped_main <- SingleR(test = as.SingleCellExperiment(seu), ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
  seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
  pred_bped_fine <- SingleR(test = as.SingleCellExperiment(seu), ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
  seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
  # scGate
  seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
  seu <- RenameCells(object = seu, add.cell.id = paste0(paste0(seu$sample, '_', seu$patient)))
  seu_list[[i]] <- seu
}
seu_list <- lapply(seu_list, function(seu){
  return(CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data))
})

qsave(seu_list, file = '/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/list.qs')
seu_list <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/list.qs')

# # check sample quality
# df <- data.frame(
#   nfeature = sapply(seu_list, function(seu_obj){nfeature <- sum(Matrix::rowSums(seu_obj@assays$RNA@counts>0)>3)
#   return(nfeature)}),
#   ncell = sapply(seu_list, function(seu_obj){ncell <- ncol(seu_obj)
#   return(ncell)})
# )
# # remove samples with <200 cells or <10000 genes 
# sample_rm <- rownames(filter(df, nfeature < 10000 | ncell < 200)); sample_rm
# seu_list <- seu_list[-which(names(seu_list) %in% sample_rm)]
# merge samples
seu <- merge(x=seu_list[[1]], y=seu_list[2:length(seu_list)], add.cell.ids = names(seu_list))
seu <- JoinLayers(seu)
seu[['dataset']] <- 'HNSC_Luoma'
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
seu$treatment <- ifelse(seu$treatment == 'Nivo+Ipi', 'aPD1+CTLA4', 'aPD1')
seu$modality <- ifelse(seu$treatment == 'aPD1+CTLA4', 'Dual', 'Mono')
# files_tcr_pre <- files[str_detect(files, 'pre-Tx_TCR_sc_tumor')]
# pattern <- "P\\d+"
# pt_tcr <- str_extract(files_tcr_pre, pattern)
# combined_pattern <- paste(paste0(pt_tcr, '_post-Tx_TCR_sc_tumor'), collapse = "|")
# files_tcr_post <- files[str_detect(files, combined_pattern)]
# contig_list <- purrr::map(paste0('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/' ,c(files_tcr_pre, files_tcr_post)), read.csv)
# combined <- combineTCR(contig_list,
#                        samples = c(paste0(pt_tcr, '_pre'), paste0(pt_tcr, '_post')),
#                        ID = rep(pt_tcr,2), cells ="T-AB", removeNA = TRUE, removeMulti = TRUE)
# seu <- combineExpression(combined, seu,
#                             cloneCall="nt",
#                             group.by = "sample",
#                             proportion = FALSE,
#                             cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
seu$response <- 'NR'
seu$treatment <- ifelse(seu$treatment == 'Nivo+Ipi', 'aPD1+CTLA4', 'aPD1')
seu$cancertype <- 'HNSC'
seu$res_metric <- 'RECIST'
seu$prior <- 'No'
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) |> 
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),
                      c('MKI67')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Prolif')
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$scGate_multi == 'panDC'] <- 'DC'
seu$celltype_major[seu$seurat_clusters == '13'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell'] <- 'Plasma'
seu$celltype_major[seu$scGate_multi == 'CD4T' & seu$celltype_major == 'unknown'] <- 'CD4+ T-cells'
seu$celltype_major[seu$scGate_multi == 'CD8T' & seu$celltype_major == 'unknown'] <- 'CD8+ T-cells'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
table(seu$celltype_major)
seu <- subset(seu, subset = (celltype_major %in% c('Fibroblasts','Epithelial cells','unknown','Endothelial cells','Neutrophils','Myocytes')) | (scGate_multi == 'Multi'), invert = T)
seu$patient <- paste0(seu$dataset, '_', seu$patient)
seu$sample <- paste0(seu$dataset, '_', seu$sample)
seu$res_metric[seu$patient == 'HNSC_Luoma_P18'] <- 'Pathology'
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/seu_r1.qs')

# PCa_Hawley
seu <- readRDS('/bigdata/zlin/Melanoma_meta/data/PCa_Hawley/primecut.integrated.geneexp.seurat.rds')
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
seu <- subset(seu, subset = patient %in% c('Patient1', 'Patient7'))
seu$dataset <- 'PCa_Hawley'
seu$cancertype <- 'PCa'
seu$res_metric <- 'PSA test'
seu$response <- 'NR'
names(seu@meta.data)[names(seu@meta.data) == 'treatment'] <- 'time_point'
seu$time_point <- as.character(seu$time_point)
seu$time_point <- ifelse(seu$time_point == 'Pre-Treatment', 'Pre', 'Post')
seu$interval <- 70
seu$prior <- 'Yes'
seu$patient <- paste0(seu$dataset, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$treatment <- 'aPD1'
seu$modality <- 'Mono'
seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu |> 
  NormalizeData() |> 
  FindVariableFeatures()|>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) |> 
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
sce <- as.SingleCellExperiment(seu)
pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('COL3A1','FAP', 'COL1A1','ATCA2'), 
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('KRT15','KRT17','EPCAM'),
                      c('MKI67')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Fibro','Endo','Mela','Epi','Proliferating')
DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'celltype_bped_main', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'scGate_multi', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '19'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' & seu$celltype_major == 'B-cells'] <- 'Plasma'
seu$celltype_major[seu$scGate_multi == 'CD4T' & seu$celltype_major == 'unknown'] <- 'CD4+ T-cells'
seu$celltype_major[seu$scGate_multi == 'CD8T' & seu$celltype_major == 'unknown'] <- 'CD8+ T-cells'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$seurat_clusters)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = (celltype_major %in% c('unknown','Neutrophils')) | (scGate_multi == 'Multi'), invert = T)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/PCa_Hawley/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/PCa_Hawley/seu_r1.qs')

# TNBC_Shiao 
# CD45+
sce <- dior::read_h5(file = "/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/adata_cd45n.h5", target.object = "singlecellexperiment")
metadata <- read.csv(file = "/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/obs_cd45n.csv", row.names = 1) |> select(batch, cohort, treatment, response_group)
seu <- CreateSeuratObject(counts = assay(sce, "X"), meta.data = metadata)
qsave(seu, file = "/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/seu_cd45n.qs")
# CD45-
sce <- dior::read_h5(file = "/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/adata_cd45p.h5", target.object = "singlecellexperiment")
metadata <- read.csv(file = "/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/obs_cd45p.csv", row.names = 1) |> select(batch, CD45_enrich, batch_num, cohort, pCR, cleared_nodes, treatment, patient_treatment, combined_tcr)
seu <- CreateSeuratObject(counts = assay(sce, "X"), meta.data = metadata)
qsave(seu, file = "/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/seu_cd45p.qs")
# 
seu <- qread("/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/seu_cd45p.qs")
seu$response <- ifelse(seu$pCR == 'NR', 'NR', 'RE')
seu$patient <- seu$cohort
seu$dataset <- 'TNBC_Shiao'
seu$prior <- 'No'
seu$time_point <- ifelse(seu$treatment == 'Base', 'Pre', ifelse(seu$treatment == 'PD1', 'Post', 'Post2'))
seu$treatment <- 'aPD1'
seu$modality <- 'Mono'
seu$res_metric <- 'Pathology'
seu$interval <- 20
seu$cancertype <- 'TNBC'
seu$sample <- paste0(seu$patient, "_", seu$time_point)
seu_list <- SplitObject(seu, split.by = 'sample')
# check sample quality
df <- data.frame(
  nfeature = sapply(seu_list, function(seu_obj){nfeature <- sum(Matrix::rowSums(seu_obj@assays$RNA$counts>0)>3)
  return(nfeature)}),
  ncell = sapply(seu_list, function(seu_obj){ncell <- ncol(seu_obj)
  return(ncell)})
)
rm(seu_list)
# remove samples with <200 cells or <9000 genes 
sample_rm <- rownames(filter(df, nfeature < 9000 | ncell < 200)); sample_rm
pt_rm <- str_split(sample_rm[!str_detect(sample_rm, '_Post2')], '_', simplify = T)[,1]
seu <- subset(seu, subset = patient %in% pt_rm | time_point == 'Post2', invert = T)
seu <- STACAS::StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
seu <- seu%>%
  NormalizeData() |>
  FindVariableFeatures() |>
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) |>
  ScaleData() |> 
  RunPCA(verbose=FALSE) |> 
  RunHarmony(group.by.vars = "sample") |> 
  RunUMAP(reduction = "harmony", dims = 1:20) |>  
  FindNeighbors(reduction = "harmony", dims = 1:20) |> FindClusters()
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes)
sce <- as.SingleCellExperiment(seu)
pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
qsave(seu, '/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/seu_r1_processing.qs')
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38','JCHAIN'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4','CXCR3'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      # c('COL3A1','FAP', 'COL1A1','ATCA2'), 
                      # c('PECAM1','VWF', 'ENG'), 
                      # c('MLANA','MITF', 'TYR'), 
                      # c('KRT15','KRT17','EPCAM'),
                      c('MKI67')
)
# names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Fibro','Endo','Mela','Epi','Proliferating')
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Proliferating')
DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters))))
DotPlot(seu, group.by = 'celltype_bped_main', features = genes_to_check) + RotatedAxis()
DotPlot(seu, group.by = 'scGate_multi', features = genes_to_check) + RotatedAxis()
table(seu$celltype_bped_main, seu$scGate_multi)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$scGate_multi == 'panDC'] <- 'DC'
seu$celltype_major[seu$scGate_multi == 'Mast' | seu$seurat_clusters == '17'] <- 'Mast'
seu$celltype_major[seu$scGate_multi == 'PlasmaCell' ] <- 'Plasma'
seu$celltype_major[seu$scGate_multi == 'CD4T' & seu$celltype_major == 'unknown'] <- 'CD4+ T-cells'
seu$celltype_major[seu$scGate_multi == 'CD8T' & seu$celltype_major == 'unknown'] <- 'CD8+ T-cells'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major)))) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = (celltype_major %in% c('unknown', 'Neutrophils', 'Endothelial cells', 'Fibroblasts', 'Epithelial cells', 'Melanocytes')) | (scGate_multi == 'Multi'), invert = T)
pt_rm <- as.data.frame.matrix(table(seu$patient, seu$time_point)) |> filter(Post == 0 | Pre == 0) |> rownames()
seu <- subset(seu, subset = patient %in% pt_rm, invert = T)
seu <- subset(seu, subset = patient %in% c('Patient03T1', 'Patient03T2'), invert = T)
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
qsave(seu, file = '/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/seu_r1.qs')
seu <- qread('/bigdata/zlin/Melanoma_meta/data/TNBC_Shiao/seu_r1.qs')

# ESCC_Ji
dir <- '/bigdata/zlin/Melanoma_meta/data/OMIX005710/'
files <- list.files(dir)[-1]
list_seu <- lapply(c(paste0(dir, files)), function(file){
  print(file)
  count_matrix <- Read10X(file)
  seu <- CreateSeuratObject(counts = count_matrix, min.features = 200)
  seu$percent_mito <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = percent_mito<=50 & nFeature_RNA < quantile(seu$nFeature_RNA, probs = 0.98) & nCount_RNA < quantile(seu$nCount_RNA, probs = 0.98))
  seu <- seu[rowSums(seu@assays$RNA$counts > 0) >= 5,]
  return(seu)
})

qc <- readxl::read_excel('/bigdata/zlin/Melanoma_meta/data/OMIX005710/13073_2024_1320_MOESM2_ESM.xlsx')
qc$pt <- str_split(qc$`Sample ID`, '_', simplify = T)[,1]
qc$tissue <- str_split(qc$`Sample ID`, '_', simplify = T)[,2]
qc$time_point <- str_split(qc$`Sample ID`, '_', simplify = T)[,3]
qc <- qc |> 
  filter(tissue == 'T') |> 
  group_by(pt) |> 
  mutate(n = n()) |> 
  filter(n == 2) |> 
  arrange(desc(`Cell count`))

list_summary <- lapply(list_seu, function(seu){
  qc_sample <- c()
  qc_sample[1] <- dim(seu)[2]
  qc_sample[2] <- sum(seu$nCount_RNA)
  qc_sample[3] <- median(seu$nCount_RNA)
  qc_sample[4] <- median(seu$nFeature_RNA)
  return(qc_sample)
})

df <- do.call(rbind, list_summary) |> data.frame() |> arrange(desc(X1))
df$sample <- 'unknown'
df$sample[2] <- 'P19_T_B'
df$sample[13] <- 'P19_T_A'
df$sample[33] <- 'P6_T_A'
df$sample[20] <- 'P18_T_A'
df$sample[24] <- 'P6_T_B'
df$sample[6] <- 'P18_T_B'
df$sample[29] <- 'P16_T_A'
df$sample[] <- 'P16_T_B'
df$sample[1] <- 'P21_T_B'



