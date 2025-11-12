pkgs <- c('Seurat','DoubletFinder','scDblFinder','tidyr','plyr','dplyr','stringr','SingleR','scGate','ggsci','BiocParallel','harmony','RColorBrewer','SingleCellExperiment','STACAS','janitor','qs','qs2','COSG','fgsea','msigdbr','infercnv','tidyplots','ggplot2','infercnv')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

# scGate model
scGate_models_DB <- get_scGateDB()
immune_cells_scgate <- c('Bcell', 'CD4T', 'CD8T', 'NK', 'Monocyte', 'Macrophage', 'PlasmaCell','panDC','Neutrophils','Mast')
non_immune_cells_scgate <- setdiff(names(scGate_models_DB$human$TME_HiRes), c(immune_cells_scgate, 'unknown', 'Multi'))
# Loading reference
bped <- celldex::BlueprintEncodeData()
bped <- bped[, bped$label.main %in% c("CD4+ T-cells", "CD8+ T-cells", "NK cells", "B-cells", "DC", "Monocytes", "Macrophages", "Neutrophils", "Endothelial cells", "Epithelial cells", "Myocytes", "Fibroblasts", "Melanocytes")]
immune_cells_bped <- c('B-cells', 'CD4+ T-cells', 'CD8+ T-cells', 'DC', 'Monocytes', 'Macrophages', 'Neutrophils', 'NK cells')
non_immune_cells_bped <- setdiff(unique(bped$label.main), c(immune_cells_bped, 'unknown'))
# cellcycle genelist
exp.mat <- read.table(file = "./tables/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, 
                      as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes 
data(EnsemblGeneTable.Hs)
preprocessing <- function(seu, sorted = F, .ncores = 10, .resolution = 0.5){
  seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data, min.cells=5, min.features=400)
  seu$percent.mito <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu$percent.ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable = EnsemblGeneTable.Hs)
  # cell cycle
  seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  # automatic annotation
  seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = .ncores)
  sce <- as.SingleCellExperiment(seu)
  pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(.ncores))
  seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
  pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(.ncores))
  seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
  seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
  seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
  seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
  # filtering
  if (sorted == F){
    seu <- seu[,! ((seu$scGate_multi %in% non_immune_cells_scgate) & (seu$celltype_bped_main %in% immune_cells_bped) |
                     (seu$scGate_multi %in% immune_cells_scgate) & (seu$celltype_bped_main %in% non_immune_cells_bped) |
                     seu$scGate_multi == 'Multi')]
  }
  # integration
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
  seu <- seu |> 
    NormalizeData() |>
    FindVariableFeatures(nfeatures = 3000)  |>
    ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) |>
    RunPCA(verbose=FALSE) |>
    IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                    new.reduction = 'harmony') |> 
    FindNeighbors(reduction = "harmony", dims = 1:20) |>
    FindClusters(resolution = .resolution) |> 
    RunUMAP(dims = 1:20, reduction = 'harmony') |> 
    JoinLayers()
  return(seu)
}

paramSweep_v5 <- function (seu, PCs = 1:10, sct = FALSE, num.cores = 1) 
{
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
    output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct, mc.cores = num.cores)
    stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep, 
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
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, ct]
      }
    }
    return(seu)
  }
}

samples <- list.files("./data/SKCM_Essen/cellranger")
samples <- samples[!samples %in% c('E19324','E23347','E23368')] # remove poor quality samples
seu_list <- lapply(samples, function(sample){
  print(sample)
  if (sample == 'E23346'){
    count_matrix <- Read10X(paste0("./data/SKCM_Essen/cellranger/", sample,"/outs/raw_feature_bc_matrix"))
    # count_matrix <- Matrix(as.matrix(count_matrix),sparse=TRUE)
  } else {
    count_matrix <- Read10X(paste0("./data/SKCM_Essen/cellranger/", sample,"/raw_feature_bc_matrix"))
  }
  seu <- CreateSeuratObject(counts = count_matrix, min.cells = 5, min.features = 400)
  if (ncol(seu) < 300) {
    print(paste("Sample", sample, "has fewer than", 300, "cells. Skipping."))
    return(NULL)  # Skip this samples
  }
  seu$sample <- sample
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
  seu <- subset(seu, subset = nFeature_RNA < 8000 & percent_mito < 20 & singlet =='yes')
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
print('Preprocessing by sample done!')
names(seu_list) <- samples
qs_save(seu_list, file = './data/SKCM_Essen/list.qs2')

seu_list <- qs_read('./data/SKCM_Essen/list.qs2')
seu_list <- Filter(Negate(is.null), seu_list)
seu <- merge(x=seu_list[[1]], y=seu_list[2:length(seu_list)])
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
seu <- seu |> 
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 3000)  |>
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) |>
  RunPCA(verbose=FALSE) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony') |> 
  FindNeighbors(reduction = "harmony", dims = 1:20) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:20, reduction = 'harmony') |> 
  JoinLayers()
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
qs_save(seu, file = './data/SKCM_Essen/processing.qs2')

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
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters))), label = T) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '19'] <- 'pDC'
seu$celltype_major[seu$seurat_clusters == '21'] <- 'Myocytes'
seu$celltype_major[seu$seurat_clusters == '17'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '20'] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters == '14'] <- 'Cycling T/NK'
seu$celltype_major[seu$seurat_clusters == '22'] <- 'Doublets'
seu$celltype_major[seu$seurat_clusters == '7'] <- 'Epithelial cells'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('CD4T','CD8T','Endothelial','Macrophage','panDC','Bcell','NK'), 
                                to = c('CD4+ T-cells','CD8+ T-cells','Endothelial cells','Macrophages','DC','B-cells','NK cells'))
marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('Epithelial cells','Epithelial','panDC','unknown','Monocyte','Myocytes','Doublets','Fibroblast','Mast'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

clin_info <- readxl::read_xlsx('./data/SKCM_Essen/clin_info.xlsx')
seu$patient <- clin_info$patient[match(seu$sample, clin_info$sample)]
seu$time_point <- clin_info$time_point[match(seu$sample, clin_info$sample)]
seu$matched <- clin_info$matched[match(seu$sample, clin_info$sample)]
seu$interval <- clin_info$interval[match(seu$sample, clin_info$sample)]
seu$response <- clin_info$response[match(seu$sample, clin_info$sample)]
seu$modality <- clin_info$treatment[match(seu$sample, clin_info$sample)]
seu$treatment <- ifelse(seu$modality == 'Mono', 'aPD1', 'aPD1+CTLA4')
seu <- subset(seu, subset = sample == 'E23359', invert = T) # IL2 treatment
seu$sample_id <- seu$sample
seu$cohort <- 'SKCM_this study'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$res_metric <- 'RECIST'
seu$prior <- 'No'
qs_save(seu, file = './data/SKCM_Essen/seu_r1.qs2')

seu <- qs_read('./data/SKCM_Essen/seu_r1.qs2')