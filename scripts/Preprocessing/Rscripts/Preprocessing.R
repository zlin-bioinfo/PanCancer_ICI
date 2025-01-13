rm(list=ls())
pkgs <- c('Seurat','DoubletFinder','scDblFinder','tidyr','plyr','dplyr','stringr','SingleR','scGate','ggsci','BiocParallel','harmony','RColorBrewer','SingleCellExperiment','STACAS','janitor','qs','qs2','COSG','fgsea','msigdbr','infercnv','tidyplots','ggplot2')
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

