#!/usr/bin/env Rscript
rm(list=ls())
pkgs <- c('Seurat','DoubletFinder','scDblFinder','tidyr','plyr','dplyr','stringr','SingleR','scGate','ggsci','qs','BiocParallel','harmony','RColorBrewer','scCustomize','SingleCellExperiment','STACAS','janitor')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")
setwd('/home/zlin/workspace/PanCancer_ICI/')
# scGate model
scGate_models_DB <- get_scGateDB()
# Loading reference
bped <- celldex::BlueprintEncodeData()
bped <- bped[, bped$label.main %in% c("CD4+ T-cells", "CD8+ T-cells", "NK cells", "B-cells", "DC", "Monocytes", "Macrophages", "Neutrophils", "Endothelial cells", "Epithelial cells", "Myocytes", "Fibroblasts", "Melanocytes")]
# cellcycle genelist
exp.mat <- read.table(file = "./tables/nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, 
                      as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes 
data(EnsemblGeneTable.Hs)
# # import cell type signature
# sign_celltype <- read.csv('/bigdata/zlin/PanCancer_ICI/tables/celltype_signature.csv') |> as.list() |> lapply(function(x) x[x != ''])
# # noise genes
# gene_rm <- readxl::read_excel('/bigdata/zlin/PanCancer_ICI/tables/41586_2022_5400_MOESM3_ESM.xlsx', sheet = '1d_Genes excluded', skip = 2)[,-1]
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
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, ct]
      }
    }
    return(seu)
  }
}

# SKCM_Plozniak
clin_info <- readxl::read_xlsx('./data/SKCM_Pozniak/1-s2.0-S0092867423013223-mmc1.xlsx') |> 
  data.frame() |> 
  select(-Abbreviations)
pt_included <- unique(clin_info$Patient.ID) |> setdiff(c('16','25','29','39'))
seu <- readRDS('./data/SKCM_Pozniak/Entire_TME.rds')
colnames(seu@meta.data)[3] <- 'patient'
seu <- subset(seu, subset = patient %in% pt_included)
seu$dataset <- 'SKCM_Plozniak'
colnames(seu@meta.data)[2] <- 'time_point'
seu$time_point <- ifelse(seu$time_point == 'BT', 'Pre', 'On')
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$response <- clin_info$Responder.[match(seu$orig.ident, clin_info$Sample.ID)]
seu$response <- ifelse(seu$response == 1, 'RE', 'NR')
seu$treatment <- clin_info$Treatment[match(seu$orig.ident, clin_info$Sample.ID)]
seu$treatment <- ifelse(seu$treatment %in% c('Nivolumab','Pembrolizumab'), 'aPD1', 'aPD1+CTLA4')
seu$modality <- ifelse(seu$treatment == 'aPD1+CTLA4', 'Dual', 'Mono')
seu$res_metric <- clin_info$Criterion
seu$res_metric <- ifelse(str_detect(seu$res_metric, 'pCR'), 'Pathology', 'RECIST')
seu$cancertype <- 'SKCM'
seu$prior <- 'No'
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)

seu_list <- SplitObject(seu, split.by = 'sample')
seu_list <- lapply(seu_list, function(seu){
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
  seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
  seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
  seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
  return(seu)
})
qsave(seu_list, './data/SKCM_Pozniak/seu_list.qs')
df <- data.frame(
  nfeature = sapply(seu_list, function(seu_obj){nfeature <- sum(Matrix::rowSums(seu_obj@assays$RNA$counts>0)>3)
  return(nfeature)}),
  ncell = sapply(seu_list, function(seu_obj){ncell <- ncol(seu_obj)
  return(ncell)})
)
print(df)









