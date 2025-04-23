
# 
# library(scRNAseq)
# sce.nest <- NestorowaHSCData()
# library(scater)
# stats <- perCellQCMetrics(sce.nest)
# qc <- quickPerCellQC(stats, percent_subsets="altexps_ERCC_percent")
# sce.nest <- sce.nest[,!qc$discard]
# library(scran)
# set.seed(101000110)
# clusters <- quickCluster(sce.nest)
# sce.nest <- computeSumFactors(sce.nest, clusters=clusters)
# sce.nest <- logNormCounts(sce.nest)
# set.seed(00010101)
# dec.nest <- modelGeneVarWithSpikes(sce.nest, "ERCC")
# top.nest <- getTopHVGs(dec.nest, prop=0.1)
# set.seed(101010011)
# sce.nest <- denoisePCA(sce.nest, technical=dec.nest, subset.row=top.nest)
# sce.nest <- runTSNE(sce.nest, dimred="PCA")
# snn.gr <- buildSNNGraph(sce.nest, use.dimred="PCA")
# colLabels(sce.nest) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
# library(slingshot)
# sce.sling <- slingshot(sce.nest, reducedDim='PCA', clusterLabels = colLabels(sce.sling))
# 
# embedded <- embedCurves(sce.sling, "TSNE")
# embedded <- slingCurves(embedded)[[1]] # only 1 path.
# embedded <- data.frame(embedded$s[embedded$ord,])
# 
# pseudo.paths <- slingPseudotime(sce.sling)
# shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
# 
# 
# gg <- plotTSNE(sce.sling, colour_by=I(shared.pseudo)) 
# seu$shared.pseudo <- shared.pseudo
# gg <- DimPlot(seu, )
# embedded <- embedCurves(sce.sling, "TSNE")
# embedded <- slingCurves(embedded)
# for (path in embedded) {
#   embedded <- data.frame(path$s[path$ord,])
#   gg <- gg + geom_path(data=embedded, aes(x=TSNE1, y=TSNE2), size=1.2)
# }
# gg
# library(Seurat)
# seu <- sce.nest |> as.Seurat(counts = 'data')
# 
# sce <- as.SingleCellExperiment(seu)
# DimPlot(seu, group.by = 'celltype_r2')
# sce.sling <- slingshot(sce, reducedDim='PCA', sce$seurat_clusters)
# embedded <- embedCurves(sce.sling, "UMAP")
# pseudo.paths <- slingPseudotime(sce.sling)
# shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
# 
# gg <- DimPlot(seu)
# seu$shared.pseudo <- shared.pseudo
# embedded <- embedCurves(sce.sling, "UMAP")
# embedded <- slingCurves(embedded)
# for (path in embedded) {
#   embedded <- data.frame(path$s[path$ord,])
#   gg <- gg + geom_path(data=embedded, aes(x=umap_1, y=umap_2), size=1.2)
# }
# gg

library(Scissor)
# setwd("/home/zlin/workspace/PanCancer_ICI")
Scissor_m <- function(bulk_dataset, sc_dataset, phenotype, tag = NULL, alpha = NULL, 
                      cutoff = 0.2, family = c("gaussian", "binomial", "cox"), 
                      Save_file = "Scissor_inputs.RData", Load_file = NULL) 
{
  library(Seurat)
  library(Matrix)
  library(preprocessCore)
  if (is.null(Load_file)) {
    common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
    if (length(common) == 0) {
      stop("There is no common genes between the given single-cell and bulk samples.")
    }
    if (class(sc_dataset) == "Seurat") {
      sc_exprs <- as.matrix(sc_dataset@assays$RNA$data)
      network <- as.matrix(sc_dataset@graphs$RNA_snn)
    }
    else {
      sc_exprs <- as.matrix(sc_dataset)
      Seurat_tmp <- CreateSeuratObject(sc_dataset)
      Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", 
                                         verbose = F)
      Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
      Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), 
                           verbose = F)
      Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, 
                                  verbose = F)
      network <- as.matrix(Seurat_tmp@graphs$RNA_snn)
    }
    diag(network) <- 0
    network[which(network != 0)] <- 1
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, 
    ])
    dataset0 <- as.matrix(dataset0)
    dataset1 <- normalize.quantiles(dataset0)
    rownames(dataset1) <- rownames(dataset0)
    colnames(dataset1) <- colnames(dataset0)
    Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
    Expression_cell <- dataset1[, (ncol(bulk_dataset) + 1):ncol(dataset1)]
    X <- cor(Expression_bulk, Expression_cell)
    quality_check <- quantile(X)
    print("|**************************************************|")
    print("Performing quality-check for the correlations")
    print("The five-number summary of correlations:")
    print(quality_check)
    print("|**************************************************|")
    if (quality_check[3] < 0.01) {
      warning("The median correlation between the single-cell and bulk samples is relatively low.")
    }
    if (family == "binomial") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        print(sprintf("Current phenotype contains %d %s and %d %s samples.", 
                      z[1], tag[1], z[2], tag[2]))
        print("Perform logistic regression on the given phenotypes:")
      }
    }
    if (family == "gaussian") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        tmp <- paste(z, tag)
        print(paste0("Current phenotype contains ", paste(tmp[1:(length(z) - 
                                                                   1)], collapse = ", "), ", and ", tmp[length(z)], 
                     " samples."))
        print("Perform linear regression on the given phenotypes:")
      }
    }
    if (family == "cox") {
      Y <- as.matrix(phenotype)
      if (ncol(Y) != 2) {
        stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
      }
      else {
        print("Perform cox regression on the given clinical outcomes:")
      }
    }
    save(X, Y, network, Expression_bulk, Expression_cell, 
         file = Save_file)
  }
  else {
    load(Load_file)
  }
  if (is.null(alpha)) {
    alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
               0.6, 0.7, 0.8, 0.9)
  }
  for (i in 1:length(alpha)) {
    set.seed(123)
    fit0 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, nlambda = 100, 
                  nfolds = min(10, nrow(X)))
    fit1 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
    if (family == "binomial") {
      Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
    }
    else {
      Coefs <- as.numeric(fit1$Beta)
    }
    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    percentage <- (length(Cell1) + length(Cell2))/ncol(X)
    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", 
                  length(Cell1), length(Cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", 
                  formatC(percentage * 100, format = "f", digits = 3)))
    if (percentage < cutoff) {
      break
    }
    cat("\n")
  }
  print("|**************************************************|")
  return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, 
                          family = family), Coefs = Coefs, Scissor_pos = Cell1, 
              Scissor_neg = Cell2))
}



# Main
datasets <- c('SKCM_Becker', 'SKCM_Plozniak', 'BCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao', 
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'NSCLC_Yan',
              'CRC_Li', 'CRC_Chen', 'PCa_Hawley','HCC_Guo','RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2'))
  seu <- seu[,!is.na(seu$celltype_r2)]
  seu <- subset(seu, subset = celltype_r2 %in% mye)
  seu$cell.id <- colnames(seu)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu <- seu |>
  subset(subset = sample %in% filter_sample) |>
  NormalizeData() 

library(ggplot2)
library(ggridges)
df <- read.csv('tables/cd4t_paga.csv')
density_data <- df |>
  group_by(celltype_r2) |>
  summarise(peak_density = median(density(dpt_pseudotime, from = 0, to = 0.1)$y)) |>
  arrange(desc(peak_density))
df |> 
  filter(time_point %in% c('Pre','On')) |> 
  mutate(
         celltype_r2 = factor(celltype_r2, levels = density_data$celltype_r2)) |> 
  ggplot(aes(x=dpt_pseudotime, fill=time_point, y=celltype_r2)) +
  geom_density_ridges(alpha=0.5) + xlim(0,0.1) + theme_ridges()


library(Seurat)
seu <- readRDS('/home/zlin/workspace/PanCancer_ICI/data/HCC_Liu/P10T_Spatial.rds.gz')



