pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','COSG','BPCells','SeuratExtend','MetBrewer','ggplot2','slingshot','CytoTRACE2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))

seu <- qs_read('data/SKCM_Becker/seu_r2.qs2')
seu <- subset(seu, subset = celltype_main == "CD4+T")
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
seu <- seu |> 
  NormalizeData() |>
  FindVariableFeatures()  |>
  ScaleData() |>
  RunPCA(verbose=FALSE) |> 
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony') |> 
  FindNeighbors(reduction = "harmony", dims = 1:20) |>
  FindClusters() |> 
  RunUMAP(dims = 1:20, reduction = 'harmony') |> 
  JoinLayers()
DimPlot(seu, group.by = 'celltype_r2', reduction = 'harmony')
sce <- as.SingleCellExperiment(seu)

sce.slingshot <- slingshot(sce, 
                           clusterLabels = sce$celltype_r2, 
                           reducedDim = reducedDims(sce)$HARMONY, 
                           start.clus = "CD4_T-naive",
                           # end.clus = c("CD4_Tctl", "CD4_Treg", "CD4_Tfh"),
                           approx_points = 150)
embedded <- embedCurves(sce.slingshot, "UMAP")
embedded <- slingCurves(embedded)[[1]] # only 1 path.
embedded <- data.frame(embedded$s[embedded$ord,])

DimPlot(seu, group.by = 'celltype_r2', reduction = 'harmony') +
geom_path(data = as.data.frame(slingCurves(embedded)[[1]]$s[slingCurves(embedded)[[1]]$ord, ]), aes(umap_1, umap_2), size = 0.5, col = "black") +
  geom_path(data = as.data.frame(slingCurves(embedded)[[2]]$s[slingCurves(embedded)[[2]]$ord, ]), aes(umap_1, umap_2), size = 0.5, col = "black") +
  geom_path(data = as.data.frame(slingCurves(embedded)[[3]]$s[slingCurves(embedded)[[3]]$ord, ]), aes(umap_1, umap_2), size = 0.5, col = "black")

pseudo.paths <- slingPseudotime(sce.slingshot)
head(pseudo.paths)

seu2 <- qs_read('data/cytotrace2.qs2')
seu2@meta.data |> head()

dim(seu2)
seu$CytoTRACE2_Score <- seu2$CytoTRACE2_Score
seu$
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap') 
FeaturePlot(seu, features = 'CytoTRACE2_Score')
seu$slingshot <- pseudo.paths[,1]
seu@meta.data |> 
  ggplot(aes(celltype_r2, CytoTRACE2_Score)) +
  geom_violin() +
  geom_boxplot() + facet_wrap(.~ time_point, ncol = 1)

plots <- plotData(seu2, pc_dims = 15, seed = 42, is_seurat = TRUE)



