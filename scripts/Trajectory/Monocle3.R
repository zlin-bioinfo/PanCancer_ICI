library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
seu <- qs_read('data/BRCA_Bassez1/seu_r2.qs2')
seu <- subset(seu, subset = celltype_main == "CD4+T")
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$patient)
seu <- seu |> 
  NormalizeData() |>
  FindVariableFeatures()  |>
  ScaleData() |>
  RunPCA(verbose=FALSE) |> 
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony') |> 
  FindNeighbors(reduction = "harmony", dims = 1:20) |>
  FindClusters() |> 
  RunUMAP(dims = 1:10, reduction = 'harmony') |> 
  JoinLayers()
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap', shuffle = T)
cds <- as.cell_data_set(seu)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
seu <- subset(seu, subset = seurat_clusters %in% c(14,16), invert=T)
cds <- as.cell_data_set(seu)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = T, label_leaves = FALSE, label_branch_points = FALSE, color_cells_by='celltype_r2')
cds <- order_cells(cds, root_cells = colnames(cds)[cds$celltype_r2 == 'CD4_T-naive'])
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)


