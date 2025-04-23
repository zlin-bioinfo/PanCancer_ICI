library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(janitor)
# MISTy
library(mistyR)
library(future)
# Distances
library(distances)
path <- 'data/SKCM_Plozniak/Melanoma2_1/outs/'
imaged <- Read10X_Image(image.dir = 'data/SKCM_Plozniak/Melanoma2_1/outs/spatial/',
                        image.name = 'tissue_lowres_image.png',
                        slice = 'Melanoma2_1')
seu <- Load10X_Spatial(path,
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial",
                       slice = "Melanoma2_1", image = imaged)
geometry <- GetTissueCoordinates(seu, cols = c("imagerow", "imagecol"), scale = NULL)
composition <- read.csv('scripts/Spatial/adata_vis.csv', row.names = 'X', check.names = F)
composition <- composition[,-c(1:6)]
seu@meta.data <- cbind(seu@meta.data, composition)
# Calculating the radius
geometry <- geometry[,-3]
geom_dist <- as.matrix(distances(geometry))  
dist_nn <- apply(geom_dist, 1, function(x) (sort(x)[2]))
paraview_radius <- ceiling(mean(dist_nn+ sd(dist_nn)))

# Create views
tumor_views <- create_initial_view(composition) %>%
  add_paraview(geometry, l= paraview_radius, family = "gaussian")

# Run misty and collect results
run_misty(tumor_views, "/home/zlin/workspace/PanCancer_ICI/scripts/Spatial/result/vignette_structural_pipeline")
misty_results <- collect_results("/home/zlin/workspace/PanCancer_ICI/scripts/Spatial/result/vignette_structural_pipeline")
misty_results %>%
  plot_improvement_stats("multi.R2") %>% 
  plot_improvement_stats("gain.R2")
misty_results %>% plot_interaction_heatmap(view = "intra", clean = TRUE)
misty_results$importances.aggregated %>%
  filter(view == "intra", Predictor == "Malignant.CNA..") %>%
  arrange(-Importance)
SpatialFeaturePlot(seu, keep.scale = NULL, features = c("Malignant.CNA..","Macro_SPP1"), image.alpha = 0)
misty_results %>% plot_interaction_heatmap(view = "para.439", clean = TRUE, 
                                           trim = 1.75, trim.measure = "gain.R2",
                                           cutoff = 0.5) 


plotseuplot1 <- VlnPlot(seu, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(seu, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
seu <- SCTransform(seu, assay = "Spatial", verbose = FALSE) |>
  RunPCA(assay = "SCT", verbose = FALSE) |>
  FindNeighbors(reduction = "pca", dims = 1:30) |>
  FindClusters(verbose = FALSE) |>
  RunUMAP(reduction = "pca", dims = 1:30)

p1 <- DimPlot(seu, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(seu, label = TRUE, label.size = 3)
p1 + p2

genes_to_check = list(c('CD3D','CD4','CD8A'), # T cells 'CD8B'
                      c('KLRD1','FCGR3A'), 
                      c('MKI67', 'TOP2A'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('MZB1','JCHAIN'),
                      c('KIT','TPSAB1'),
                      c('LILRA4','PLD4'),
                      c('CLEC9A','CD1C','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('PECAM1','VWF', 'ENG'),
                      c("RGS5",'ACTA2'),
                      c('COL1A1','FAP'),
                      c('KRT19', 'EPCAM'),
                      c('MLANA','TYR'),
                      c('RPL11','RPL10A')
)
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  RotatedAxis()
cd4tnaive_features <- list(c(
  'TCF7', 'CD4', 'CCR7', 'IL7R', 'FHIT', 'LEF1', 'MAL', 'NOSIP', 'LDHB', 'PIK3IP1'
))
seu <- AddModuleScore(
  object = seu,
  features = cd4tnaive_features,
  ctrl = 5,
  name = 'CD4_T-naive'
)
SpatialFeaturePlot(seu, features = c('CD4_T-naive1','C1QC','FOXP3'), min.cutoff = c(0.3,1,0.5), alpha = c(0.2, 1))

imaged <- Read10X_Image(image.dir = 'data/GSE175540/GSE175540_RAW/',
                        image.name = 'GSM5924030_ffpe_c_2_tissue_lowres_image.png.gz',
                        slice = 'ffpe_c_2')







