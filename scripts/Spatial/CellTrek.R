options(stringsAsFactors = F)
library("qs2")
library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")
library(BPCells)
options(future.globals.maxSize = 8 * 1024^3)  # 8 GiB limit

matrix_dir = 'filtered_feature_bc_matrix/' 
counts = Read10X(data.dir = 'data/NSCLC_Yan/14728962/ST-visium/1A/outs/filtered_feature_bc_matrix/')  
data = CreateSeuratseuect(counts = counts, project = 'NSCLC_Yan', assay = 'Spatial')

seu_st1 <- Load10X_Spatial(data.dir = 'data/SKCM_Plozniak/Melanoma2_1/outs/')  
seu_st <- Load10X_Spatial(data.dir = 'data/NSCLC_Yan/14728962/ST-visium/1A/outs/')  
seu_st <- seu_st |> 
  SCTransform(assay = "Spatial", verbose = FALSE) |> 
  RunPCA(assay = "SCT", verbose = FALSE) |> 
  FindNeighbors(reduction = "pca", dims = 1:30) |> 
  FindClusters(verbose = FALSE) |> 
  RunUMAP(reduction = "pca", dims = 1:30)
  
seu <- qs_read("data/NSCLC_Yan/seu_final.qs2")
seu <- NormalizeData(seu) |> 
  FindVariableFeatures() |> 
  SketchData(
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch") 

seu <- seu |> 
  FindVariableFeatures() |> 
  ScaleData() |> 
  RunPCA() |> 
  FindNeighbors(reduction = "pca", dims = 1:30) |> 
  FindClusters() |> 
  RunUMAP(reduction = "pca", dims = 1:30)

## Rename the cells/spots with syntactically valid names
seu_st <- RenameCells(seu_st, new.names=make.names(Cells(seu_st)))
seu <- RenameCells(seu, new.names=make.names(Cells(seu)))

## Visualize the ST data
SpatialDimPlot(seu_st)
DimPlot(seu_st)

seu_traint <- CellTrek::traint(st_data=seu_st, sc_data=seu, sc_assay='sketch', st_assay = 'Spatial',cell_names='celltype_r2')


imgpath = "spatial"
img = Seurat::Read10X_Image(image.dir = imgpath)  

Seurat::DefaultAssay(seuect = img) <- 'Spatial'  

img = img[colnames(x = data)]  
data[['image']] = img  

SpatialFeaturePlot(data, features = "nCount_Spatial")