source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
# IMCISION
matrix_count <- read.table('./data/HNSC_vanderLeun/GSM7324294_Count_data_IMCISION.txt') |> as.sparse()
matrix_meta <- read.table('./data/HNSC_vanderLeun/GSM7324295_Meta_data_IMCISION.txt', header = TRUE, sep = '\t', row.names = 1)
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta, min.cells=5, min.features=400)
seu$cohort <- 'HNSC_vanderLeun'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$time_point <- mapvalues(seu$timepoint, from = c('pre','post'), to = c('Pre','On'))
seu$sample <- paste0(seu$patient, '_', seu$time_point)

seu <- preprocessing(seu, sorted = T)
qs_save(seu, './data/HNSC_vanderLeun/processing.qs2')

seu <- qs_read('./data/HNSC_vanderLeun/processing.qs2')
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
seu@meta.data |> tabyl(celltype_bped_main, seurat_clusters)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '13'] <- 'pDC'
seu$celltype_major[seu$seurat_clusters == '17'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters %in% c('15','16')] <- 'DC'
seu$celltype_major[seu$seurat_clusters == '11'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '9'] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu@meta.data |> tabyl(celltype_major)
seu$celltype_major <- mapvalues(seu$celltype_major,
                                from = c('Epithelial','CD4T','CD8T','NK','Bcell','panDC','Fibroblast','Endothelial','Macrophage','Monocyte'),
                                to = c('Epithelial cells','CD4+ T-cells','CD8+ T-cells','NK cells','B-cells','DC','Fibroblasts','Endothelial cells','Macrophages','Monocytes'))

marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('unknown','Epithelial cells','Endothelial cells'), invert = T)



seu <- qs_read('./data/HNSC_vanderLeun/seu_r1.qs2')

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

seu$treatment <- 'aPD1+CTLA4'
seu$interval <- 28
seu$cohort <- 'HNSC_IMCISION'
seu$res_metric <- 'RECIST'
seu$prior <- 'No'
seu$modality <- 'Dual'

qs_save(seu, file = './data/HNSC_vanderLeun/seu_r1.qs2')










