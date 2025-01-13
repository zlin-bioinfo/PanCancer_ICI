source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
seu <- readRDS('./data/PCa_Hawley/primecut.integrated.geneexp.seurat.rds')
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data)
seu <- subset(seu, subset = patient %in% c('Patient1', 'Patient7'))
seu$cohort <- 'PCa_Hawley'
seu$cancertype <- 'PC'
seu$res_metric <- 'PSA test'
seu$response <- 'NR'
names(seu@meta.data)[names(seu@meta.data) == 'treatment'] <- 'time_point'
seu$time_point <- as.character(seu$time_point)
seu$time_point <- ifelse(seu$time_point == 'Pre-Treatment', 'Pre', 'On')
seu$interval <- 70
seu$prior <- 'Yes'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$cohort, '_', seu$time_point)
seu$treatment <- 'aPD1'
seu$modality <- 'Mono'

seu <- preprocessing(seu)
qs_save(seu, './data/PCa_Hawley/processing.qs2')

seu <- qs_read('./data/PCa_Hawley/processing.qs2')
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
seu$celltype_major[seu$seurat_clusters == '14'] <- 'pDC'
seu$celltype_major[seu$seurat_clusters == '13'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '11'] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major,
                                from = c('CD4T','CD8T','panDC','Bcell','Fibroblast','Endothelial','Macrophage'),
                                to = c('CD4+ T-cells','CD8+ T-cells','DC','B-cells','Fibroblasts','Endothelial cells','Macrophages'))
# marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('unknown','Monocyte','Epithelial','NK','DC','Neutrophils'), invert = T)

seu <- qs_read('./data/PCa_Hawley/seu_r1.qs2')

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]
qs_save(seu, file = './data/PCa_Hawley/seu_r1.qs2')


