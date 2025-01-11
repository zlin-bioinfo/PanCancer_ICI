source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
clin_info <- readxl::read_xlsx('./data/SKCM_Plozniak/1-s2.0-S0092867423013223-mmc1.xlsx') |> 
  data.frame() |> 
  select(-Abbreviations)
pt_included <- unique(clin_info$Patient.ID) |> setdiff(c('16','25','29','39'))
seu <- readRDS('./data/SKCM_Plozniak/Entire_TME.rds')
seu <- CreateSeuratObject(counts = seu@assays$RNA$counts, meta.data = seu@meta.data, min.cells = 5, min.features = 400)
seu$percent_mito <- PercentageFeatureSet(seu, pattern = "^MT-")
seu$percent_ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
colnames(seu@meta.data)[5] <- 'patient'
seu <- subset(seu, subset = patient %in% pt_included)
seu$cohort <- 'SKCM_Plozniak'
colnames(seu@meta.data)[4] <- 'time_point'
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

seu <- preprocessing(seu)
qs_save(seu, './data/SKCM_Plozniak/processing.qs2')
gc()
seu <- qs_read('./data/SKCM_Plozniak/processing.qs2')
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
                      c('ACTA2', "PDGFRB","RGS5", "COX4I2","DCN"),
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
seu$celltype_major[seu$seurat_clusters == '16'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters == '14'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '7'] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters == '8'] <- 'Cycling T/NK'
seu$celltype_major[seu$seurat_clusters %in% c(0,3,9,11,13)] <- 'Melanocytes'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('Epithelial','CD4T','CD8T','Endothelial','Macrophage'), 
                                to = c('Epithelial cells','CD4+ T-cells','CD8+ T-cells','Endothelial cells','Macrophages'))
marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('panDC','unknown','Monocyte','Fibroblast','Myocytes'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$res_metric <- ifelse(seu$res_metric == 'RECIST', 'RECISTv1.1',seu$res_metric)
seu$interval <- round(2.5*7)
seu <- subset(seu, subset = percent.mito < 20)

qs_save(seu, file = './data/SKCM_Plozniak/seu_r1.qs2')

seu <- qs_read('./data/SKCM_Plozniak/seu_r1.qs2')





