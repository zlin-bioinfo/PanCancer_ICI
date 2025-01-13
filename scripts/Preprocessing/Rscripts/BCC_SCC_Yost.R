source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
matrix_count <- data.table::fread('./data/BCC_Yost/GSE123813_bcc_scRNA_counts.txt') |> 
  tibble::column_to_rownames(var = 'V1') |> as.sparse()
matrix_meta <- data.table::fread('./data/BCC_Yost/GSE123813_bcc_all_metadata.txt') |> 
  tibble::column_to_rownames(var = 'cell.id') 
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta, min.cells=5, min.features=400) |> 
  subset(subset = orig.ident %in% c('bcc.su010.pre.tcell','bcc.su010.post.tcell'), invert = T)
seu$cancertype <- 'BCC'
seu$res_metric <- 'RECIST'
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'treatment')] <- "time_point"
seu$time_point <- mapvalues(seu$time_point, from = c('pre','post'), to = c('Pre','On'))
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'cluster')] <- "celltype_orig"
seu$cohort <- 'BCC_Yost'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$prior <- 'Yes'
seu$prior[seu$patient == 'BCC_Yost_su004'] <- 'No'
seu$treatment <- 'aPD1'
seu$modality <- 'Mono'
rm(matrix_count);rm(matrix_meta);rm(matrix_tcr)
seu <- preprocessing(seu)
qs_save(seu, './data/BCC_Yost/processing.qs2')

seu <- qs_read('./data/BCC_Yost/processing.qs2')
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
                      c('KRT15','KRT17','EPCAM'), 
                      c('MKI67','TOP2A')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Neu','Fibro','PC','SMC','Endo','Mela','Epi','Proliferating')
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters))), label = T) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '11'] <- 'pDC'
seu$celltype_major[seu$seurat_clusters == '18'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters %in% c('7','17')] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '13'] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters == '12'] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('Epithelial','CD4T','CD8T','Bcell','PlasmaCell','Macrophage','Monocyte'), 
                                to = c('Epithelial cells','CD4+ T-cells','CD8+ T-cells','B-cells','Plasma cells','Macrophages','Monocytes'))
marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
seu <- subset(seu, subset = celltype_major %in% c('panDC','Endothelial','Fibroblast','unknown'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

matrix_tcr <- read.table('./data/BCC_Yost/GSE123813_bcc_tcr.txt')
seu$cdr3s_nt <- matrix_tcr[colnames(seu), 'cdr3s_nt']
seu$cdr3s_aa <- matrix_tcr[colnames(seu), 'cdr3s_aa']
seu$UMAP1 <- NULL
seu$UMAP2 <- NULL
# clinical
clinical <- readxl::read_xlsx('./data/BCC_Yost/41591_2019_522_MOESM2_ESM.xlsx', skip = 2, n_max = 16)
clinical <- clinical[clinical$Patient %in% unique(str_replace(seu$patient, 'BCC_Yost_','')),]
clinical$`scRNA days pre treatment`[clinical$Patient == 'su001'] <- -78
clinical$`scRNA days post treatment`[clinical$Patient == 'su003'] <- 121
clinical$interval <- as.numeric(clinical$`scRNA days post treatment`) 
seu$interval <- clinical$interval[match(str_replace(seu$patient, 'BCC_Yost_',''), clinical$Patient)]
seu$response <- clinical$Response[match(str_replace(seu$patient, 'BCC_Yost_',''), clinical$Patient)]
seu$response <- ifelse(seu$response == 'Yes', 'RE', 'NR')
qs_save(seu, file = './data/BCC_Yost/seu_r1.qs2')
seu <- qs_read('./data/BCC_Yost/seu_r1.qs2')

# SCC_Yost
matrix_count <- data.table::fread('./data/SCC_Yost/GSE123813_scc_scRNA_counts.txt') |>
  tibble::column_to_rownames(var = 'V1') |> as.sparse()
matrix_meta <- data.table::fread('./data/SCC_Yost/GSE123813_scc_metadata.txt') |>
  tibble::column_to_rownames(var = 'cell.id')
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta, min.cells=5, min.features=400)
seu$subtype <- 'SCC'
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'treatment')] <- "time_point"
seu$time_point <- mapvalues(seu$time_point, from = c('pre','post'), to = c('Pre','On'))
colnames(seu@meta.data)[which(colnames(seu@meta.data) == 'cluster')] <- "celltype_orig"
seu$cohort <- 'SCC_Yost'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$patient,'_', seu$time_point)
# su010 has 2 samples with distinct responses
# su013 doesn't have enough cell number in pre-Tx sample
seu <- subset(seu, subset = patient %in% c('SCC_Yost_su010', 'SCC_Yost_su013'), invert = T) 
seu$cancertype <- 'SCC'
seu$res_metric <- 'RECIST'
seu$prior <- 'Yes'
seu$treatment <- 'aPD1'
clinical <- readxl::read_xlsx('./data/BCC_Yost/41591_2019_522_MOESM2_ESM.xlsx', skip = 2, n_max = 16)
clinical <- clinical[clinical$Patient %in% unique(str_replace(seu$patient, 'SCC_Yost_','')),]
seu$interval <- clinical$interval[match(str_replace(seu$patient, 'SCC_Yost_',''), clinical$Patient)]
seu$response <- 'NR'
seu$response[seu$response == 'SCC_Yost_su011'] <- 'RE'
seu <- seu |> 
  StandardizeGeneSymbols(slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs) |> 
  NormalizeData() |> 
  CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
matrix_tcr <- read.table('./data/SCC_Yost/GSE123813_scc_tcr.txt')
seu$cdr3s_nt <- matrix_tcr[colnames(seu), 'cdr3s_nt']
seu$cdr3s_aa <- matrix_tcr[colnames(seu), 'cdr3s_aa']
clinical$interval <- as.numeric(clinical$`scRNA days post treatment`) 
seu$interval <- clinical$interval[match(str_replace(seu$patient, 'SCC_Yost_',''), clinical$Patient)]
qs_save(seu, file = './data/SCC_Yost/seu_r1.qs2')


