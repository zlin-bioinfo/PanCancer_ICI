source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
input_dir <- "./data/TNBC_Zhang/"
count_matrix <- Read10X(input_dir,gene.column=1)
clinical_data <- readxl::read_xlsx("./data/TNBC_Zhang/1-s2.0-S1535610821004992-mmc2.xlsx", skip = 1) |>
  janitor::clean_names() |>
  slice(2:23) |>
  tidyr::fill(treatment) 
patient_id <- filter(clinical_data, tumor == 'Y', x11 == 'Y')$patient_id
count_matrix <- count_matrix[, str_split(colnames(count_matrix),'\\.', simplify = T)[,2] %in% paste0(rep(c('Pre_','Post_'),each=12),rep(paste0(patient_id, '_t'),2))]
seu <- CreateSeuratObject(counts = count_matrix, min.cells=5, min.features=400)
seu$sample <- paste0(str_split(str_split(rownames(seu@meta.data),'\\.', simplify = T)[,2],'_',simplify=T)[,1], '_',
                     str_split(str_split(rownames(seu@meta.data),'\\.', simplify = T)[,2],'_',simplify=T)[,2])
seu <- preprocessing(seu)
qs_save(seu, file = './data/TNBC_Zhang/processing.qs2')

seu <- qs_read('./data/TNBC_Zhang/processing.qs2')
seu <- seu |> FindClusters(resolution = 0.5)
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
seu$celltype_major[seu$seurat_clusters == '18'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters %in% c('6','15')] <- 'Plasma cells'
# seu$celltype_major[seu$seurat_clusters == '10' & seu$celltype_major == 'B-cells'] <- 'Cycling B'
# seu$celltype_major[seu$seurat_clusters == '10' & seu$celltype_major %in% c('CD4+ T-cells', 'CD8+ T-cells', 'NK cells')] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('CD4T','CD8T','NK'), 
                                to = c('CD4+ T-cells','CD8+ T-cells','NK cells'))
# marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('Macrophage','unknown','Monocyte','panDC','Bcell','Epithelial cells','Endothelial cells'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()

seu <- qs_read('./data/TNBC_Zhang/seu_r1.qs2')

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

seu$cohort <- 'TNBC_Zhang'
seu$time_point <- str_split(seu$sample, '_', simplify = T)[,1]
seu$patient <- paste0(seu$cohort, '_', str_split(seu$sample, '_', simplify = T)[,2])
seu$sample <- paste0(seu$patient, '_', seu$time_point)
clinical_data <- readxl::read_xlsx("./data/TNBC_Zhang/1-s2.0-S1535610821004992-mmc2.xlsx", skip = 1) |>
  janitor::clean_names() |>
  slice(2:23) |>
  tidyr::fill(treatment) 

seu$response <- clinical_data$clinical_efficacy_number[match(str_split(seu$patient, '_', simplify = T)[,3], clinical_data$patient_id)]
seu$treatment <- clinical_data$treatment[match(str_split(seu$patient, '_', simplify = T)[,3], clinical_data$patient_id)]
seu$treatment[seu$treatment == 'Anti-PD-L1+ Chemo'] <- 'aPDL1+Chemo'
seu$interval <- 28
seu$response <- ifelse(seu$response %in% c('PR', 'CR'), 'RE', 'NR')
seu$res_metric <- 'RECIST'
seu$modality <- 'Mono'
seu$prior <- 'No'

qs_save(seu, file = './data/TNBC_Zhang/seu_r1.qs2')


