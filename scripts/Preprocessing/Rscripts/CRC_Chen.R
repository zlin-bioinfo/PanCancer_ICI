source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
count_matrix <- Read10X('./data/CRC_Chen')
metadata <- read.table('./data/CRC_Chen/CRC-ICB_metadata.txt.gz')
metadata$time_point <- str_split(metadata$Ident, '-', simplify = T)[,3]
metadata$interval[metadata$Patient == 'P01'] <- 26
metadata$interval[metadata$Patient == 'P02'] <- 20
metadata$interval[metadata$Patient == 'P03'] <- 106
metadata$interval[metadata$Patient == 'P04'] <- 21
metadata$interval[metadata$Patient == 'P05'] <- 110
metadata$interval[metadata$Patient == 'P08'] <- 21
metadata$interval[metadata$Patient == 'P09'] <- 110
metadata$interval[metadata$Patient == 'P11'] <- 72
metadata$interval[metadata$Patient == 'P12'] <- 29
metadata$interval[metadata$Patient == 'P14'] <- 62
metadata$interval[metadata$Patient == 'P15'] <- 65
metadata$interval[metadata$Patient == 'P16'] <- 87
metadata$interval[metadata$Patient == 'P17'] <- 28
metadata$interval[metadata$Patient == 'P18'] <- 96
metadata$interval[metadata$Patient == 'P19'] <- 71
metadata$interval[metadata$Patient == 'P20'] <- 41
metadata$interval[metadata$Patient == 'P21'] <- 100
metadata$interval[metadata$Patient == 'P22'] <- 40
metadata$interval[metadata$Patient == 'P23'] <- 51
metadata$interval[metadata$Patient == 'P24'] <- 68
metadata$interval[metadata$Patient == 'P25'] <- 57
metadata$interval[metadata$Patient == 'P26'] <- 60
metadata$time_point[metadata$time_point == 'I'] <- 'Pre'
metadata$time_point[metadata$time_point == 'II'] <- 'On'
metadata$time_point[metadata$Patient == 'P21' & metadata$time_point == 'III'] <- 'On'
clin_info <- readxl::read_xlsx('./data/CRC_Chen/1-s2.0-S1535610824002344-mmc2.xlsx', sheet = 1, skip = 1)
library(tibble)
metadata <- metadata |> 
  rownames_to_column(var = 'cell') |> 
  left_join(clin_info, by = join_by('Patient' == `Patient ID`)) |> 
  column_to_rownames(var = 'cell')
seu <- CreateSeuratObject(counts = count_matrix, meta.data = metadata, min.cells = 5, min.features = 400)
seu <- subset(seu, subset = Tissue == 'Tumor' & time_point %in% c('Pre', 'On'))
seu <- subset(seu, subset = Ident != 'CRC01-T2-II')
seu$cohort <- 'CRC_Chen'
seu$patient <- seu$Patient
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$treatment <- seu$`Treatment Regimen`
seu$treatment <- ifelse(str_detect(seu$treatment, 'CapeOx'), 'aPD1+Chemo', 'aPD1')
seu$modality <- 'Mono'
seu$res_metric <- 'iRECIST+pathology'
seu$response <- ifelse(seu$Response == 'CR', 'RE', 'NR')
seu$cancertype <- 'CRC'
seu$prior <- 'No'
seu$prior[seu$patient == 'P02'] <- 'Yes'

seu <- preprocessing(seu)
qs_save(seu, './data/CRC_Chen/processing.qs2')

seu <- qs_read('./data/CRC_Chen/processing.qs2')
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
seu$celltype_major[seu$seurat_clusters == '17'] <- 'pDC'
seu$celltype_major[seu$seurat_clusters == '14'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters == '3'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters %in% c('11','15')] <- 'Fibroblasts'
seu$celltype_major[seu$seurat_clusters == '13'] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
# seu@meta.data |> tabyl(celltype_major)
seu$celltype_major <- mapvalues(seu$celltype_major,
                                from = c('Epithelial','CD4T','CD8T','Bcell','Fibroblast','Endothelial','Macrophage','Monocyte'),
                                to = c('Epithelial cells','CD4+ T-cells','CD8+ T-cells','B-cells','Fibroblasts','Endothelial cells','Macrophages','Monocytes'))

# marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('unknown','Myocytes','Melanocytes'), invert = T)


seu <- qs_read('./data/CRC_Chen/seu_r1.qs2')
seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

seu <- subset(seu, subset = percent.mito < 20)
qs_save(seu, file = './data/CRC_Chen/seu_r1.qs2')




