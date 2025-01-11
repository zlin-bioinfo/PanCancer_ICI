source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
files <- list.files('./data/HNSC_Luoma/')
files_scrna_pre <- files[str_detect(files, 'pre-Tx_GEX_sc_tumor')]
pattern <- "P\\d+"
pt_rna <- str_extract(files_scrna_pre, pattern)
combined_pattern <- paste(paste0(pt_rna, '_post-Tx_GEX_sc_tumor'), collapse = "|")
files_scrna_post <- files[str_detect(files, combined_pattern)]
count_list <- purrr::map(paste0('./data/HNSC_Luoma/' ,c(files_scrna_pre, files_scrna_post)), Read10X_h5)
clin_df <- readxl::read_xlsx('./data/HNSC_Luoma/1-s2.0-S0092867422007231-mmc1.xlsx') |> 
  filter(`Pat. ID` %in% pt_rna) |> 
  tibble::column_to_rownames(var = 'Pat. ID')
clin_df_dup <- clin_df
rownames(clin_df_dup) <- paste0(rownames(clin_df), '_')
clin_df <- rbind(clin_df, clin_df_dup)
pt_rna_dup <- rep(pt_rna, times=2)
seu_list <- list()
timepoint <- rep(c('Pre','On'), each=6)
treatment <- rep(c('Nivo+Ipi', 'Nivo', 'Nivo+Ipi', 'Nivo', 'Nivo', 'Nivo+Ipi'), times = 2)
for (i in 1:length(count_list)){
  seu <- CreateSeuratObject(counts = count_list[[i]], min.cells=5, min.features=400)
  seu[['time_point']] <- rep(c('Pre','On'), each=6)[[i]]
  seu[['patient']] <- rep(pt_rna_dup, times=2)[[i]]
  seu[['sample']] <- paste0(seu$patient, '_', seu$time_point)
  seu[['treatment']] <- treatment[[i]]
  seu[['modality']] <- clin_df[i, 'Cohort']
  seu[['interval']] <- 28
  seu[['response']] <- clin_df[i, 'RECIST response excluding non measurable']
  seu_list[[i]] <- seu
}
seu <- merge(x=seu_list[[1]], y=seu_list[2:length(seu_list)])
seu <- JoinLayers(seu)

seu <- preprocessing(seu, sorted = T)
qs_save(seu, './data/HNSC_Luoma/processing.qs2')
seu <- qs_read('./data/HNSC_Luoma/processing.qs2')
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
seu$celltype_major[seu$seurat_clusters == '14'] <- 'pDC'
seu$celltype_major[seu$seurat_clusters == '15'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters == '7'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '13','18'] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major,
                                from = c('Epithelial','CD4T','CD8T','NK','Bcell','panDC','Fibroblast','Endothelial','Macrophage','Monocyte'),
                                to = c('Epithelial cells','CD4+ T-cells','CD8+ T-cells','NK cells','B-cells','DC','Fibroblasts','Endothelial cells','Macrophages','Monocytes'))
# marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('unknown','Epithelial cells','Endothelial cells','Fibroblasts','Myocytes','Melanocytes'), invert = T)

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

seu$treatment <- ifelse(seu$treatment == 'Nivo', 'aPD1','aPD1+CTLA4')
seu$response <- ifelse(seu$response == 'stable/progress', 'NR', NA)
seu$cohort <- 'HNSC_Luoma'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$prior <- 'No'
seu$res_metric <- 'RECIST'

qs_save(seu, file = './data/HNSC_Luoma/seu_r1.qs2')






