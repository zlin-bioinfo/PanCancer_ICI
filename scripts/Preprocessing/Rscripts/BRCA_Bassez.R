setwd("/home/zlin/workspace/PanCancer_ICI")
source("scripts/Preprocessing/Rscripts/Preprocessing.R")
# Bassez1(anti-PD1)
matrix_count <- readRDS('data/BRCA_Bassez1/1863-counts_cells_cohort1.rds')
matrix_meta <- read.csv('data/BRCA_Bassez1/1872-BIOKEY_metaData_cohort1_web.csv',row.names = 1)
matrix_tcr <- read.csv('data/BRCA_Bassez1/1879-BIOKEY_barcodes_vdj_combined_cohort1.csv')
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta, min.cells = 5, min.features = 400)
seu <- subset(seu, subset = expansion == "n/a", invert = TRUE)
seu[['cdr3s_nt']] <- matrix_tcr$cdr3_nt[match(str_split(colnames(seu),'-', simplify = T)[,1], matrix_tcr$barcode)]
extract_unique_tcr <- function(cdr3s_nt) {
  if (is.na(cdr3s_nt)) return(NA)
  tra <- str_extract_all(cdr3s_nt, "TRA:[^;]*")
  trb <- str_extract_all(cdr3s_nt, "TRB:[^;]*")
  if (length(tra[[1]]) == 1 & length(trb[[1]]) == 1) {
    return(cdr3s_nt)
  } else {
    return(NA)
  }
}
seu[['cdr3s_nt_unique']] <-  sapply(seu$cdr3s_nt, extract_unique_tcr)

colnames(seu@meta.data)[which(colnames(seu@meta.data)=='patient_id')] <- "patient"
seu$patient <- str_replace(seu$patient, 'BIOKEY', 'BRCA_Bassez1')
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='timepoint')] <- "time_point"
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='BC_type')] <- "subtype"
seu$response <- ifelse(seu$expansion == 'E', 'RE', 'NR')
seu$sample <- paste0(seu$patient,'_', seu$time_point)
seu$treatment <- 'aPD1'
seu$cohort <- 'BRCA_Bassez1'
seu$cancertype <- 'BRCA'
seu$res_metric <- 'T-cell expansion'
seu$site <- 'n/a'
seu$orig.ident <- NULL
seu$nCount_RNA <- NULL
seu$nFeature_RNA <- NULL
seu <- preprocessing(seu)
qs_save(seu, file = 'data/BRCA_Bassez1/processing.qs2')

seu <- qs_read('data/BRCA_Bassez1/processing.qs2')
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
                      c('KRT15','KRT17','EPCAM'),
                      c('MKI67','TOP2A')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Neu','Fibro','PC','SMC','Endo','Mela','Epi','Proliferating')
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters))), label = T) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '19'] <- 'pDC'
seu$celltype_major[seu$seurat_clusters == '18'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters == '21'] <- 'Doublets'
seu$celltype_major[seu$seurat_clusters == '11'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '10'] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters == '14'] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu@meta.data |> tabyl(celltype_major)
seu$celltype_major <- mapvalues(seu$celltype_major,
                                from = c('CD4T','CD8T','Fibroblast','Bcell'),
                                to = c('CD4+ T-cells','CD8+ T-cells','Fibroblasts','B-cells'))
# marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
seu <- subset(seu, subset = celltype_major %in% c('Macrophage','Epithelial','Endothelial','unknown','panDC','Melanocytes','Doublets','Myocytes'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()

seu <- qs_read('data/BRCA_Bassez1/seu_r1.qs2')
seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

seu[["percent.mito"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
seu$interval <- 9
seu$metric <- 'T cell Expansion'
seu$modality <- 'Mono'
seu$prior <- 'Yes'
qs_save(seu, file = 'data/BRCA_Bassez1/seu_r1.qs2')

#Bassez2(Chemo+anti-PD1)
matrix_count <- readRDS('data/BRCA_Bassez2/1867-counts_cells_cohort2.rds')
matrix_meta <- read.csv('data/BRCA_Bassez2/1871-BIOKEY_metaData_cohort2_web.csv',row.names = 1)
matrix_tcr <- read.csv('data/BRCA_Bassez2/1880-BIOKEY_barcodes_vdj_combined_cohort2.csv')
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta, min.cells = 5, min.features = 400)
seu[['cdr3s_nt']] <- matrix_tcr$cdr3_nt[match(str_split(colnames(seu),'-', simplify = T)[,1], matrix_tcr$barcode)]
extract_unique_tcr <- function(cdr3s_nt) {
  if (is.na(cdr3s_nt)) return(NA)
  tra <- str_extract_all(cdr3s_nt, "TRA:[^;]*")
  trb <- str_extract_all(cdr3s_nt, "TRB:[^;]*")
  if (length(tra[[1]]) == 1 & length(trb[[1]]) == 1) {
    return(cdr3s_nt)
  } else {
    return(NA)
  }
}
seu[['cdr3s_nt_unique']] <-  sapply(seu$cdr3s_nt, extract_unique_tcr)

colnames(seu@meta.data)[which(colnames(seu@meta.data)=='patient_id')] <- "patient"
seu$patient <- str_replace(seu$patient, 'BIOKEY', 'BRCA_Bassez2')
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='timepoint')] <- "time_point"
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='BC_type')] <- "subtype"
seu$response <- ifelse(seu$expansion == 'E', 'RE', 'NR')
seu$sample <- paste0(seu$patient,'_', seu$time_point)
seu$treatment <- 'aPD1 (NACT)'
seu$cohort <- 'BRCA_Bassez2'
seu$cancertype <- 'BRCA'
seu$res_metric <- 'T-cell expansion'

seu <- preprocessing(seu)
qs_save(seu, file = 'data/BRCA_Bassez2/processing.qs2')

seu <- qs_read('data/BRCA_Bassez2/processing.qs2')
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
                      c("KRT7","KRT19",'KRT15','KRT17','EPCAM'),
                      c('MKI67','TOP2A')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Neu','Fibro','PC','SMC','Endo','Mela','Epi','Proliferating')
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters))), label = T) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '18'] <- 'pDC'
seu$celltype_major[seu$scGate_multi == 'Mast'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters == '13'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '6'] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters == '12'] <- 'Cycling T/NK'
seu$celltype_major[(seu$seurat_clusters %in% c('3','11','14','15','16')) & (seu$celltype_bped_main == 'DC')] <- 'Epithelial cells'
# seu$celltype_major[(seu$seurat_clusters == '13') & (seu$celltype_bped_main == 'DC')] <- 'Plasma cells'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major,
                                from = c('CD4T','CD8T','Fibroblast','Bcell','Endothelial','Epithelial'),
                                to = c('CD4+ T-cells','CD8+ T-cells','Fibroblasts','B-cells','Endothelial cells','Epithelial cells'))


# marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
seu <- subset(seu, subset = celltype_major %in% c('unknown','Fibroblast','Macrophage') , invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()


seu <- qs_read(file = 'data/BRCA_Bassez2/seu_r1.qs2')
seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]
seu[["percent.mito"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
seu$treatment <- 'aPD1'
seu$interval <- 9
seu$res_metric <- 'T-cell expansion'
seu$prior <- 'Yes'
qs_save(seu, file = 'data/BRCA_Bassez2/seu_r1.qs2')


seu <- qs_read('data/BRCA_Bassez1/seu_r1.qs2')
celltype_ref <- c("Fibroblasts", "Monocytes", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils", "DC", "Mast","Mural cells")
gene_order <- read.table('data/hg38_gencode_v27.txt', header = F,row.names = 1)
lapply(unique(seu$patient), function(pt){
  seu_sub <- seu |>
    subset(subset = patient == pt) |>
    subset(subset = celltype_major %in% c("Fibroblasts", "Monocytes", "Epithelial cells", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils", "DC", "Mast","Mural cells"))
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=seu_sub@assays$RNA$counts,
                                      annotations_file=data.frame(row.names = colnames(seu_sub), 'Celltype' = seu_sub$celltype_major),
                                      delim="\t",
                                      gene_order_file=gene_order,
                                      ref_group_names=unique(seu_sub$celltype_major)[!unique(seu_sub$celltype_major) == "Epithelial cells"]
  )
  output_dir_full = paste0('data/BRCA_Bassez1/infercnv/', pt)
  infercnv_obj = suppressWarnings(infercnv::run(infercnv_obj,
                                                cutoff=0.1,
                                                out_dir=output_dir_full,
                                                cluster_by_groups=T,
                                                cluster_references = F,
                                                analysis_mode="subclusters",
                                                HMM=T,
                                                HMM_type='i3',
                                                denoise=T,
                                                plot_steps = F,
                                                num_threads = 20))
})
print('done')

make_seurat_from_infercnv_obj <- function(infercnv_obj) {
  return(CreateSeuratObject(counts = infercnv_obj@count.data, project="infercnv"))
}
folders <- list.files('data/BRCA_Bassez1/infercnv')
infercnv_output <- lapply(folders, function(folder){
  print(folder)
  infercnv_obj <- readRDS(paste0('data/BRCA_Bassez1/infercnv/',folder,'/run.final.infercnv_obj'))
  seu <- add_to_seurat(make_seurat_from_infercnv_obj(infercnv_obj), 
                       infercnv_output_path = paste0('data/BRCA_Bassez1/infercnv/',folder), assay_name="RNA", top_n=10)
  cnv_cols <- grep('proportion_cnv_chr', names(seu@meta.data), value = T)
  cnvs <- seu@meta.data[, cnv_cols]
  seu$proportion_cnv_avg <- rowMeans(cnvs)
  cnv_cols <- grep('has_cnv_chr', names(seu@meta.data), value = T)
  cnvs <- seu@meta.data[, cnv_cols]
  seu$has_cnv_avg <- rowMeans(cnvs)
  seu$celltype_major <- str_replace(seu$infercnv_subcluster, '_s\\d+','')
  seu$malignant <- 'no'
  seu$malignant[seu$celltype_major %in% c("Epithelial cells") & 
                  seu$has_cnv_avg > quantile(seu$has_cnv_avg[seu$celltype_major %in% celltype_ref], 0.9) & 
                  seu$proportion_cnv_avg > quantile(seu$proportion_cnv_avg[seu$celltype_major %in% celltype_ref], 0.9)] <- 'yes'
  # visualization
  seu@meta.data |>
    select(celltype_major, infercnv_subcluster, proportion_cnv_avg, has_cnv_avg) |>
    mutate(Celltype = case_when(celltype_major %in% celltype_ref ~ 'Ref',
                                celltype_major %in% c("Epithelial cells") ~ celltype_major)) |>
    tidyplot(x = Celltype, y = has_cnv_avg, color = Celltype) |>
    add_boxplot() |>
    add_test_pvalue(ref.group = 3) + RotatedAxis(45)
  ggsave(paste0('data/BRCA_Bassez1/infercnv/',folder,'/boxplot.pdf'), height = 4, width = 5)
  return(data.frame('Malignant'=seu$malignant))
})
infercnv_output <- do.call(rbind, infercnv_output)
write.csv(infercnv_output, 'data/BRCA_Bassez1/infercnv/infercnv_output.csv', row.names = T)


seu <- qs_read('data/BRCA_Bassez2/seu_r1.qs2')
if (dir.exists('data/BRCA_Bassez2/infercnv')==F){
  dir.create('data/BRCA_Bassez2/infercnv')
}
celltype_ref <- c("Fibroblasts", "Monocytes", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils", "DC", "Mast","Mural cells")
gene_order <- read.table('data/hg38_gencode_v27.txt', header = F,row.names = 1)
lapply(unique(seu$patient), function(pt){
  seu_sub <- seu |>
    subset(subset = patient == pt) |>
    subset(subset = celltype_major %in% c("Fibroblasts", "Monocytes", "Epithelial cells", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils", "DC", "Mast","Mural cells"))
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=seu_sub@assays$RNA$counts,
                                      annotations_file=data.frame(row.names = colnames(seu_sub), 'Celltype' = seu_sub$celltype_major),
                                      delim="\t",
                                      gene_order_file=gene_order,
                                      ref_group_names=unique(seu_sub$celltype_major)[!unique(seu_sub$celltype_major) == "Epithelial cells"]
  )
  output_dir_full = paste0('data/BRCA_Bassez2/infercnv/', pt)
  infercnv_obj = suppressWarnings(infercnv::run(infercnv_obj,
                                                cutoff=0.1,
                                                out_dir=output_dir_full,
                                                cluster_by_groups=T,
                                                cluster_references = F,
                                                analysis_mode="subclusters",
                                                HMM=T,
                                                HMM_type='i3',
                                                denoise=T,
                                                plot_steps = F,
                                                num_threads = 20))
})
print('done')

make_seurat_from_infercnv_obj <- function(infercnv_obj) {
  return(CreateSeuratObject(counts = infercnv_obj@count.data, project="infercnv"))
}
folders <- list.files('data/BRCA_Bassez2/infercnv')
infercnv_output <- lapply(folders, function(folder){
  print(folder)
  infercnv_obj <- readRDS(paste0('data/BRCA_Bassez2/infercnv/',folder,'/run.final.infercnv_obj'))
  seu <- add_to_seurat(make_seurat_from_infercnv_obj(infercnv_obj), 
                       infercnv_output_path = paste0('data/BRCA_Bassez2/infercnv/',folder), assay_name="RNA", top_n=10)
  cnv_cols <- grep('proportion_cnv_chr', names(seu@meta.data), value = T)
  cnvs <- seu@meta.data[, cnv_cols]
  seu$proportion_cnv_avg <- rowMeans(cnvs)
  cnv_cols <- grep('has_cnv_chr', names(seu@meta.data), value = T)
  cnvs <- seu@meta.data[, cnv_cols]
  seu$has_cnv_avg <- rowMeans(cnvs)
  seu$celltype_major <- str_replace(seu$infercnv_subcluster, '_s\\d+','')
  seu$malignant <- 'no'
  seu$malignant[seu$celltype_major %in% c("Epithelial cells") & 
                  seu$has_cnv_avg > quantile(seu$has_cnv_avg[seu$celltype_major %in% celltype_ref], 0.9) & 
                  seu$proportion_cnv_avg > quantile(seu$proportion_cnv_avg[seu$celltype_major %in% celltype_ref], 0.9)] <- 'yes'
  # visualization
  seu@meta.data |>
    select(celltype_major, infercnv_subcluster, proportion_cnv_avg, has_cnv_avg) |>
    mutate(Celltype = case_when(celltype_major %in% celltype_ref ~ 'Ref',
                                celltype_major %in% c("Epithelial cells") ~ celltype_major)) |>
    tidyplot(x = Celltype, y = has_cnv_avg, color = Celltype) |>
    add_boxplot() |>
    add_test_pvalue(ref.group = 3) + RotatedAxis(45)
  ggsave(paste0('data/BRCA_Bassez2/infercnv/',folder,'/boxplot.pdf'), height = 4, width = 5)
  return(data.frame('Malignant'=seu$malignant))
})
infercnv_output <- do.call(rbind, infercnv_output)
write.csv(infercnv_output, 'data/BRCA_Bassez2/infercnv/infercnv_output.csv', row.names = T)







