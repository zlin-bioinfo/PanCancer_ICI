setwd("/home/zlin/workspace/PanCancer_ICI")
source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
matrix_count <- Read10X('./data/HNSC_Franken/HNSCC/')
matrix_meta <- read.csv('./data/HNSC_Franken/HNSCC_metadata.csv',row.names = 1)
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta, min.cells = 5, min.features = 400)
seu <- subset(seu, subset = SampleType %in% c('On-treatment', 'Pre-treatment'))
seu$Treatment <- mapvalues(seu$Treatment,from = c('Durvalumab','Durvalumab-Tremelimumab'), to = c('aPDL1','aPDL1+CTLA4'))
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='Treatment')] <- "treatment"
seu$SampleType <- mapvalues(seu$SampleType, from = c('Pre-treatment','On-treatment'), to = c('Pre','Post'))
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='SampleType')] <- "time_point"
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='Patient')] <- "patient"
seu$cohort <- 'HNSC_Franken'
seu$patient <- paste0(seu$cohort, '_', 'P', seu$patient)
seu$sample <- paste0(seu$patient,'_', seu$time_point)
seu$cancertype <- 'HNSC'
seu$res_metric <- 'T-cell expansion&MRI'
seu$subtype <- 'HNSC'
e <- paste0('HNSC_Franken_P', c(1,7,8,9,13,14,16,17,19,20))
ne <- paste0('HNSC_Franken_P', c(2,3,4,10,11,12,15,18))
seu$expansion <- NA
seu$expansion[seu$patient %in% e] <- 'Yes'
seu$expansion[seu$patient %in% ne] <- 'No'
df_mri <- data.frame(pt = paste0('HNSC_Franken_P', 1:20),
                     mri_pre = c(43,37,36,34,28,14,43,43,30,33,44,25,64,41,47,46,25,37,32,49),
                     mri_post = c(42,45,52,43,26,13,39,45,27,35,43,27,55,35,NA,27,25,49,34,47))
seu$mri_pre <- df_mri$mri_pre[match(seu$patient, df_mri$pt)]
seu$mri_post <- df_mri$mri_post[match(seu$patient, df_mri$pt)]
seu$mri_change <- seu$mri_post - seu$mri_pre
seu$mri_fc <- seu$mri_post/seu$mri_pre
seu$response <- 'NR'
seu$response[seu$expansion == 'Yes'] <- 'RE'
seu$res_metric <- 'T-cell expansion'

seu$interval <- 14
seu$modality <- ifelse(seu$treatment == 'aPDL1+CTLA4', 'Dual', 'Mono')
seu <- preprocessing(seu)
qs_save(seu, file = './data/HNSC_Franken/processing.qs2')

seu <- qs_read('./data/HNSC_Franken/processing.qs2')
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
seu$celltype_major[seu$seurat_clusters == '19'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters == '21'] <- 'Myocytes'
seu$celltype_major[seu$seurat_clusters == '12'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '13'] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters == '10'] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('Epithelial','CD4T','CD8T','Fibroblast'), 
                                to = c('Epithelial cells','CD4+ T-cells','CD8+ T-cells','Fibroblasts'))
marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('panDC','unknown','Melanocytes','Endothelial','Macrophage','Monocyte'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()

seu <- qs_read('./data/HNSC_Franken/seu_r1.qs2')

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

seu$response <- 'NR'
seu$response[seu$expansion == 'Yes'] <- 'RE'
seu$res_metric <- 'T-cell expansion'
seu$res_metric[seu$expansion == 'Yes' & seu$mri_change <= 0] <- 'T-cell expansion+MRI'
seu$response[seu$expansion == 'Yes' & seu$mri_change > 0] <- 'NR'
unique(seu$interval)

seu$prior <- 'Yes'
seu$treatment
qs_save(seu, file = './data/HNSC_Franken/seu_r1.qs2')

seu <- qs_read('./data/HNSC_Franken/seu_r1.qs2')
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
  output_dir_full = paste0('data/HNSC_Franken/infercnv/', pt)
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
folders <- list.files('data/HNSC_Franken/infercnv')
infercnv_output <- lapply(folders, function(folder){
  print(folder)
  infercnv_obj <- readRDS(paste0('data/HNSC_Franken/infercnv/',folder,'/run.final.infercnv_obj'))
  seu <- add_to_seurat(make_seurat_from_infercnv_obj(infercnv_obj), 
                       infercnv_output_path = paste0('data/HNSC_Franken/infercnv/',folder), assay_name="RNA", top_n=10)
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
  ggsave(paste0('data/HNSC_Franken/infercnv/',folder,'/boxplot.pdf'), height = 4, width = 5)
  return(data.frame('Malignant'=seu$malignant))
})
infercnv_output <- do.call(rbind, infercnv_output)
write.csv(infercnv_output, 'data/HNSC_Franken/infercnv/infercnv_output.csv', row.names = T)






