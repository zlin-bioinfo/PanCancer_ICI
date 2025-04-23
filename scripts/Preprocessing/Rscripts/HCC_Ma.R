setwd("/home/zlin/workspace/PanCancer_ICI")
source("scripts/Preprocessing/Rscripts/Preprocessing.R")
sample_info <- read.table('data/HCC_Ma/GSE151530_Info.txt', sep = '\t', header = T)
sample_info |> tabyl(Sample,S_ID)
sample_info <- tibble::column_to_rownames(sample_info, var = 'Cell')
count_mat <- Read10X('data/HCC_Ma/')
seu <- CreateSeuratObject(counts = count_mat, min.cells=5, min.features=400)
seu@meta.data <- cbind(seu@meta.data, sample_info[colnames(seu),])
seu@meta.data |> tabyl(Sample,S_ID)
seu <- subset(seu, subset = Sample %in% c('H73a','H73b','H68a','H68b','H58a','H58b','H49a','H49b','H34a','H34c'))
seu$sample <- seu$Sample
seu$patient <- gsub("[a-z]$", "", seu$sample)
seu$interval <- NA
seu$time_point <- ifelse(grepl('a',seu$sample), 'Pre', 'On')
seu <- preprocessing(seu)
qs_save(seu, 'data/HCC_Ma/processing.qs2')

seu <- qs_read('data/HCC_Ma/processing.qs2')
seu <- FindClusters(seu, resolution = 1)
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
                      c('ACTA2', "RGS5", "COX4I2","DCN",'SCN4B'),
                      c("DES", "TNNT3", "COX6A2", "ACTC1",  "MYL1"),
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('HP','ADH1B','KRT7','EPCAM','AFP'),
                      c('MKI67','TOP2A')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Neu','Fibro','PC','SMC','Endo','Mela','Epi','Proliferating')
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters))), label = T) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
marker_cosg <- COSG::cosg(seu, groups='all', assay='RNA', slot='data', mu=10, n_genes_user=50)
seu@meta.data |> tabyl(orig.ident, Tissues)
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters %in% c(20,21)] <- 'Doublets'
seu$celltype_major[seu$seurat_clusters == 23] <- 'low quality'
seu$celltype_major[seu$seurat_clusters == 22] <- 'pDC'
seu$celltype_major[seu$seurat_clusters == 13] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters %in% c(9,18)] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters %in% c(3,17,19)] <- 'Endothelial cells'
seu$celltype_major[seu$seurat_clusters == 15] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
seu$celltype_major[seu$seurat_clusters %in% c(6,8)] <- 'Epithelial cells'
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('unknown','low quality','Doublets'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]
seu$cohort <- 'HCC_Ma'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$treatment <- 'aPDL1/CTLA4'
seu$response <- 'NR'
seu$response[seu$patient %in% c('HCC_Ma_H73','HCC_Ma_H68','HCC_Ma_H34')] <- 'R'
seu$res_metric <- 'Unspecified'
seu$subset <- 'All TME'
seu$cancertype <- 'HCC'
seu <- qs_read('data/HCC_Ma/seu_r1.qs2')
qs_save(seu, file = 'data/HCC_Ma/seu_r1.qs2')

seu <- qs_read('data/HCC_Ma/seu_r1.qs2')
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
  output_dir_full = paste0('data/HCC_Ma/infercnv/', pt)
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
folders <- list.files('data/HCC_Ma/infercnv')
infercnv_output <- lapply(folders, function(folder){
  print(folder)
  infercnv_obj <- readRDS(paste0('data/HCC_Ma/infercnv/',folder,'/run.final.infercnv_obj'))
  seu <- add_to_seurat(make_seurat_from_infercnv_obj(infercnv_obj), 
                       infercnv_output_path = paste0('data/HCC_Ma/infercnv/',folder), assay_name="RNA", top_n=10)
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
  ggsave(paste0('data/HCC_Ma/infercnv/',folder,'/boxplot.pdf'), height = 4, width = 5)
  return(data.frame('Malignant'=seu$malignant))
})
infercnv_output <- do.call(rbind, infercnv_output)
write.csv(infercnv_output, 'data/HCC_Ma/infercnv/infercnv_output.csv', row.names = T)



