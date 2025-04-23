source("scripts/Preprocessing/Rscripts/Preprocessing.R")
count_mat <- Read10X('data/RCC_Bi/')
meta_data <- read.table('data/RCC_Bi/metadata.txt', header = T, sep = '\t', row.names = 'NAME')
meta_data <- meta_data[-1,-which(colnames(meta_data) %in% c('organ','organ__ontology_label','library_preparation_protocol__ontology_label','library_preparation_protocol','Initial_Louvain_Cluster'))]
meta_data$cohort <- 'RCC_Bi'
meta_data$patient <- paste0(meta_data$cohort, '_', meta_data$donor_id)
meta_data$sample <- meta_data$patient
meta_data <- meta_data |> distinct(patient, .keep_all = T)
seu <- CreateSeuratObject(counts = count_mat, min.cells=5, min.features=400)
seu$interval <- NA
seu$res_metric <- 'RECIST'
seu$patient <- paste0('RCC_Bi_', str_split(colnames(seu),'\\.', simplify=T)[,2] |> str_replace('p','P'))
seu$sample <- seu$patient
seu$response <- meta_data$ICB_Response[match(seu$patient, meta_data$patient)]
seu$sex <- meta_data$sex[match(seu$patient, meta_data$patient)]
seu$time_point <- ifelse(seu$response == 'NoICB', 'NoICB', 'ICB')
seu <- preprocessing(seu)
seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]
qs_save(seu, 'data/RCC_Bi/processing.qs2')

seu <- qs_read('data/RCC_Bi/processing.qs2')
seu <- seu |> FindClusters(resolution = 1)
genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('CD27','CD38','JCHAIN'), # Plasma cells 
                      c('LILRA4','IL3RA','PLD4'),
                      c('KIT','TPSAB1','CPA3'),
                      c('CLEC9A','FCER1A','CD1C','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('CXCR1', 'CXCR2', 'PTGS2','OLR1', 'VEGFA'),
                      c('COL3A1', 'FAP', 'COL1A1'), 
                      c('ACTA2', "RGS5", "COX4I2","DCN"),
                      c("DES", "TNNT3", "COX6A2", "ACTC1",  "MYL1"),
                      c('PECAM1','VWF', 'ENG'), 
                      c('MLANA','MITF', 'TYR'), 
                      c('KRT15','KRT17','KRT19','EPCAM', "ACSM2A",  "SLC22A5", "ACSM2B"),
                      c('MKI67','TOP2A')
)
names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Neu','Fibro','PC','SMC','Endo','Mela','Epi','Proliferating')
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters))), label = T) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters %in% c(4,7,9,10,12,17,19)] <- 'Malignant'
seu$celltype_major[seu$seurat_clusters == '20'] <- 'DC'
seu$celltype_major[seu$seurat_clusters == '22'] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters == '16'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '13'] <- 'Cycling T/NK'
seu$celltype_major[seu$seurat_clusters == '21'] <- 'Cycling myeloids'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('CD4T'), 
                                to = c('CD4+ T-cells'))
seu <- subset(seu, subset = celltype_major %in% c('CD8T','unknown','Melanocytes','Endothelial','Macrophage','Monocyte','Fibroblasts','Epithelial cells','Myocytes'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$cohort <- 'RCC_Bi'
seu$treatment <- meta_data$ICB_Exposed[match(seu$patient, meta_data$patient)]
seu$treatment <- ifelse(seu$treatment == 'ICB', 'aPD1', 'NoICI')
seu$time_point <- ifelse(seu$time_point == 'ICB', 'ICI_exposed', 'NoICI')
qs_save(seu, 'data/RCC_Bi/seu_r1.qs2')

seu <- qs_read('./data/RCC_Bi/seu_r1.qs2')
seu$treatment[seu$treatment == 'aPD1/L1'] <- 'aPD1'
seu$treatment[seu$patient == 'RCC_Bi_P915'] <- 'aPD1+CTLA4'
celltype_ref <- c("Fibroblasts", "Monocytes", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils", "DC", "Mast","Mural cells")
gene_order <- read.table('data/hg38_gencode_v27.txt', header = F,row.names = 1)
lapply(unique(seu$patient), function(pt){
  seu_sub <- seu |>
    subset(subset = patient == pt) |>
    subset(subset = celltype_major %in% c("Fibroblasts", "Monocytes", "Malignant", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils", "DC", "Mast","Mural cells"))
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=seu_sub@assays$RNA$counts,
                                      annotations_file=data.frame(row.names = colnames(seu_sub), 'Celltype' = seu_sub$celltype_major),
                                      delim="\t",
                                      gene_order_file=gene_order,
                                      ref_group_names=unique(seu_sub$celltype_major)[!unique(seu_sub$celltype_major) == "Malignant"]
  )
  output_dir_full = paste0('data/RCC_Bi/infercnv/', pt)
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
folders <- list.files('data/RCC_Bi/infercnv')
infercnv_output <- lapply(folders, function(folder){
  print(folder)
  infercnv_obj <- readRDS(paste0('data/RCC_Bi/infercnv/',folder,'/run.final.infercnv_obj'))
  seu <- add_to_seurat(make_seurat_from_infercnv_obj(infercnv_obj), 
                       infercnv_output_path = paste0('data/RCC_Bi/infercnv/',folder), assay_name="RNA", top_n=10)
  cnv_cols <- grep('proportion_cnv_chr', names(seu@meta.data), value = T)
  cnvs <- seu@meta.data[, cnv_cols]
  seu$proportion_cnv_avg <- rowMeans(cnvs)
  cnv_cols <- grep('has_cnv_chr', names(seu@meta.data), value = T)
  cnvs <- seu@meta.data[, cnv_cols]
  seu$has_cnv_avg <- rowMeans(cnvs)
  seu$celltype_major <- str_replace(seu$infercnv_subcluster, '_s\\d+','')
  seu$malignant <- 'no'
  seu$malignant[seu$celltype_major %in% c("Malignant") & 
                  seu$has_cnv_avg > quantile(seu$has_cnv_avg[seu$celltype_major %in% celltype_ref], 0.9) & 
                  seu$proportion_cnv_avg > quantile(seu$proportion_cnv_avg[seu$celltype_major %in% celltype_ref], 0.9)] <- 'yes'
  # visualization
  seu@meta.data |>
    select(celltype_major, infercnv_subcluster, proportion_cnv_avg, has_cnv_avg) |>
    mutate(Celltype = case_when(celltype_major %in% celltype_ref ~ 'Ref',
                                celltype_major %in% c("Malignant") ~ celltype_major)) |>
    tidyplot(x = Celltype, y = has_cnv_avg, color = Celltype) |>
    add_boxplot() |>
    add_test_pvalue(ref.group = 3) + RotatedAxis(45)
  ggsave(paste0('data/RCC_Bi/infercnv/',folder,'/boxplot.pdf'), height = 4, width = 5)
  return(data.frame('Malignant'=seu$malignant))
})
infercnv_output <- do.call(rbind, infercnv_output)
write.csv(infercnv_output, 'data/RCC_Bi/infercnv/infercnv_output.csv', row.names = T)




