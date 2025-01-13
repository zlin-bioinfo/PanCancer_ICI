setwd("/home/zlin/workspace/PanCancer_ICI")
source("scripts/Preprocessing/Rscripts/Preprocessing.R")
clin_info <- readxl::read_xlsx('data/SKCM_Plozniak/1-s2.0-S0092867423013223-mmc1.xlsx') |> 
  data.frame() |> 
  select(-Abbreviations)
pt_included <- unique(clin_info$Patient.ID) |> setdiff(c('16','25','29','39'))
seu <- readRDS('data/SKCM_Plozniak/Entire_TME.rds')
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
qs_save(seu, 'data/SKCM_Plozniak/processing.qs2')
gc()
seu <- qs_read('data/SKCM_Plozniak/processing.qs2')
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

qs_save(seu, file = 'data/SKCM_Plozniak/seu_r1.qs2')

seu <- qs_read('data/SKCM_Plozniak/seu_r1.qs2')

gene_order <- read.table('data/hg38_gencode_v27.txt', header = F,row.names = 1)
lapply(unique(seu$patient)[unique(seu$patient) != 'SKCM_Plozniak_41'], function(pt){
  tryCatch({
    seu_sub <- seu |>
      subset(subset = patient == pt) |>
      subset(subset = celltype_major %in% c("Fibroblasts", "Monocytes", "Melanocytes", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils"))
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=seu_sub@assays$RNA$counts,
                                        annotations_file=data.frame(row.names = colnames(seu_sub), 'Celltype' = seu_sub$celltype_major),
                                        delim="\t",
                                        gene_order_file=gene_order,
                                        ref_group_names=unique(seu_sub$celltype_major)[!unique(seu_sub$celltype_major) == "Melanocytes"]
    )
    output_dir_full = paste0('data/SKCM_Plozniak/infercnv/', pt)
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
  }, error = function(e) {
    message("Skipping iteration for x = ", x, ": ", e$message)
    return(NULL)  # Skip this iteration
  })
})
print('done')

make_seurat_from_infercnv_obj <- function(infercnv_obj) {
  return(CreateSeuratObject(counts = infercnv_obj@count.data, project="infercnv"))
}
folders <- list.files('data/SKCM_Plozniak/infercnv')
folders <- folders[folders != 'SKCM_Plozniak_41']
infercnv_output <- lapply(folders, function(folder){
  tryCatch({
    print(folder)
    infercnv_obj <- readRDS(paste0('data/SKCM_Plozniak/infercnv/',folder,'/run.final.infercnv_obj'))
    seu <- add_to_seurat(make_seurat_from_infercnv_obj(infercnv_obj), 
                         infercnv_output_path = paste0('data/SKCM_Plozniak/infercnv/',folder), assay_name="RNA", top_n=10)
    cnv_cols <- grep('proportion_cnv_chr', names(seu@meta.data), value = T)
    cnvs <- seu@meta.data[, cnv_cols]
    seu$proportion_cnv_avg <- rowMeans(cnvs)
    cnv_cols <- grep('has_cnv_chr', names(seu@meta.data), value = T)
    cnvs <- seu@meta.data[, cnv_cols]
    seu$has_cnv_avg <- rowMeans(cnvs)
    seu$celltype_major <- str_replace(seu$infercnv_subcluster, '_s\\d+','')
    seu$malignant <- 'no'
    seu$malignant[seu$celltype_major %in% c('Melanocytes') & 
                    seu$has_cnv_avg > quantile(seu$has_cnv_avg[seu$celltype_major %in% c("Fibroblasts", "Monocytes", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils")], 0.9) & 
                    seu$proportion_cnv_avg > quantile(seu$proportion_cnv_avg[seu$celltype_major %in% c("Fibroblasts", "Monocytes", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils")], 0.9)] <- 'yes'
    # visualization
    seu@meta.data |>
      select(celltype_major, infercnv_subcluster, proportion_cnv_avg, has_cnv_avg) |>
      mutate(Celltype = case_when(celltype_major %in% c("Fibroblasts", "Monocytes", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils") ~ 'Ref',
                                  celltype_major %in% c('Melanocytes') ~ celltype_major)) |>
      tidyplot(x = Celltype, y = has_cnv_avg, color = Celltype) |>
      add_boxplot() |>
      add_test_pvalue(ref.group = 3) + RotatedAxis(45)
    ggsave(paste0('data/SKCM_Plozniak/infercnv/',folder,'/boxplot.pdf'), height = 4, width = 5)
    return(data.frame('Malignant'=seu$malignant))
  }, error = function(e) {
    message("Skipping iteration for x = ", x, ": ", e$message)
    return(NULL)  # Skip this iteration
  })
})
infercnv_output <- do.call(rbind, infercnv_output)
write.csv(infercnv_output, 'data/SKCM_Plozniak/infercnv/infercnv_output.csv', row.names = T)




