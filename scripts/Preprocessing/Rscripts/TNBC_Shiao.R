setwd("/home/zlin/workspace/PanCancer_ICI")
source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
# CD45-
seu_45neg <- schard::h5ad2seurat('./data/TNBC_Shiao/GSE246613_PembroRT_non_immune_cells.h5ad')
seu_45neg <- DietSeurat(seu_45neg, layers = 'counts', assays = 'RNA')
seu_45neg <- subset(seu_45neg, subset = treatment %in% c("Base", "PD1"))
seu_45neg$time_point <- ifelse(seu_45neg$treatment == 'Base', 'Pre', 'On')
seu_45neg$response <- ifelse(seu_45neg$response_group == 'NR', 'NR', 'RE')
seu_45neg@meta.data <- seu_45neg@meta.data |> select(time_point, response, cohort)

seu_45pos <- schard::h5ad2seurat('./data/TNBC_Shiao/GSE246613_PembroRT_immune_R100_final.h5ad')
seu_45pos <- DietSeurat(seu_45pos, layers = 'counts', assays = 'RNA')
seu_45pos <- subset(seu_45pos, subset = treatment %in% c("Base", "PD1"))
seu_45pos$time_point <- ifelse(seu_45pos$treatment == 'Base', 'Pre', 'On')
seu_45pos$response <- ifelse(seu_45pos$pCR == 'NR', 'NR', 'RE')
seu_45pos@meta.data <- seu_45pos@meta.data |> select(time_point, response, cohort)
seu <- merge(seu_45neg, y = seu_45pos, add.cell.ids = c("cd45n", "cd45p"))
rm(seu_45neg);rm(seu_45pos)

seu$patient <- seu$cohort
seu$cohort <- 'TNBC_Shiao'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$prior <- 'No'
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$treatment <- 'aPD1'
seu$modality <- 'Mono'
seu$res_metric <- 'Pathology'
seu$interval <- 20
seu$cancertype <- 'TNBC'
seu <- preprocessing(seu)
qs_save(seu, './data/TNBC_Shiao/processing.qs2')

seu <- qs_read('./data/TNBC_Shiao/processing.qs2')
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
marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
msig_df <- msigdbr(species = "Homo sapiens", category = 'C8')
msig_list <- split(x=geneset$gene_symbol, f=geneset$gs_name)
ora_res <- fora(pathways = msig_list, genes = head(marker_cosg$names$`21`, n=20), universe = rownames(seu), minSize = 1, maxSize = Inf)

seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '17'] <- 'Mast'
seu$celltype_major[seu$seurat_clusters %in% c('9','14')] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters %in% c('3','18')] <- 'Fibroblasts'
seu$celltype_major[seu$seurat_clusters == '15'] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters == '13'] <- 'Cycling T/NK'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('Epithelial','CD4T','CD8T','Bcell','PlasmaCell','panDC','Endothelial'), 
                                to = c('Epithelial cells','CD4+ T-cells','CD8+ T-cells','B-cells','Plasma cells','DC','Endothelial cells'))
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('Myocytes','Melanocytes'), 
                                to = c('Epithelial cells','Epithelial cells'))
marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('Macrophage','unknown','Monocyte'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()

seu <- qs_read('./data/TNBC_Shiao/seu_r1.qs2')

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

seu$percent.mito <- PercentageFeatureSet(seu, pattern = "^MT-")
seu$percent.ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

qs_save(seu, file = './data/TNBC_Shiao/seu_r1.qs2')

if (dir.exists('data/TNBC_Shiao/infercnv')==F){
  dir.create('data/TNBC_Shiao/infercnv')
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
  output_dir_full = paste0('data/TNBC_Shiao/infercnv/', pt)
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
folders <- list.files('data/TNBC_Shiao/infercnv')
infercnv_output <- lapply(folders, function(folder){
  print(folder)
  infercnv_obj <- readRDS(paste0('data/TNBC_Shiao/infercnv/',folder,'/run.final.infercnv_obj'))
  seu <- add_to_seurat(make_seurat_from_infercnv_obj(infercnv_obj), 
                       infercnv_output_path = paste0('data/TNBC_Shiao/infercnv/',folder), assay_name="RNA", top_n=10)
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
  ggsave(paste0('data/TNBC_Shiao/infercnv/',folder,'/boxplot.pdf'), height = 4, width = 5)
  return(data.frame('Malignant'=seu$malignant))
})
infercnv_output <- do.call(rbind, infercnv_output)
write.csv(infercnv_output, 'data/TNBC_Shiao/infercnv/infercnv_output.csv', row.names = T)




