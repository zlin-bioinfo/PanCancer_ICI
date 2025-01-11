# in terminal
# ls GSM* | awk -F '_' '{print $1}' | uniq| while read i; do mkdir $i; mv *$i*gz $i; done
# find -name "*matrix.mtx.gz" | while read i; do mv $i $(dirname $i)/matrix.mtx.gz; done
# find -name "*features.tsv.gz" | while read i; do mv $i $(dirname $i)/features.tsv.gz; done
# find -name "*barcodes.tsv.gz" | while read i; do mv $i $(dirname $i)/barcodes.tsv.gz; done
source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
sample_list = list.files("./data/CRC_Li/", pattern = 'GSM')
input_dir_list = paste0(rep("./data/CRC_Li/", 14), sample_list)
seu_list <- lapply(sample_list, function(sample) {
  input_dir <- input_dir_list[str_detect(input_dir_list, sample)]
  count_matrix <- Read10X(input_dir)
  seu <- CreateSeuratObject(counts = count_matrix, min.cells=5, min.features=400)
  seu$SampleID <- sample
  return(seu)
})
names(seu_list) <- sample_list
qs_save(seu_list,'./data/CRC_Li/raw.qs2')

seu_list <- qs_read('./data/CRC_Li/raw.qs2')
sample_list = list.files("./data/CRC_Li/", pattern = 'GSM')
meta_df <- data.frame(row.names = sample_list, 
                      patient = rep(c('P21','P24','P25','P27','P28','P30','P31'), each=2), 
                      time_point = rep(c('On','Pre'),7),
                      treatment = rep('aPD1',14),
                      response = rep('pCR', 14))
meta_df$response[meta_df$patient == 'P31'] <- 'non-pCR'
meta_df$sample <- paste0(meta_df$patient, '_', meta_df$time_point)
# merge samples
for (i in 1:14){
  seu_list[[i]]$sample <- meta_df[i,'sample']
  seu_list[[i]]$response <- meta_df[i,'response']
  seu_list[[i]]$treatment <- meta_df[i,'treatment']
  seu_list[[i]]$patient <- meta_df[i,'patient']
  seu_list[[i]]$time_point <- meta_df[i,'time_point']
}
seu <- merge(x=seu_list[[1]], y=seu_list[2:length(seu_list)]) |> JoinLayers()
seu <- preprocessing(seu)
qs_save(seu, './data/CRC_Li/processing.qs2')
seu <- qs_read('./data/CRC_Li/processing.qs2')
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

seu$celltype_major <- seu$celltype_bped_main
seu$celltype_major[seu$seurat_clusters == '13'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '14'] <- 'Fibroblasts'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('Epithelial'), 
                                to = c('Epithelial cells'))
# marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('unknown','Neutrophils','Myocytes','Melanocytes'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()

qs_save(seu, file = './data/CRC_Li/seu_r1.qs2')

seu <- qs_read('./data/CRC_Li/seu_r1.qs2')
seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

seu$cohort <- 'CRC_Li'
seu$patient <- paste0(seu$cohort, '_',seu$patient)
seu$sample <- paste0(seu$cohort, '_',seu$sample)
seu$response <- ifelse(seu$response == 'pCR', 'RE','NR')
seu$interval <- 84
seu$res_metric <- 'Pathology'
seu$modality <- 'Mono'
seu$prior <- 'Unknown'
seu <- seu |> subset(subset = percent.mito < 20)

seu <- qs_save(seu, './data/CRC_Li/seu_r1.qs2')



