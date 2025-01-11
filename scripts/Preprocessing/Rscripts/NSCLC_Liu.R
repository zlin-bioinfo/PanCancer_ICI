source("./scripts/Preprocessing/Rscripts/Preprocessing.R")
matrix_count <- readRDS('./data/NSCLC_Liu/GSE179994_all.Tcell.rawCounts.rds')
matrix_meta <- data.table::fread('./data/NSCLC_Liu/GSE179994_Tcell.metadata.tsv') |> 
  tibble::column_to_rownames('cellid') 
seu <- CreateSeuratObject(matrix_count, meta.data = matrix_meta, min.cells = 5, min.features = 400)
matrix_tcr <- read.delim('./data/NSCLC_Liu/GSE179994_all.scTCR.tsv', row.names = 1)
seu$percent.mito <- PercentageFeatureSet(seu, pattern = "^MT-")
seu$percent.ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")

seu[['patient']] <- str_split(seu$sample, '\\.', simplify = T)[,1]
seu[['time_point']] <- str_split(seu$sample, '\\.', simplify = T)[,2]

seu$time_point <- ifelse(seu$time_point == 'pre', 'Pre', 'On')
seu$cohort <- 'NSCLC_Liu'

seu <- subset(seu, subset = percent.mito <20)
seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable = EnsemblGeneTable.Hs)
seu <- seu |> NormalizeData() 
# cell cycle
seu <- CellCycleScoring(seu, s.features = s.genes, g2m.features = g2m.genes)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
# automatic annotation
seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 10)
sce <- as.SingleCellExperiment(seu)
pred_bped_main <- SingleR(test = sce, ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(10))
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = sce, ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(10))
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
seu <- seu |> 
  FindVariableFeatures(nfeatures = 3000)  |>
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) |>
  RunPCA(verbose=FALSE) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony') |> 
  FindNeighbors(reduction = "harmony", dims = 1:20) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:20, reduction = 'harmony') |> 
  JoinLayers()

genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
                      c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('MKI67','TOP2A')
)
names(genes_to_check) <- c('T','NK','Proliferating')
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(seu, group.by = 'seurat_clusters', cols = getPalette(length(unique(seu$seurat_clusters))), label = T) /
  DotPlot(seu, group.by = 'seurat_clusters', features = genes_to_check) + RotatedAxis()

marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)

seu$celltype_major <- 'T cells'
seu$celltype_major[seu$celltype_bped_main == 'NK cells' & seu$scGate_multi ==  'NK'] <- 'NK cells'

seu <- subset(seu, subset = seurat_clusters %in% c(6,8), invert=T)
seu <- subset(seu, subset = scGate_multi == 'Multi' | 
                scGate_multi %in% c("panDC", "Fibroblast", "NK") | 
                celltype_bped_main == c('NK cells', "Epithelial cells", "B-cells", "Monocytes", "DC", "Fibroblasts",  "Macrophages", "Endothelial cells", "Melanocytes") |
                (celltype_bped_main == 'unknown' & scGate_multi == 'unknown'),
              invert=T)

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

qs_save(seu, '')
marker_genes <- c("CD4","CD8A","NKG7","CCR7", "RGCC", "MYADM", "LMNA", "ANXA1", "CCL5", "GZMK", "GZMA", "DUSP1","CD69", "FOS", "LY6E", "IFI6", "ISG15", "IFI27", "RPS27", "RPL41", "RPS29", "TCF7", "CXCL13", "IFNG", "PDCD1", "TOX", "CCR8", "LXN", "FOXP3", "TNFRSF9", "CRIP1", "UBA52", "TUBA1B", "TUBA1", "TYMS", "STMN1", "MKI67")
seu$cluster[is.na(seu$cluster)] <- 'NA'
DotPlot(seu, group.by = 'cluster', features = marker_genes) + RotatedAxis()

seu <- qs_read('./data/NSCLC_Liu/seu_r1.qs2')

pt_keep <- c('P1', 'P10', 'P13', 'P19', 'P29', 'P30', 'P33', 'P35')
seu <- subset(seu, subset = patient %in% pt_keep)
seu$interval <-  30
seu$interval[seu$patient == 'P1'] <- 213-14
seu$interval[seu$patient == 'P10'] <- 61-7
seu$interval[seu$patient == 'P13'] <- 94-17
seu$interval[seu$patient == 'P19'] <- 62-14
seu$interval[seu$patient == 'P29'] <- 55-8
seu$interval[seu$patient == 'P30'] <- 57-8
seu$interval[seu$patient == 'P33'] <- 81-14
seu$interval[seu$patient == 'P35'] <- 151-7
seu$cohort <- 'NSCLC_Liu'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$treatment <- 'aPD1'
seu$modality <- 'Mono'
seu$res_metric <- 'RECIST'
seu$response <- 'RE'
seu$prior <- 'No'

qs_save(seu, file = './data/NSCLC_Liu/seu_r1.qs2')








