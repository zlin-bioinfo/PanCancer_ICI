#!/usr/bin/env Rscript
# setwd("/home/zlin/workspace/PanCancer_ICI")
source("./scripts/Preprocessing/Rscripts/Preprocessing.R")

samples <- list.files("./data/SKCM_Becker/cellranger")
samples <- samples[!samples %in% c('E19324','E23347','E23368')] # remove poor quality samples
seu_list <- lapply(samples, function(sample){
  print(sample)
  if (sample == 'E23346'){
    count_matrix <- Read10X(paste0("./data/SKCM_Becker/cellranger/", sample,"/outs/raw_feature_bc_matrix"))
    # count_matrix <- Matrix(as.matrix(count_matrix),sparse=TRUE)
  } else {
    count_matrix <- Read10X(paste0("./data/SKCM_Becker/cellranger/", sample,"/raw_feature_bc_matrix"))
  }
  seu <- CreateSeuratObject(counts = count_matrix, min.cells = 5, min.features = 400)
  if (ncol(seu) < 300) {
    print(paste("Sample", sample, "has fewer than", 300, "cells. Skipping."))
    return(NULL)  # Skip this samples
  }
  seu$sample <- sample
  seu$percent_mito <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu$percent_ribo <- PercentageFeatureSet(seu, pattern = "^RP[SL]")
  seu <- seu |> NormalizeData() |>
    FindVariableFeatures() |>
    ScaleData(vars.to.regress=c('nCount_RNA')) |>
    RunPCA() |>
    RunUMAP(dims = 1:20)
  # identify doublets
  # scDblFinder
  sce <- as.SingleCellExperiment(seu)
  sce <- scDblFinder(sce, dbr.sd=1)
  seu[['scDblFinder.class']] <- sce$scDblFinder.class
  # DoubletFinder
  sweep.list <- paramSweep_v5(seu, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn |>
    arrange(desc(BCmetric))
  pK <- pK[1, 'pK']
  pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
  nExp <- round(ncol(seu) * 0.08)
  seu <- doubletFinder_v5(seu, pN = 0.25, pK = pK, nExp = nExp, PCs = 1:10)
  seu[['DoubletFinder.class']] <- seu@meta.data[, colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]
  seu[[colnames(seu@meta.data)[grepl("DF.classification", colnames(seu@meta.data))]]] <- NULL
  seu[[colnames(seu@meta.data)[grepl("pANN", colnames(seu@meta.data))]]] <- NULL
  seu[['singlet']] <- ifelse(seu$scDblFinder.class=='doublet' & seu$DoubletFinder.class=='Doublet', 'no', 'yes')
  # filter low-quality cells and doublets
  seu <- subset(seu, subset = nFeature_RNA < 8000 & percent_mito < 20 & singlet =='yes')
  seu <- StandardizeGeneSymbols(seu, slot = 'counts', EnsemblGeneTable=EnsemblGeneTable.Hs)
  seu <- seu%>%
    NormalizeData() |>
    FindVariableFeatures() |>
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  # Preliminary annotation
  pred_bped_main <- SingleR(test = as.SingleCellExperiment(seu), ref = bped, labels = bped$label.main, BPPARAM=MulticoreParam(30))
  seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
  pred_bped_fine <- SingleR(test = as.SingleCellExperiment(seu), ref = bped, labels = bped$label.fine, BPPARAM=MulticoreParam(30))
  seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels
  # scGate
  seu <- scGate(seu, model = scGate_models_DB$human$TME_HiRes, ncores = 30)
  return(seu)
})
print('Preprocessing by sample done!')
names(seu_list) <- samples
qs_save(seu_list, file = './data/SKCM_Becker/list.qs2')

seu_list <- qs_read('./data/SKCM_Becker/list.qs2')
seu_list <- Filter(Negate(is.null), seu_list)
seu <- merge(x=seu_list[[1]], y=seu_list[2:length(seu_list)])
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
seu <- seu |> 
  NormalizeData() |>
  FindVariableFeatures(nfeatures = 3000)  |>
  ScaleData(vars.to.regress = c("S.Score", "G2M.Score")) |>
  RunPCA(verbose=FALSE) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony') |> 
  FindNeighbors(reduction = "harmony", dims = 1:20) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:20, reduction = 'harmony') |> 
  JoinLayers()
seu$celltype_bped_main[is.na(seu$celltype_bped_main)] <- 'unknown'
seu$celltype_bped_fine[is.na(seu$celltype_bped_fine)] <- 'unknown'
seu$scGate_multi[is.na(seu$scGate_multi)] <- 'unknown'
qs_save(seu, file = './data/SKCM_Becker/processing.qs2')

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
seu$celltype_major[seu$seurat_clusters == '19'] <- 'pDC'
seu$celltype_major[seu$seurat_clusters == '21'] <- 'Myocytes'
seu$celltype_major[seu$seurat_clusters == '17'] <- 'Plasma cells'
seu$celltype_major[seu$seurat_clusters == '20'] <- 'Pericytes'
seu$celltype_major[seu$seurat_clusters == '14'] <- 'Cycling T/NK'
seu$celltype_major[seu$seurat_clusters == '22'] <- 'Doublets'
seu$celltype_major[seu$seurat_clusters == '7'] <- 'Epithelial cells'
seu$celltype_major[seu$celltype_major == 'unknown'] <- seu$scGate_multi[seu$celltype_major == 'unknown']
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu$celltype_major <- mapvalues(seu$celltype_major, 
                                from = c('CD4T','CD8T','Endothelial','Macrophage','panDC','Bcell','NK'), 
                                to = c('CD4+ T-cells','CD8+ T-cells','Endothelial cells','Macrophages','DC','B-cells','NK cells'))
marker_cosg <- cosg(seu |> JoinLayers(), groups='all', assay='RNA', slot='data', mu=1, n_genes_user=100)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()
seu <- subset(seu, subset = celltype_major %in% c('Epithelial cells','Epithelial','panDC','unknown','Monocyte','Myocytes','Doublets','Fibroblast','Mast'), invert = T)
DimPlot(seu, group.by = 'celltype_major', cols = getPalette(length(unique(seu$celltype_major))), label = T) /
  DotPlot(seu, group.by = 'celltype_major', features = genes_to_check) + RotatedAxis()

seu@meta.data <- seu@meta.data[, !grepl("UCell", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("is.pure_", colnames(seu@meta.data))]
seu@meta.data <- seu@meta.data[, !grepl("CellOntology", colnames(seu@meta.data))]

clin_info <- readxl::read_xlsx('./data/SKCM_Becker/clin_info.xlsx')
seu$patient <- clin_info$patient[match(seu$sample, clin_info$sample)]
seu$time_point <- clin_info$time_point[match(seu$sample, clin_info$sample)]
seu$matched <- clin_info$matched[match(seu$sample, clin_info$sample)]
seu$interval <- clin_info$interval[match(seu$sample, clin_info$sample)]
seu$response <- clin_info$response[match(seu$sample, clin_info$sample)]
seu$modality <- clin_info$treatment[match(seu$sample, clin_info$sample)]
seu$treatment <- ifelse(seu$modality == 'Mono', 'aPD1', 'aPD1+CTLA4')
seu <- subset(seu, subset = sample == 'E23359', invert = T) # IL2 treatment
seu$sample_id <- seu$sample
seu$cohort <- 'SKCM_this study'
seu$patient <- paste0(seu$cohort, '_', seu$patient)
seu$sample <- paste0(seu$patient, '_', seu$time_point)
seu$res_metric <- 'RECIST'
seu$prior <- 'No'
qs_save(seu, file = './data/SKCM_Becker/seu_r1.qs2')

seu <- qs_read('./data/SKCM_Becker/seu_r1.qs2')

# DA(miloR)
# library(miloR)
# milo <- seu |> 
#   as.SingleCellExperiment() |> 
#   Milo() 
# milo
# milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "HARMONY")
# milo <- makeNhoods(milo, prop = 0.1, k = 30, d=30, refined = TRUE, reduced_dims = "HARMONY")
# plotNhoodSizeHist(milo)
# milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample="sample")
# 
# milo_design <- data.frame(colData(milo))[,c("sample", "patient", "time_point","response")]
# 
# ## Convert batch info from integer to factor
# milo_design$time_point <- factor(milo_design$time_point, levels = c('Pre','On'))
# milo_design$response <- as.factor(milo_design$response)
# milo_design$patient <- as.factor(milo_design$patient)
# milo_design$sample <- as.factor(milo_design$sample) 
# milo_design <- distinct(milo_design)
# rownames(milo_design) <- milo_design$sample
# 
# milo_design
# 
# milo <- calcNhoodDistance(milo, d=30, reduced.dim = "HARMONY")
# da_results <- testNhoods(milo, design = ~ response + time_point, design.df = milo_design, reduced.dim="HARMONY")
# head(da_results)
# 
# da_results <- annotateNhoods(milo, da_results, coldata_col = "celltype_major")
# plotDAbeeswarm(da_results, group.by = "celltype_major")
gene_order <- read.table('data/hg38_gencode_v27.txt', header = F,row.names = 1)
lapply(unique(seu$patient), function(pt){
  seu_sub <- seu |>
    subset(subset = patient == pt) |>
    subset(subset = celltype_major %in% c("Fibroblasts", "Monocytes", "Melanocytes", "Endothelial cells", "Macrophages", "NK cells", "pDC", "Pericytes", "Neutrophils"))
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=seu_sub@assays$RNA$counts,
                                      annotations_file=data.frame(row.names = colnames(seu_sub), 'Celltype' = seu_sub$celltype_major),
                                      delim="\t",
                                      gene_order_file=gene_order,
                                      ref_group_names=unique(seu_sub$celltype_major)[!unique(seu_sub$celltype_major) == "Melanocytes"]
  )
  output_dir_full = paste0('data/SKCM_Becker/infercnv/', pt)
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
folders <- list.files('data/SKCM_Becker/infercnv')
infercnv_output <- lapply(folders, function(folder){
  print(folder)
  infercnv_obj <- readRDS(paste0('data/SKCM_Becker/infercnv/',folder,'/run.final.infercnv_obj'))
  seu <- add_to_seurat(make_seurat_from_infercnv_obj(infercnv_obj), 
                       infercnv_output_path = paste0('data/SKCM_Becker/infercnv/',folder), assay_name="RNA", top_n=10)
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
  ggsave(paste0('data/SKCM_Becker/infercnv/',folder,'/boxplot.pdf'), height = 4, width = 5)
  return(data.frame('Malignant'=seu$malignant))
})
infercnv_output <- do.call(rbind, infercnv_output)
write.csv(infercnv_output, 'data/SKCM_Becker/infercnv/infercnv_output.csv', row.names = T)




