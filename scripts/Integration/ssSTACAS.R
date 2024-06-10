#!/usr/bin/env Rscript
rm(list=ls())
pkgs <- c('Seurat','dplyr','stringr','ggsci','qs','RColorBrewer','STACAS','SignatuR','SeuratDisk','SCP','scRNAtoolVis')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

# Loading dataset
SKCM_Becker <- qread('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/seu_r2.qs')
BRCA_Bassez1 <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r2.qs')
BRCA_Bassez2 <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/seu_r2.qs')
TNBC_Zhang <- qread('/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/seu_r2.qs')
BCC_Yost <- qread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/seu_r2.qs')
SCC_Yost <- qread('/bigdata/zlin/Melanoma_meta/data/SCC_Yost/seu_r2.qs')
BCC_SCC_Yost <- merge(x=BCC_Yost, y=SCC_Yost)
BCC_SCC_Yost[['RNA']] <- JoinLayers(BCC_SCC_Yost[['RNA']])
HNSC_IMCISION <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_IMCISION/seu_r2.qs')
HNSC_Luoma <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/seu_r2.qs')
NSCLC_Liu <- qread('/bigdata/zlin/Melanoma_meta/data/NSCLC_Liu/seu_r2.qs')
CRC_Li <- qread('/bigdata/zlin/Melanoma_meta/data/CRC_Li/seu_r2.qs')
PCa_Hawley <- qread('/bigdata/zlin/Melanoma_meta/data/PCa_Hawley/seu_r2.qs')
datasets <- list(SKCM_Becker, BRCA_Bassez1, BRCA_Bassez2, TNBC_Zhang, BCC_SCC_Yost, HNSC_IMCISION, HNSC_Luoma, NSCLC_Liu, CRC_Li, PCa_Hawley)

# seu <- merge(x=datasets[[1]], y=datasets[2:length(datasets)], add.cell.ids = names(datasets)) 
# nfeatures <- 2000
# ndim <- 20
# genes.blocklist <- c(GetSignature(SignatuR$Hs$Blocklists),
#                      GetSignature(SignatuR$Hs$Compartments))
# print('start integration')
# seu <- seu %>% 
#   FindVariableFeatures(nfeatures = nfeatures) %>%
#   NormalizeData() %>% ScaleData() %>%
#   # RunPCA(npcs=ndim) %>% RunUMAP(dims=1:ndim) %>% 
#   SplitObject(split.by = "dataset") %>% 
#   Run.STACAS(genesBlockList = genes.blocklist, cell.labels = "celltype_main",
#              dims = 1:ndim, anchor.features = 2000) %>% RunUMAP(dims=1:ndim)
# 
# qsave(seu, '/bigdata/zlin/Melanoma_meta/data/ssSTACAS_integrated.qs')
# print('done')
# 
# seu <- qread('/bigdata/zlin/Melanoma_meta/data/ssSTACAS_integrated.qs')
# DefaultAssay(seu) <- 'RNA'
# seu <- JoinLayers(seu)
# seu <- seu %>% 
#   FindVariableFeatures(nfeatures = 2000) %>%
#   NormalizeData() %>% ScaleData()
# genes_to_check = list(c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B'), # T cells 'CD8B'
#                       c('KLRD1','KLRB1', 'KLRC1', 'NCAM1'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
#                       c('CD79A','CD19', 'MS4A1'),  # B cells 
#                       c('CD27','CD38'), # Plasma cells 
#                       c('LILRA4','IL3RA','PLD4'),
#                       c('KIT','TPSAB1','CPA3'),
#                       c('CLEC9A','FCER1A','LAMP3'), 
#                       c('CD68', 'LYZ', 'CD14'),  
#                       c('COL3A1','FAP', 'COL1A1','ATCA2'), 
#                       c('PECAM1','VWF', 'ENG')
# )
# names(genes_to_check) <- c('T','NK','B','Plasma','pDC','Mast','cDC','Mo/Mac','Fibro','Endo')
# colors <- c('CD4+T' = '#E31A1C',
#             'CD8+T' = '#1F78B4',
#             'NK' = '#A6CEE3',
#             'B' = '#CAB2D6',
#             'Plasma' = '#6A3D9A',
#             'pDC' = '#FF7F00',
#             'Mast' = '#B15928',
#             'cDC' = '#FFFF99',
#             'Mono' = '#B2DF8A',
#             'Macro' = '#33A02C',
#             'Endo' = '#FB9A99',
#             'CAF' = '#FDBF6F')
# Idents(seu) <- factor(seu$celltype_main, levels = c('CD4+T', 'CD8+T', 'NK', 'B', 'Plasma', 'pDC', 'Mast', 'cDC', 'Mono', 'Macro', 'Endo', 'CAF'))
# CellDimPlot(seu, group.by = "celltype_main", reduction = "umap", theme_use = "theme_blank", palcolor = colors)
# DimPlot(seu, group.by = "celltype_main", reduction = "umap", theme_use = "theme_blank", palcolor = unlist(colors))
# jjDotPlot(object = seu,
#           gene = c('CD3D','CD8A','GNLY','MS4A1','JCHAIN','LILRA4','TPSAB1','FCER1A','LYZ','C1QA','VWF','COL1A1'),
#           id = 'celltype_main', ytree = F,
#           cluster.order = rev(c('CD4+T', 'CD8+T', 'NK', 'B', 'Plasma', 'pDC', 'Mast', 'cDC', 'Mono', 'Macro', 'Endo', 'CAF')))
# 
# seu <- JoinLayers(seu)
# seu <- seu %>% 
#   FindVariableFeatures(nfeatures = 2000) %>%
#   NormalizeData() %>% ScaleData()
# DotPlot(seu, group.by = "celltype_main",  
#         features = c('CD3D','CD8A','GNLY','MS4A1','JCHAIN','LILRA4','TPSAB1','FCER1A','LYZ','C1QA','VWF','COL1A1'))
# data <- GetAssayData(seu_int[["integrated"]], slot = "scale.data")
# pcs <- prcomp(x = data)
# umap.dr <- CreateDimReducObject(
#   embeddings = pcs$rotation,
#   loadings = pcs$x,
#   stdev = pcs$sdev,
#   key = "UMAP",
#   assay = "RNA"
# )
# seu[['umap']] <- umap.dr
# 

nfeatures <- 1000
ndim <- 20

seu <- merge(x=datasets[[1]], y=datasets[2:length(datasets)], add.cell.ids = names(datasets)) %>% subset(subset = celltype_main == 'CD8+T')
genes.blocklist <- c(GetSignature(SignatuR$Hs$Blocklists),
                     GetSignature(SignatuR$Hs$Compartments))
print('start integration')
seu <- seu %>% 
  FindVariableFeatures(nfeatures = nfeatures) %>%
  NormalizeData() %>% ScaleData() %>%
  # RunPCA(npcs=ndim) %>% RunUMAP(dims=1:ndim) %>% 
  SplitObject(split.by = "dataset") %>% 
  Run.STACAS(genesBlockList = genes.blocklist, cell.labels = "celltype_r2",
             dims = 1:ndim, anchor.features = nfeatures) %>% RunUMAP(dims=1:ndim)
DefaultAssay(seu) <- 'RNA'
seu <- JoinLayers(seu)
seu <- seu %>% 
  FindVariableFeatures(nfeatures = nfeatures) %>%
  NormalizeData() %>% ScaleData()
qsave(seu, '/bigdata/zlin/Melanoma_meta/data/ssSTACAS_cd8t.qs')
print('done')

seu <- merge(x=datasets[[1]], y=datasets[2:length(datasets)], add.cell.ids = names(datasets)) %>% subset(subset = celltype_main == 'CD4+T')
genes.blocklist <- c(GetSignature(SignatuR$Hs$Blocklists),
                     GetSignature(SignatuR$Hs$Compartments))
print('start integration')
seu <- seu %>% 
  FindVariableFeatures(nfeatures = nfeatures) %>%
  NormalizeData() %>% ScaleData() %>%
  # RunPCA(npcs=ndim) %>% RunUMAP(dims=1:ndim) %>% 
  SplitObject(split.by = "dataset") %>% 
  Run.STACAS(genesBlockList = genes.blocklist, cell.labels = "celltype_r2",
             dims = 1:ndim, anchor.features = nfeatures) %>% RunUMAP(dims=1:ndim)
DefaultAssay(seu) <- 'RNA'
seu <- JoinLayers(seu)
seu <- seu %>% 
  FindVariableFeatures(nfeatures = nfeatures) %>%
  NormalizeData() %>% ScaleData()
qsave(seu, '/bigdata/zlin/Melanoma_meta/data/ssSTACAS_cd4t.qs')
print('done')

seu <- merge(x=datasets[[1]], y=datasets[2:length(datasets)], add.cell.ids = names(datasets)) %>% subset(subset = celltype_main %in% c('pDC','Mast','cDC','Mono','Macro'))
genes.blocklist <- c(GetSignature(SignatuR$Hs$Blocklists),
                     GetSignature(SignatuR$Hs$Compartments))
print('start integration')
seu <- seu %>% 
  FindVariableFeatures(nfeatures = nfeatures) %>%
  NormalizeData() %>% ScaleData() %>%
  # RunPCA(npcs=ndim) %>% RunUMAP(dims=1:ndim) %>% 
  SplitObject(split.by = "dataset") %>% 
  Run.STACAS(genesBlockList = genes.blocklist, cell.labels = "celltype_func",
             dims = 1:ndim, anchor.features = nfeatures) %>% RunUMAP(dims=1:ndim)
DefaultAssay(seu) <- 'RNA'
seu <- JoinLayers(seu)
seu <- seu %>% 
  FindVariableFeatures(nfeatures = nfeatures) %>%
  NormalizeData() %>% ScaleData()
qsave(seu, '/bigdata/zlin/Melanoma_meta/data/ssSTACAS_myeloids.qs')
print('done')
# 
# seu <- merge(x=datasets[[1]], y=datasets[2:length(datasets)], add.cell.ids = names(datasets)) %>% subset(subset = celltype_main %in% c('B','Plasma'))
# genes.blocklist <- c(GetSignature(SignatuR$Hs$Blocklists),
#                      GetSignature(SignatuR$Hs$Compartments))
# print('start integration')
# seu <- seu %>% 
#   FindVariableFeatures(nfeatures = nfeatures) %>%
#   NormalizeData() %>% ScaleData() %>%
#   # RunPCA(npcs=ndim) %>% RunUMAP(dims=1:ndim) %>% 
#   SplitObject(split.by = "dataset") %>% 
#   Run.STACAS(genesBlockList = genes.blocklist, cell.labels = "celltype_r2",
#              dims = 1:ndim, anchor.features = nfeatures) %>% RunUMAP(dims=1:ndim)
# DefaultAssay(seu) <- 'RNA'
# seu <- JoinLayers(seu)
# seu <- seu %>% 
#   FindVariableFeatures(nfeatures = nfeatures) %>%
#   NormalizeData() %>% ScaleData()
# qsave(seu, '/bigdata/zlin/Melanoma_meta/data/ssSTACAS_bplasma.qs')
# print('done')
# 
# seu <- merge(x=datasets[[1]], y=datasets[2:length(datasets)], add.cell.ids = names(datasets)) %>% subset(subset = celltype_main %in% c('Endo','CAF'))
# genes.blocklist <- c(GetSignature(SignatuR$Hs$Blocklists),
#                      GetSignature(SignatuR$Hs$Compartments))
# print('start integration')
# seu <- seu %>% 
#   FindVariableFeatures(nfeatures = nfeatures) %>%
#   NormalizeData() %>% ScaleData() %>%
#   # RunPCA(npcs=ndim) %>% RunUMAP(dims=1:ndim) %>% 
#   SplitObject(split.by = "dataset") %>% 
#   Run.STACAS(genesBlockList = genes.blocklist, cell.labels = "celltype_r2",
#              dims = 1:ndim, anchor.features = nfeatures) %>% RunUMAP(dims=1:ndim)
# DefaultAssay(seu) <- 'RNA'
# seu <- JoinLayers(seu)
# seu <- seu %>% 
#   FindVariableFeatures(nfeatures = nfeatures) %>%
#   NormalizeData() %>% ScaleData()
# qsave(seu, '/bigdata/zlin/Melanoma_meta/data/ssSTACAS_endocaf.qs')
# print('done')
# 
# seu <- qread('/bigdata/zlin/Melanoma_meta/data/ssSTACAS_cd4t.qs')
# CellDimPlot(seu, group.by = c("celltype_r2","dataset"), reduction = "umap", theme_use = "theme_blank")
# 
seu <- qread('/bigdata/zlin/Melanoma_meta/data/ssSTACAS_cd4t.qs')
CellDimPlot(seu, group.by = c("celltype_r2","celltype_func","dataset"), reduction = "umap", theme_use = "theme_blank")

library(harmony)
seu <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r2.qs') %>% subset(subset = celltype_main %in% c('pDC','Mast','cDC','Mono','Macro'))
seu <- seu %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% 
  RunPCA(verbose=FALSE) %>% 
  RunHarmony(group.by.vars = "dataset") %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) 
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.1)
CellDimPlot(seu, group.by = c("celltype_r2","dataset"), reduction = "umap", theme_use = "theme_blank")
CellDimPlot(seu, group.by = 'seurat_clusters', reduction = "umap", theme_use = "theme_blank")
seu <- AddModuleScore(seu,
                       features = list(GetSignature(SignatuR$Hs$Compartments$Ribo), GetSignature(SignatuR$Hs$Compartments$Mito)),
                       name=c("ribo","mito"))
# Plot scores
FeaturePlot(seu,
            features = c("C1QC","SPP1","CD14"))

marker_cosg <- cosg(
  seu,
  groups=c('0','1'),
  mu=100, ## If you would like to identify more specific marker genes, you could assign mu to larger values, such as mu=10 or mu=100
  n_genes_user=100)

x_BP = compareCluster(marker_cosg[[1]], fun='enrichGO', 
                      OrgDb = 'org.Hs.eg.db', keyType = 'SYMBOL', ont="BP")

p<-dotplot(x_BP) + theme(axis.text.x = element_text(angle=45, hjust=1),
                         axis.text.y = element_text(size=12),
                         panel.spacing = unit(5, "mm"))+
  scale_colour_gradientn(colours =c("#E54924","#498EA4"));p

