pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','COSG','BPCells','SeuratExtend','MetBrewer','ggplot2','CytoTRACE2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(max.print = 10000)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)  # Set to 2 GiB
setwd("/home/zlin/workspace/PanCancer_ICI")

# Malignant cells
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao',
              'HNSC_Franken',  
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 
              'PCa_Hawley', 'RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) 
  seu$celltype_main[seu$celltype_main == 'Malignant(CNA+)'] <- 'Epithelial(CNA+)'
  seu <- seu |>
    subset(subset = celltype_main %in% c('Melanocytes(CNA+)','Melanocytes(CNA-)','Epithelial(CNA+)','Epithelial(CNA-)')) |> 
    subset(subset = celltype_bped_main %in% c('B-cells', 'CD4+ T-cells', 'CD8+ T-cells', 'DC', 'Monocytes', 'Macrophages', 'Neutrophils', 'NK cells'), invert = T) |> 
    subset(subset = scGate_multi %in% c('Bcell', 'CD4T', 'CD8T', 'NK', 'Monocyte', 'Macrophage', 'PlasmaCell','panDC','Neutrophils','Mast'), invert = T)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu <- seu |> 
  subset(subset = sample %in% filter_sample) |> 
  NormalizeData() |> 
  FindVariableFeatures() |> 
  SketchData(ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
seu <- seu |> 
  FindVariableFeatures() |> 
  ScaleData() |>
  RunPCA(verbose=T) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony', verbose = T) |> 
  FindNeighbors(reduction = "harmony", dims = 1:30) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:30, reduction = 'harmony')
seu <- ProjectIntegration(seu, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
seu <- ProjectData(seu, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                   full.reduction = "harmony.full", dims = 1:30, refdata = list(seurat_clusters_full = 'seurat_clusters'))
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",
               reduction.key = "UMAP_full_")
qs_save(seu, 'data/seu_malignant.qs2')
seu <- qs_read('data/seu_malignant.qs2')
seu[["sketch"]] <- JoinLayers(seu[["sketch"]])
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis()
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 0.5, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)),
        label = T)

DimPlot(seu, group.by = 'celltype_main', reduction = 'umap.full', alpha = 0.5, 
        cols = rev(color_pro(length(unique(seu$celltype_main)), 1)),
        label = F) +
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1]) , y = min(seu@reductions$umap.full@cell.embeddings[,2]),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])) + 3, yend = min(seu@reductions$umap.full@cell.embeddings[,2])),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1]) , y = min(seu@reductions$umap.full@cell.embeddings[,2]),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])) , yend =min(seu@reductions$umap.full@cell.embeddings[,2]) + 3),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) + 1.5, y = min(seu@reductions$umap.full@cell.embeddings[,2]) -1, label = "UMAP_1",
           color="black", size = 3) + 
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) -1, y = min(seu@reductions$umap.full@cell.embeddings[,2]) + 1.5, label = "UMAP_2",
           color="black",size = 3, angle=90)
ggsave('figures/UMAP/UMAP_Malignant.png', height = 4, width = 5, dpi = 300)
DefaultAssay(seu) <- 'sketch'
seu <- JoinLayers(seu)
marker_genes <- FindAllMarkers(seu, min.pct = 0.25,logfc.threshold = 0.25, only.pos = T)
top30 <- marker_genes |> 
  group_by(cluster) |> 
  arrange(desc(avg_log2FC), .by_group=T) |> 
  filter(avg_log2FC>1, p_val_adj<0.05) |> 
  slice_head(n=30) |> 
  ungroup() 
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
gg <- bitr(top30$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db') 
top30 <- merge(top30, gg, by.x='gene', by.y = 'SYMBOL')

kk <- compareCluster(ENTREZID~cluster, data = top30, fun=enrichGO, OrgDb = org.Hs.eg.db, ont ='')
g <- dotplot(kk, label_format=100) + 
  aes(x=sub("\n.*", "", Cluster)) + 
  xlab("Cell Clusters") +
  ggtitle(NULL) +
  theme(axis.text.x = element_text(angle=30, hjust=1))
p3 <- p2 + coord_cartesian() + 
  ggfun::theme_noxaxis() + 
  xlab(NULL) 
insert_top(g, p3, height=.2)


# CD4+T cells
cd4_t <- c('CD4_Tcm','CD4_Treg','CD4_T-ISG','CD4_Tfh','CD4_Tstr','CD4_Tctl','CD4_Th17','CD4_T-naive')
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 'SCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'NSCLC_Liu',
              'PCa_Hawley', 'HCC_Guo', 'HCC_Ma','RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) |> 
    subset(subset = celltype_r2 %in% cd4_t) 
    # cytotrace2(ncores = 48, is_seurat = T, slot_type = "counts", species = "human", seed = 123)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu <- seu |> 
  NormalizeData() |> 
  FindVariableFeatures() |> 
  SketchData(ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
seu <- seu |> 
  FindVariableFeatures() |> 
  ScaleData() |>
  RunPCA(verbose=T) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony', verbose = T) |> 
  FindNeighbors(reduction = "harmony", dims = 1:30) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:30, reduction = 'harmony')
seu <- ProjectIntegration(seu, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
seu <- ProjectData(seu, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                   full.reduction = "harmony.full", dims = 1:30, refdata = list(seurat_clusters_full = 'seurat_clusters'))
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",
               reduction.key = "UMAP_full_")
qs_save(seu, 'data/seu_CD4T.qs2')
seu <- qs_read('data/seu_CD4T.qs2')
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 0.5, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)),
        label = T) 
seu <- JoinLayers(seu)
genes_to_check = list(c('CD3D','CD4','CD8A'), # T cells 'CD8B'
                      c('KLRD1','FCGR3A'), 
                      c('MKI67', 'TOP2A'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('MZB1','JCHAIN'),
                      c('KIT','TPSAB1'),
                      c('LILRA4','PLD4'),
                      c('CLEC9A','CD1C','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('PECAM1','VWF', 'ENG'),
                      c("RGS5",'ACTA2'),
                      c('COL1A1','FAP'),
                      c('KRT19', 'EPCAM'),
                      c('MLANA','TYR'),
                      c('RPL11','RPL10A')
)
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis()
seu$celltype_r2[seu$seurat_clusters_full == 6] <- 'Cycling'
seu$celltype_r2[seu$seurat_clusters_full %in% c(11,12,13,15)] <- 'Doublet'
seu$celltype_r2[seu$seurat_clusters_full == 14] <- 'RP-high'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$celltype_bped_main), c('CD4+ T-cells', 'CD8+ T-cells', 'NK cells','unknown', NA))] <- 'Doublet'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$scGate_multi), c('CD8T', 'CD4T', 'NK','unknown', NA))] <- 'Doublet'
adjusted_cd4t <- seu@meta.data |> 
  filter(celltype_r2 %in% c('Cycling','Doublet','RP-high')) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(adjusted_cd4t, 'tables/adjusted_cd4t.csv')
genes_to_check <- c('CD8A','CD3D','CD3E','CD4','TCF7','LEF1','CCR7','SELL',
                    'IL7R','CXCR4','CD55','GPR183','CD69','LMNA','ANXA1',
                    'DUSP2','NR4A2',"HSPA1A","HSPA6", 'FOS','JUN',
                    'TNF','GZMA','GZMB','GZMK','TBX21','CCL4','CCL5','IFNG','PRF1','CD40LG','STAT4','IL12RB2', #'CX3CR1','KLRG1','EOMES',IL2	CD40LG	TBX21	STAT4	IL12RB2
                    'CXCR5','IL21','BCL6','CXCL13',"GNG4", "CD200", "IGFL2",
                    'PDCD1','LAG3','HAVCR2','IL2RA','CTLA4','LAYN','RTKN2', 'FOXP3','TNFRSF9',
                    'ISG15',"IFIT1","IFIT3",'IFI44L',
                    'CXCR6','CCR6','KLRB1','RORA','RORC','IL17A','IL26','MKI67','TOP2A','RPL10A','RPL11')
DotPlot(seu, rev(genes_to_check), group.by = 'seurat_clusters_full', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  RotatedAxis() + coord_flip() + 
  xlab('') + ylab('') 
seu <- subset(seu, subset = celltype_r2 %in% c('Doublet','RP-high'), invert = T)
# rerun UMAP
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
seu$celltype_r2 <- factor(seu$celltype_r2, levels =  c('CD4_T-naive','CD4_Tcm','CD4_Tstr','CD4_Tctl','CD4_Tfh','CD4_Treg','CD4_T-ISG','CD4_Th17','Cycling'))
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=1, override.aes = list(size = 3)))
ggsave('figures/UMAP/UMAP_CD4+T.png', height = 4, width = 5.5, dpi = 300)

seu <- CreateSeuratObject(counts = GetAssayData(seu, layer = 'counts'), meta.data = seu@meta.data) |> NormalizeData()
seu$celltype_r2 <- factor(seu$celltype_r2, levels =  c('CD4_T-naive','CD4_Tcm','CD4_Tctl','CD4_Tfh','CD4_Treg','CD4_Th17','CD4_T-ISG','CD4_Tstr','Cycling'))
genes_to_check <- c('TCF7','LEF1','CCR7','SELL',
                    'IL7R','CXCR4','CD55','GPR183','CD69','LMNA','ANXA1',
                    'DUSP2','NR4A2',
                    'TNF','GZMA','GZMB','GZMK','TBX21','CCL4','CCL5','IFNG','PRF1','CX3CR1','KLRG1','EOMES', #IL2	CD40LG	TBX21	STAT4	IL12RB2
                    'CD40LG','STAT4','CXCR5','IL21','BCL6','CXCL13',"GNG4", "CD200", "IGFL2",
                    'PDCD1','LAG3','HAVCR2','IL2RA','CTLA4','LAYN','RTKN2', 'FOXP3','TNFRSF9',
                    'CXCR6','CCR6','KLRB1','RORA','RORC','IL17A','IL26',
                    'ISG15',"IFIT1","IFIT3",'IFI44L',
                    "HSPA1B","HSPA6", 'FOS','JUN','MKI67','TOP2A')
DotPlot(seu, rev(genes_to_check), group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() + coord_flip() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9)) +
  guides(color = guide_colorbar(frame.colour = "black",
                                frame.linewidth = 0.5,
                                ticks.colour = "black",
                                ticks.linewidth = 0.1,
                                barwidth = 5,
                                barheight = 1,
                                title = "Expression \n(scaled)",
                                title.position = "top",
                                title.hjust = 0.5,
                                direction = 'horizontal')) +
  scale_color_gradientn(colors = rev(brewer.pal(8, "RdBu")),  
                        breaks = c(-1, 1),  
                        labels = c("Low", "High"), 
                        limits = c(-1, 1)) + 
  xlab('') + ylab('') 
ggsave('figures/Dotplot/CD4T.pdf', height = 10, width = 5)
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_bplasma.csv', row.names = F)

# CD8+T cells
cd8_t <- c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',"CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',"CD8_NK-like",'gdT','MAIT')
nk <- c('NK_CD56hiCD16lo', 'NK_CD56loCD16hi')
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 'SCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'NSCLC_Liu',
              'PCa_Hawley','HCC_Guo','HCC_Ma','RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) |> 
    subset(subset = celltype_r2 %in% c(cd8_t,nk)) 
    # cytotrace2(ncores = 48, is_seurat = T, slot_type = "counts", species = "human", seed = 123)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu <- seu |> 
  NormalizeData() |> 
  FindVariableFeatures() |> 
  SketchData(ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
seu <- seu |> 
  FindVariableFeatures() |> 
  ScaleData() |>
  RunPCA(verbose=T) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony', verbose = T) |> 
  FindNeighbors(reduction = "harmony", dims = 1:30) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:30, reduction = 'harmony')
seu <- ProjectIntegration(seu, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
seu <- ProjectData(seu, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                   full.reduction = "harmony.full", dims = 1:30, refdata = list(seurat_clusters_full = 'seurat_clusters'))
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",
               reduction.key = "UMAP_full_")
qs_save(seu, 'data/seu_CD8TNK.qs2')
seu <- qs_read('data/seu_CD8TNK.qs2')
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 0.5, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)),
        label = T) 
seu <- JoinLayers(seu)
genes_to_check = list(c('CD3D','CD4','CD8A'), # T cells 'CD8B'
                      c('KLRD1','FCGR3A'), 
                      c('MKI67', 'TOP2A'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('MZB1','JCHAIN'),
                      c('KIT','TPSAB1'),
                      c('LILRA4','PLD4'),
                      c('CLEC9A','CD1C','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('PECAM1','VWF', 'ENG'),
                      c("RGS5",'ACTA2'),
                      c('COL1A1','FAP'),
                      c('KRT19', 'EPCAM'),
                      c('MLANA','TYR'),
                      c('RPL11','RPL10A')
)
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis()
genes_to_check = list(c('CD3D','CD4','CD8A'), # T cells 'CD8B'
                      c('KLRD1','FCGR3A'), 
                      c('MKI67', 'TOP2A'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1'),  # B cells 
                      c('MZB1','JCHAIN'),
                      c('KIT','TPSAB1'),
                      c('LILRA4','PLD4'),
                      c('CLEC9A','CD1C','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('PECAM1','VWF', 'ENG'),
                      c("RGS5",'ACTA2'),
                      c('COL1A1','FAP'),
                      c('KRT19', 'EPCAM'),
                      c('MLANA','TYR'),
                      c('RPL11','RPL10A')
)
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters_full', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  RotatedAxis()
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)), label = T)
seu$celltype_r2[seu$seurat_clusters_full == 5] <- 'Cycling'
seu$celltype_r2[seu$seurat_clusters_full %in% c(9,10)] <- 'Doublet'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$celltype_bped_main), c('CD4+ T-cells', 'CD8+ T-cells', 'NK cells','unknown', NA))] <- 'Doublet'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$scGate_multi), c('CD8T', 'CD4T', 'NK','unknown', NA))] <- 'Doublet'
adjusted_cd8tnk <- seu@meta.data |> 
  filter(celltype_r2 %in% c('Cycling','Doublet')) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(adjusted_cd8tnk, 'tables/adjusted_cd8tnk.csv')
seu <- subset(seu, subset = celltype_r2 == 'Doublet', invert = T)
# rerun UMAP
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
seu$celltype_r2 <- factor(seu$celltype_r2, levels = c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
                          "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
                          "CD8_NK-like",'MAIT','gdT','NK_CD56loCD16hi','NK_CD56hiCD16lo','Cycling'))
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=1, override.aes = list(size = 3)))
ggsave('figures/UMAP/UMAP_CD8+T.png', height = 4, width = 5.5, plot = p)

genes_to_check <- c('CD3D','CD8A','CD8B','CD4', 'TCF7','LEF1','CCR7','SELL',
                    'IL7R', 'ANXA1', 'CD55', 'CD27','LMNA','ZNF683','CXCR6','ITGAE','CXCR5',
                    'GZMB','GZMH','GZMK','EOMES','CXCR3','IFNG','PRF1','CXCL13',
                    'CTLA4','PDCD1','LAG3','TIGIT','LAYN','ISG15','IFIT1','IFIT3', "HSPA1A","HSPA6","NR4A1",# 'MYL12A','MYL12B',
                    'SLC4A10','KLRB1','RORA',
                    'TRDV2','TRGV9','CX3CR1','TBX21','FGFBP2','KLRF1','FCGR3A',
                    'GNLY','TYROBP','XCL1','XCL2','NCAM1','MKI67','TOP2A','RPL11','RPL10A','CD19','CD79A')
DotPlot(seu, rev(genes_to_check), group.by = 'seurat_clusters_full', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() + coord_flip() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9)) +
  guides(color = guide_colorbar(frame.colour = "black",
                                frame.linewidth = 0.5,
                                ticks.colour = "black",
                                ticks.linewidth = 0.1,
                                barwidth = 5,
                                barheight = 1.1,
                                title = "Expression \n(scaled)",
                                title.position = "top",
                                title.hjust = 0.5,
                                direction = 'horizontal'),
         size = guide_legend(direction = 'horizontal',
                             title = "Percentage",
                             title.position = "top",
                             title.hjust = 0.5)) +
  scale_color_gradientn(colors = rev(brewer.pal(8, "RdBu")),  
                        breaks = c(-1, 1),  
                        labels = c("Low", "High"), 
                        limits = c(-1, 1)) + 
  xlab('') + ylab('')
DimPlot(seu, group.by = 'celltype_main', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$cohort)), 1)),
        label = T)
excluded_cd8t <- seu@meta.data |> 
  filter(seurat_clusters_full %in% c(7,8,12)) |> 
  mutate(celltype_r2 = case_when(seurat_clusters_full %in% c(7,8) ~ 'Cycling',
                                 seurat_clusters_full == 12 ~ 'Doublet')) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(excluded_cd8t, 'tables/excluded_cd8t.csv')
seu <- subset(seu, subset = seurat_clusters_full %in% c(7,8,12), invert = T)
seu$celltype_r2 <- factor(seu$celltype_r2, levels =  c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',"CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',"CD8_NK-like"))
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)),
        label = F) +
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2,
                              override.aes = list(size = 3))) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1]) , y = min(seu@reductions$umap.full@cell.embeddings[,2]),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])) + 3, yend = min(seu@reductions$umap.full@cell.embeddings[,2])),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1]) , y = min(seu@reductions$umap.full@cell.embeddings[,2]),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])) , yend =min(seu@reductions$umap.full@cell.embeddings[,2]) + 3),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) + 1.5, y = min(seu@reductions$umap.full@cell.embeddings[,2]) -1, label = "UMAP_1",
           color="black", size = 3) + 
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) -1, y = min(seu@reductions$umap.full@cell.embeddings[,2]) + 1.5, label = "UMAP_2",
           color="black",size = 3, angle=90) +
  guides(color = guide_legend(ncol=2,
                              override.aes = list(size = 3)))
ggsave('figures/UMAP/UMAP_CD8T.png', height = 4, width = 7, dpi = 300)
seu[["sketch"]] <- JoinLayers(seu[["sketch"]])
genes_to_check <- c('CD3D','CD8A','CD8B',
                    'TCF7','LEF1','CCR7','SELL',
                    'IL7R', 'ANXA1', 'CD55', 'CD27','LMNA','ZNF683','CXCR6','ITGAE','CXCR5',
                    'GZMB','GZMH','GZMK','EOMES','CXCR3','IFNG','PRF1','CXCL13',
                    'CTLA4','PDCD1','LAG3','TIGIT','LAYN','ISG15','IFIT1','IFIT3', "HSPA1A","HSPA6","NR4A1",# 'MYL12A','MYL12B',
                    'SLC4A10','KLRB1','RORA',
                    'TRDV2','TRGV9','CX3CR1','TBX21','FGFBP2','KLRF1','FCGR3A',
                    'GNLY','TYROBP','XCL1','XCL2','NCAM1','MKI67','TOP2A')
DotPlot(seu, rev(genes_to_check), group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() + coord_flip() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9)) +
  guides(color = guide_colorbar(frame.colour = "black",
                                frame.linewidth = 0.5,
                                ticks.colour = "black",
                                ticks.linewidth = 0.1,
                                barwidth = 5,
                                barheight = 1.1,
                                title = "Expression \n(scaled)",
                                title.position = "top",
                                title.hjust = 0.5,
                                direction = 'horizontal'),
         size = guide_legend(direction = 'horizontal',
                             title = "Percentage",
                             title.position = "top",
                             title.hjust = 0.5)) +
  scale_color_gradientn(colors = rev(brewer.pal(8, "RdBu")),  
                        breaks = c(-1, 1),  
                        labels = c("Low", "High"), 
                        limits = c(-1, 1)) + 
  xlab('') + ylab('')
ggsave('figures/Dotplot/CD8TNK.pdf', height = 12, width = 7)

# B/plasma cells
bplasma <- c("B-naive", "B-ISG", "B-HSP", "B_MT2A", 
             "ACB_EGR1", "ACB_NR4A2", "ACB_CCR7", "B-memory", "B-AtM", 
             "GCB-pre", "GCB-DZ_SUGCT", "GCB-LZ_LMO2",
             "GCB-cycling", "PC-cycling",
             "PC-early_RGS13", "PC-early_LTB", "PC_IGHG", "PC_IGHA")
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 
              'PCa_Hawley','HCC_Guo','HCC_Ma','RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) |> 
    subset(subset = celltype_r2 %in% bplasma) 
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu <- seu |> 
  NormalizeData() |> 
  FindVariableFeatures() |> 
  SketchData(ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
seu <- seu |> 
  FindVariableFeatures() |> 
  ScaleData() |>
  RunPCA(verbose=T) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony', verbose = T) |> 
  FindNeighbors(reduction = "harmony", dims = 1:30) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:30, reduction = 'harmony')
seu <- ProjectIntegration(seu, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
seu <- ProjectData(seu, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                   full.reduction = "harmony.full", dims = 1:30, refdata = list(seurat_clusters_full = 'seurat_clusters'))
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",
               reduction.key = "UMAP_full_")
qs_save(seu, 'data/seu_Bplasma.qs2')
seu <- qs_read('data/seu_Bplasma.qs2')
seu[["sketch"]] <- JoinLayers(seu[["sketch"]])
genes_to_check = list(c('CD3D','CD4','CD8A'), # T cells 'CD8B'
                      c('KLRD1','FCGR3A'), 
                      c('MKI67', 'TOP2A'), # NK cells 'KLRB1', 'KLRC1', 'CD16', 'CD56', 'CD11b', 'CD11c'
                      c('CD79A','CD19', 'MS4A1','IGHG1','IGHA1'),  # B cells 
                      c('MZB1','JCHAIN'),
                      c('KIT','TPSAB1'),
                      c('LILRA4','PLD4'),
                      c('CLEC9A','CD1C','LAMP3'), 
                      c('CD68', 'LYZ', 'CD14'),  
                      c('PECAM1','VWF', 'ENG'),
                      c("RGS5",'ACTA2'),
                      c('COL1A1','FAP'),
                      c('KRT19', 'EPCAM'),
                      c('MLANA','TYR'),
                      c('RPL11','RPL10A')
)
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  RotatedAxis()
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)), label = T)
seu$celltype_r2[seu$seurat_clusters_full != 11 & seu$celltype_r2 == "PC-cycling"] <- "PC-trans"
seu$celltype_r2[seu$seurat_clusters_full == 0 & seu$celltype_r2 %in% c("PC-cycling", "PC-early_RGS13", "PC-early_LTB", "PC_IGHG", "PC_IGHA")] <- "Mixed"
seu$celltype_r2[seu$seurat_clusters_full %in% c(8,9,12,15)] <- "Doublet"
adjusted_bplasma <- seu@meta.data |> 
  filter(celltype_r2 %in% c("Mixed", "PC-trans", "Doublet")) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(adjusted_bplasma, 'tables/adjusted_bplasma.csv')
seu <- subset(seu, subset = seurat_clusters_full %in% c(8,9,12,15), invert = T)
seu <- subset(seu, subset = celltype_r2 == 'Mixed', invert = T)
# rerun UMAP
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
seu$celltype_r2[seu$celltype_r2 == 'PC-early_LTB'] <- "PC_IGHG"
seu$celltype_r2 <- factor(seu$celltype_r2, levels =  c("B-naive", "B-ISG", "B-HSP", "B_MT2A", 
                                                       "ACB_EGR1", "ACB_NR4A2", "ACB_CCR7", "B-memory", "B-AtM", 
                                                       "GCB-pre", "GCB-DZ_SUGCT", "GCB-LZ_LMO2",
                                                       "GCB-cycling", "PC-cycling","PC-trans",
                                                       "PC-early_RGS13", "PC_IGHG", "PC_IGHA"))
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
ggsave('figures/UMAP/UMAP_Bplasma.png', height = 4, width = 6, dpi = 300)

seu <- CreateSeuratObject(counts = GetAssayData(seu, layer = 'counts'), meta.data = seu@meta.data) |> NormalizeData()
genes_to_check <- c('TCL1A','YBX3','FCER2',
                    'IFIT3', 'IFI44L', 'STAT1', 'ISG15',
                    'HSPA1A','DNAJB1','HSPA1B',
                    'MT1X','MT2A',
                    'TNF','EGR1',
                    'ZNF331','VPS37B','NR4A2','CD69','CD86',
                    'CD83','CCR7','NFKBID',
                    'CRIP1','S100A10','S100A4', 'ITGB1',
                    'DUSP4','FCRL4','CCR1','ZEB2','ITGAX',
                    'ENO1','PSME2','NME1',
                    'SUGCT','NEIL1','MEF2B','RGS13',
                    'LMO2','GMDS','MARCKSL1',
                    'MKI67','TOP2A',
                    'XBP1','JCHAIN','SSR4',
                    "CD38", "MZB1", "PRDM1", "IRF4", 
                    'CD19','CD27',"MS4A1","LTB", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DQA1",
                    'IGHM','IGHD',"IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2")
DotPlot(seu, genes_to_check, group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() +# coord_flip() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9)) +
  guides(color = guide_colorbar(frame.colour = "black",
                                frame.linewidth = 0.5,
                                ticks.colour = "black",
                                ticks.linewidth = 0.1,
                                barwidth = 5,
                                barheight = 0.9,
                                title = "Expression \n(scaled)",
                                title.position = "top",
                                title.hjust = 0.5,
                                direction = 'horizontal'),
         size = guide_legend(direction = 'vertical',
                             title = "Percentage",
                             title.position = "top",
                             title.hjust = 0.5)) +
  scale_color_gradientn(colors = rev(brewer.pal(8, "RdBu")),  
                        breaks = c(-1, 1),  
                        labels = c("Low", "High"), 
                        limits = c(-1, 1)) + 
  xlab('') + ylab('') 
ggsave('figures/Dotplot/Bplasma.pdf', height = 4.5, width = 14)
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_bplasma.csv', row.names = F)

# Myeloids
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 
              'PCa_Hawley','HCC_Guo','HCC_Ma','RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) |> 
    subset(subset = celltype_main %in% c('Mast','pDC','cDC','Mono/macro')) 
  seu$cell.id <- colnames(seu)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu <- seu |> 
  NormalizeData() |> 
  FindVariableFeatures() |> 
  SketchData(ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")
seu <- seu |> 
  FindVariableFeatures() |> 
  ScaleData() |>
  RunPCA(verbose=T) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony', verbose = T) |> 
  FindNeighbors(reduction = "harmony", dims = 1:30) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:30, reduction = 'harmony')
seu <- ProjectIntegration(seu, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
seu <- ProjectData(seu, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                   full.reduction = "harmony.full", dims = 1:30, refdata = list(seurat_clusters_full = 'seurat_clusters'))
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",
               reduction.key = "UMAP_full_")
qs_save(seu, 'data/seu_Myeloids.qs2')
seu <- qs_read('data/seu_Myeloids.qs2')
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 0.5, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)),
        label = T)
seu[["sketch"]] <- JoinLayers(seu[["sketch"]])
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') +
  RotatedAxis()
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)), label = T)
seu$cell.id <- paste0(seu$cohort,'-',seu$cell.id)
adjusted_myeloids <- seu@meta.data |> 
  filter(seurat_clusters_full %in% c(8,9,14,15)) |> 
  mutate(celltype_r2 = case_when(seurat_clusters_full == 9 ~ 'Cycling',
                                 seurat_clusters_full %in% c(8,14,15) ~ 'Doublet')) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(adjusted_myeloids, 'tables/adjusted_myeloids.csv')
seu <- subset(seu, subset = seurat_clusters_full %in% c(8,9,14,15), invert = T)
# rerun UMAP
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
myeloids <- c('Mast','pDC','cDC1', 
              'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'mregDC', 'MoDC', 
              'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
              'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro-ISG', 
              'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2')
seu$celltype_r2 <- factor(seu$celltype_r2, levels = myeloids)
p <- DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
qs_save(p, 'data/umap_myeloids.qs2')
p <- qs_read('data/umap_myeloids.qs2')
p
ggsave('figures/UMAP/UMAP_myeloids.png', height = 4, width = 7, dpi = 300)



genes_to_check <- c('TPSAB1','LILRA4','SELL',
                    'CLEC9A','XCR1',
                    'CD1C','FCER1A','IL1B','NLRP3','NFKBIA','ICAM1',
                    'CD1A','CD207','LAMP3','CCR7','IDO1','NFKB1','DUSP4','MARCKSL1','BIRC3',
                    'S100A8','S100A9','FCN1','CD14','CD300E','FCGR3A','LST1','LILRB2',
                    'CXCL2','TREM1','AREG','EREG',
                    'VCAN','VEGFA','INHBA','SPP1','FN1','OLR1',
                    'ISG15','IFIT1','IFIT2','IFI6','CXCL9','CXCL10','TNFRSF14','IL4I1',
                    'CCL2','CCL3','CCL4','DUSP1','TNF','JUNB','ATF3','CCL3L1',
                    'NR4A2','LYVE1','FOLR2','MRC1','MSR1',
                    'C1QC',"MARCKS",'APOE','TREM2','GPNMB','MMP9','CD36','HLA-DQA1','CD274','PDCD1LG2','TOP2A','MKI67','RPL10A','RPL11') # 'TGFB1','CD80','CD86'
DotPlot(seu, rev(genes_to_check), group.by = 'seurat_clusters_full', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() + coord_flip() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9)) +
  guides(color = guide_colorbar(frame.colour = "black",
                                frame.linewidth = 0.5,
                                ticks.colour = "black",
                                ticks.linewidth = 0.1,
                                barwidth = 5,
                                barheight = 1.1,
                                title = "Expression \n(scaled)",
                                title.position = "top",
                                title.hjust = 0.5,
                                direction = 'horizontal')) +
  scale_color_gradientn(colors = rev(brewer.pal(8, "RdBu")),  
                        breaks = c(-1, 1),  
                        labels = c("Low", "High"), 
                        limits = c(-1, 1)) + 
  xlab('') + ylab('') 
seu$celltype_r2[seu$seurat_clusters_full == 9] <- "Cycling"
mye_cycling <- seu@meta.data |> select(cell.id, celltype_r2) |> filter(celltype_r2 == 'Cycling')
write.csv(mye_cycling, 'tables/mye_cycling.csv')
myeloids <- c('Mast','pDC','cDC1', 
              'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'mregDC', 'MoDC', 
              'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
              'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro-ISG', 
              'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2','Cycling')
seu$celltype_r2 <- factor(seu$celltype_r2, levels = myeloids)
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 0.5, 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)),
        label = F) +
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1]-1) , y = min(seu@reductions$umap.full@cell.embeddings[,2]-1),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])) + 2, yend = min(seu@reductions$umap.full@cell.embeddings[,2])-1),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1])-1, y = min(seu@reductions$umap.full@cell.embeddings[,2]-1),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])-1) , yend =min(seu@reductions$umap.full@cell.embeddings[,2]) + 2),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) + 0.5, y = min(seu@reductions$umap.full@cell.embeddings[,2]) -2, label = "UMAP_1",
           color="black", size = 3) + 
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) -2, y = min(seu@reductions$umap.full@cell.embeddings[,2]) + 0.5, label = "UMAP_2",
           color="black",size = 3, angle=90)
ggsave('figures/UMAP/UMAP_Myeloids.png', height = 4, width = 7, dpi = 300)

seu[["sketch"]] <- JoinLayers(seu[["sketch"]])
genes_to_check <- c('TPSAB1','LILRA4','SELL',
                    'CLEC9A','XCR1',
                    'CD1C','FCER1A','IL1B','NLRP3','NFKBIA','ICAM1',
                    'CD1A','CD207','LAMP3','CCR7','IDO1','NFKB1','DUSP4','MARCKSL1','BIRC3',
                    'S100A8','S100A9','FCN1','CD14','CD300E','FCGR3A','LST1','LILRB2',
                    'CXCL2','TREM1','AREG','EREG',
                    'VCAN','VEGFA','INHBA','SPP1','FN1','OLR1',
                    'ISG15','IFIT1','IFIT2','IFI6','CXCL9','CXCL10','TNFRSF14','IL4I1',
                    'CCL2','CCL3','CCL4','DUSP1','TNF','JUNB','ATF3','CCL3L1',
                    'NR4A2','LYVE1','FOLR2','MRC1','MSR1',
                    "MARCKS",'C1QC','APOE','TREM2','GPNMB','MMP9','CD36','TOP2A','MKI67','HLA-DQA1') # 'TGFB1','CD80','CD86'
DotPlot(seu, rev(genes_to_check), group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() + coord_flip() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9)) +
  guides(color = guide_colorbar(frame.colour = "black",
                                frame.linewidth = 0.5,
                                ticks.colour = "black",
                                ticks.linewidth = 0.1,
                                barwidth = 5,
                                barheight = 1.1,
                                title = "Expression \n(scaled)",
                                title.position = "top",
                                title.hjust = 0.5,
                                direction = 'horizontal')) +
  scale_color_gradientn(colors = rev(brewer.pal(8, "RdBu")),  
                        breaks = c(-1, 1),  
                        labels = c("Low", "High"), 
                        limits = c(-1, 1)) + 
  xlab('') + ylab('') 
ggsave('figures/Dotplot/Myeloids.pdf', height = 12, width = 7)
# DCs
seu_sub <- subset(seu, subset = celltype_r2 %in% c('pDC','cDC1', 'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'mregDC', 'MoDC'))
seu_sub$celltype_r2[seu_sub$celltype_r2 %in% c('cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'MoDC')] <- 'cDC2(CD1C+)'
seu_sub$celltype_r2 <- factor(seu_sub$celltype_r2, levels = c('pDC','cDC1','cDC2(CD1C+)','mregDC'))
genes_to_check <- list(
  Antigen = c("HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "CD74"), 
  Costimulation = c("CD40", "CD80", "CD83", "CD70", "CD86", "RELB", "TNFSF4", "ICOSLG", "TNFSF9", "TNFSF18"),     
  Immune_checkpoints = c("CD274", "PDCD1LG2", "PVR", "LGALS9", "IDO1", "FAS", "CD200", "HAVCR2", "LILRB1"),  
  Soluble_factors = c("IL10", "IL12B", "CCL2", "CCL4", "CCL5", "XCL1", "CXCL9", "CXCL10", "IFI6", "ISG15"),  
  TLRs_and_adaptors = c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6","TLR7", "TLR8", "TLR9", "MYD88", "TICAM1"),  
  Migration = c("ICAM1", "MARCKS", "MARCKSL1", "MYO1G", "CCR7")       
)

DotPlot(seu_sub, c("CD40", "CD80", "CD83", "CD70", "CD86", "RELB"), 
        group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  coord_flip() + RotatedAxis() + ggtitle('Maturation') + ylab('') + xlab('') + guides(color="none", size = 'none') + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        plot.title = element_text(hjust = 0.5, face = 'plain'))
ggsave('figures/Functional_score_Mac/cDC_maturation.pdf', height = 4, width = 2)
DotPlot(seu_sub, c("CD274", "PDCD1LG2","FAS", "CD200","SOCS1","SOCS2"), 
        group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') +
  theme_minimal() + 
  coord_flip() + RotatedAxis() + ggtitle('Regulatory') + ylab('') + xlab('') + guides(color="none", size = 'none') + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        plot.title = element_text(hjust = 0.5, face = 'plain'))
ggsave('figures/Functional_score_Mac/cDC_regulatory.pdf', height = 4, width = 2.5)
DotPlot(seu_sub, c("ICAM1", "MARCKS", "MARCKSL1", "MYO1G", "CCR7","CXCL16","FSCN1"), 
        group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  coord_flip() + RotatedAxis() + ggtitle('Migration') + ylab('') + xlab('') + guides(color="none", size = 'none') + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        plot.title = element_text(hjust = 0.5, face = 'plain'))
ggsave('figures/Functional_score_Mac/cDC_migration.pdf', height = 4, width = 2.5)
DotPlot(seu_sub, c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6","TLR7", "TLR8", "TLR9", "MYD88", "TICAM1","MAVS"), 
        group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  coord_flip() + RotatedAxis() + ggtitle('TRLs') + ylab('') + xlab('') + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 10),
        plot.title = element_text(hjust = 0.5, face = 'plain')) +
  guides(color = guide_colorbar(frame.colour = "black",
                                frame.linewidth = 0.5,
                                ticks.colour = "black",
                                ticks.linewidth = 0.1,
                                barwidth = 5,
                                barheight = 0.8,
                                title = "Expression \n(scaled)",
                                title.position = "top",
                                title.hjust = 0.5,
                                direction = 'horizontal')) +
  scale_color_gradientn(colors = rev(brewer.pal(8, "RdBu")),  
                        breaks = c(-1, 1),  
                        labels = c("Low", "High"), 
                        limits = c(-1, 1))
ggsave('figures/Functional_score_Mac/cDC_tlrs.pdf', height = 4, width = 4.2)


non_immune <- metadata |> filter(celltype_main %in% c('Endo','CAF','Mural')) |> pull(celltype_r2) |> unique()
non_immune <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein","Pericytes","SMC", "Myofibroblasts", "CAF_SFRP2", 
                "CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap")
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao',
              'HNSC_Franken',  
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 
              'PCa_Hawley','RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) |> 
    subset(subset = celltype_r2 %in% non_immune) 
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu <- seu |> 
  subset(subset = sample %in% filter_sample) |> 
  NormalizeData() |> 
  FindVariableFeatures() |> 
  SketchData(ncells = 3000, method = "LeverageScore", sketched.assay = "sketch")
seu <- seu |> 
  FindVariableFeatures() |> 
  ScaleData() |>
  RunPCA(verbose=T) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony', verbose = T) |> 
  FindNeighbors(reduction = "harmony", dims = 1:30) |>
  FindClusters(resolution = 0.5) |> 
  RunUMAP(dims = 1:30, reduction = 'harmony')
seu <- ProjectIntegration(seu, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
seu <- ProjectData(seu, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                   full.reduction = "harmony.full", dims = 1:30, refdata = list(seurat_clusters_full = 'seurat_clusters'))
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full",
               reduction.key = "UMAP_full_")
qs_save(seu, 'data/seu_nonimmune.qs2')
seu <- qs_read('data/seu_nonimmune.qs2')
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 0.5, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)),
        label = T)
seu[["sketch"]] <- JoinLayers(seu[["sketch"]])
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') +
  RotatedAxis()
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)), label = T)
seu$cell.id <- paste0(seu$cohort,'-',seu$cell.id)
adjusted_nonimmune <- seu@meta.data |> 
  filter(seurat_clusters_full %in% c(7,8,9,10,11,15)) |> 
  mutate(celltype_r2 = case_when(seurat_clusters_full == 9 ~ 'Cycling',
                                 seurat_clusters_full %in% c(7,10,11,15) ~ 'Doublet',
                                 seurat_clusters_full == 8 ~ 'low-quality')) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(adjusted_nonimmune, 'tables/adjusted_nonimmune.csv')
seu <- subset(seu, subset = seurat_clusters_full %in% c(7,8,9,10,11,15), invert = T)
# rerun UMAP
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
non_immune <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein","Pericytes","SMC", "Myofibroblasts", "CAF_SFRP2", 
                "CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap")
seu$celltype_r2 <- factor(seu$celltype_r2, levels = non_immune)
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
ggsave('figures/UMAP/UMAP_nonimmune.png', height = 4, width = 7, dpi = 300)



non_immune <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein","Pericytes","SMC", "Myofibroblasts", "CAF_SFRP2", 
                "CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap")
seu$celltype_r2 <- factor(seu$celltype_r2, levels = non_immune)
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)),
        label = F) +
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1]) , y = min(seu@reductions$umap.full@cell.embeddings[,2]),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])) + 3, yend = min(seu@reductions$umap.full@cell.embeddings[,2])),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  geom_segment(aes(x = min(seu@reductions$umap.full@cell.embeddings[,1]) , y = min(seu@reductions$umap.full@cell.embeddings[,2]),
                   xend = min(min(seu@reductions$umap.full@cell.embeddings[,1])) , yend =min(seu@reductions$umap.full@cell.embeddings[,2]) + 3),
               colour = "black", size = 0.5,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) + 1.5, y = min(seu@reductions$umap.full@cell.embeddings[,2]) -1, label = "UMAP_1",
           color="black", size = 3) + 
  annotate("text", x = min(seu@reductions$umap.full@cell.embeddings[,1]) -1, y = min(seu@reductions$umap.full@cell.embeddings[,2]) + 1.5, label = "UMAP_2",
           color="black",size = 3, angle=90) +
  guides(color = guide_legend(ncol=1,
                              override.aes = list(size = 3)))
ggsave('figures/UMAP/UMAP_nonimmune.png', height = 4, width = 5.5)

seu[["sketch"]] <- JoinLayers(seu[["sketch"]])
genes_to_check <- c('PROX1', 'LYVE1','FLT4','CCL21','COL9A3',
                    'GJA5','FBLN5','GJA4',
                    'CA4','CD36','RGCC',
                    'COL4A1','KDR','ESM1','CXCR4',
                    'ACKR1','SELP','CLU',
                    'RGS5','ACTA2','MYH11',
                    'HOPX','TGFB1','HMGA2','HMGA1','CDKN2A','CDH2',
                    'FAP','SFRP2','SFRP4','IGF1',
                    'PI16','CD34','CD55','MFAP5',
                    'LRRC15','POSTN',
                    'WNT5A','GREM1',
                    'MMP1','MMP3','CXCL13','ISG15','IL7R','IL6','CXCL1','CXCL2','CEBPD','NFKB1',
                    'C7','ADAMDEC1','CD74','B2M') 
DotPlot(seu, rev(genes_to_check), group.by = 'seurat_clusters_full', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() + coord_flip() + 
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        # axis.text.x = element_text(color = "black", size = 9, hjust = 0),
        axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9)) +
  guides(color = guide_colorbar(frame.colour = "black",
                                frame.linewidth = 0.5,
                                ticks.colour = "black",
                                ticks.linewidth = 0.1,
                                barwidth = 5,
                                barheight = 0.9,
                                title = "Expression \n(scaled)",
                                title.position = "top",
                                title.hjust = 0.5,
                                direction = 'horizontal'),
         size = guide_legend(direction = 'horizontal',
                             title = "Percentage",
                             title.position = "top",
                             title.hjust = 0.5)) +
  scale_color_gradientn(colors = rev(brewer.pal(8, "RdBu")),  
                        breaks = c(-1, 1),  
                        labels = c("Low", "High"), 
                        limits = c(-1, 1)) + 
  xlab('') + ylab('') 
ggsave('figures/Dotplot/Non-immune.pdf', height = 10, width = 6)



# Bhattacharrya distance
library(distdimscr)
library(tidyverse)
DefaultAssay(seu) <- 'RNA'
celltype_main <- c("CD4+T", "CD8+T", "NK", "B", "Plasma", 
                   "Mast", "pDC", "cDC", "Mono/macro", 
                   "Endo", "Mural", "CAF")
bhatt.dist <- bhatt.dist.rand <- as.data.frame(matrix(NA,ncol=length(celltype_main),nrow = 100))
names(bhatt.dist)<-celltype_main
names(bhatt.dist.rand)<-celltype_main

for (CT in celltype_main){
  print(CT)
  for (j in 1:10){
    print(j)
    set.seed(j)
    n.pre <- length(which(seu$time_point=="Pre" & seu$celltype_main==CT))
    n.on <- length(which(seu$time_point=="On" & seu$celltype_main==CT))
    if (n.pre < n.on){
      cells.on <-sample(colnames(seu)[which(seu$time_point=="On" & seu$celltype_main==CT)], n.pre)
      cells.pre <- colnames(seu)[seu$celltype_main==CT & seu$time_point=="Pre"]
    } else {
      cells.pre <-sample(colnames(seu)[which(seu$time_point=="Pre" & seu$celltype_main==CT)], n.on)
      cells.on <- colnames(seu)[seu$celltype_main==CT & seu$time_point=="On"]
    }
    # Harmony embeddings
    tmp<-seu@reductions$harmony.full@cell.embeddings
    cells.pre.embedding <- tmp[cells.pre,]
    cells.on.embedding <- tmp[cells.on,]
    for (i in 1:10) {
      d <- (j-1)*10+i
      bhatt.dist[d,CT] <- dim_dist(embed_mat_x=cells.pre.embedding,embed_mat_y=cells.on.embedding,dims_use=1:20,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=FALSE)
      bhatt.dist.rand[d,CT] <- dim_dist(embed_mat_x=cells.pre.embedding,embed_mat_y=cells.on.embedding,dims_use=1:20,num_cells_sample=100,distance_metric="bhatt_dist",random_sample=TRUE)
    }
  }
}
bhatt.dist.relative <- bhatt.dist - bhatt.dist.rand
write.csv(bhatt.dist.relative, 'tables/bhatt_dist_main_relative.csv', row.names = F)
write.csv(bhatt.dist, 'tables/bhatt_dist_main.csv', row.names = F)
write.csv(bhatt.dist.rand, 'tables/bhatt_dist_main_rand.csv', row.names = F)

bhatt.dist.relative <- read.csv('tables/bhatt_dist_main_relative.csv', check.names = F)
bhatt.dist <- read.csv('tables/bhatt_dist_main.csv', check.names = F)
bhatt.dist.rand <- read.csv('tables/bhatt_dist_main_rand.csv', check.names = F)
bhatt.dist.relative |> 
  pivot_longer(values_to = 'relative_distance', names_to = 'cell_type', cols = c(1:ncol(bhatt.dist.relative))) |> 
  group_by(cell_type) |> 
  mutate(median_distance = median(relative_distance, na.rm = TRUE)) |> 
  ungroup() |> 
  mutate(cell_type = fct_reorder(cell_type, median_distance)) |> 
  ggplot(aes(cell_type, relative_distance)) +
  geom_boxplot()

library(effsize) # for Cliff's Delta
library(transport) # for Wasserstein distance if needed

# List of cell types
cell_types <- c("CD4+T", "CD8+T", "NK", "Cycling T/NK", "B", "Plasma", 
                "Mast", "pDC", "cDC", "Mono", "Macro", 
                "Endo", "Mural", "CAF")

# Initialize result storage
results <- data.frame(CellType = character(),
                      Wilcoxon_p = numeric(),
                      KS_p = numeric(),
                      Permutation_p = numeric(),
                      CliffDelta = numeric(),
                      MedianDifference = numeric(),
                      HodgesLehmann = numeric())

# Permutation Test Function
perm_test <- function(dist1, dist2, n_perm = 1000) {
  observed_diff <- mean(dist1) - mean(dist2)
  combined <- c(dist1, dist2)
  labels <- rep(c("obs", "rand"), c(length(dist1), length(dist2)))
  
  perm_diffs <- replicate(n_perm, {
    perm_labels <- sample(labels)
    mean(combined[perm_labels == "obs"]) - mean(combined[perm_labels == "rand"])
  })
  
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
  return(p_value)
}

# Loop through each cell type
for (CT in cell_types) {
  # Extract Bhattacharyya distances for current cell type
  bhatt_obs <- bhatt.dist[, CT]
  bhatt_rand <- bhatt.dist.rand[, CT]
  
  # Skip if data is missing
  if (any(is.na(bhatt_obs)) | any(is.na(bhatt_rand))) next
  
  # Wilcoxon Test
  wilcox_p <- wilcox.test(bhatt_obs, bhatt_rand, alternative = "two.sided")$p.value
  
  # Kolmogorov-Smirnov Test
  ks_p <- ks.test(bhatt_obs, bhatt_rand)$p.value
  
  # Permutation Test
  perm_p <- perm_test(bhatt_obs, bhatt_rand)
  
  # Effect Sizes
  cliff_delta <- as.numeric(cliff.delta(bhatt_obs, bhatt_rand)$estimate)
  median_diff <- median(bhatt_obs) - median(bhatt_rand)
  hl_shift <- as.numeric(wilcox.test(bhatt_obs, bhatt_rand, conf.int = TRUE)$estimate)
  
  # Append results
  results <- rbind(results, data.frame(
    CellType = CT, 
    Wilcoxon_p = wilcox_p, 
    KS_p = ks_p, 
    Permutation_p = perm_p,
    CliffDelta = cliff_delta,
    MedianDifference = median_diff,
    HodgesLehmann = hl_shift
  ))
}

# Multiple Testing Corrections
results$Wilcoxon_adj <- p.adjust(results$Wilcoxon_p, method = "BH")  
results$KS_adj <- p.adjust(results$KS_p, method = "BH")
results$Permutation_adj <- p.adjust(results$Permutation_p, method = "BH")

# View Final Results
print(results)




