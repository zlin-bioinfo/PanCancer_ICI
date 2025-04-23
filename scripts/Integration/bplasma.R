pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','COSG','BPCells','SeuratExtend','MetBrewer','ggplot2','CytoTRACE2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(max.print = 10000)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)  # Set to 2 GiB
setwd("/home/zlin/workspace/PanCancer_ICI")

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
seu <- JoinLayers(seu)
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
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$celltype_bped_main), c('B-cells', 'unknown', NA))] <- 'Doublet'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$scGate_multi), c("Bcell", "PlasmaCell",'unknown', NA))] <- 'Doublet'
adjusted_bplasma <- seu@meta.data |> 
  filter(celltype_r2 %in% c("Mixed", "PC-trans", "Doublet")) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(adjusted_bplasma, 'tables/adjusted_bplasma.csv')
seu <- subset(seu, subset = celltype_r2 %in% c('Mixed', 'Doublet'), invert = T)
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


