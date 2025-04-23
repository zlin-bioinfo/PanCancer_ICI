pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','COSG','BPCells','SeuratExtend','MetBrewer','ggplot2','CytoTRACE2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(max.print = 10000)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)  # Set to 2 GiB
setwd("/home/zlin/workspace/PanCancer_ICI")

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
  RotatedAxis()
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 0.5, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)),
        label = T)

DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = T)

seu$celltype_r2[seu$seurat_clusters_full == 9] <- 'Cycling'
seu$celltype_r2[seu$seurat_clusters_full %in% c(7,14,15)] <- 'Doublet'
seu$celltype_r2[seu$seurat_clusters_full %in% c(6, 16)] <- 'pDC'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$celltype_bped_main), c("Monocytes","DC", "Macrophages", "Neutrophils", 'unknown', NA))] <- 'Doublet'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$scGate_multi), c("Macrophage",  "Monocyte", "panDC","Mast","Neutrophils",'unknown', NA))] <- 'Doublet'
adjusted_myeloids <- seu@meta.data |> 
  filter(celltype_r2 %in% c('Cycling','Doublet','pDC')) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(adjusted_myeloids, 'tables/adjusted_myeloids.csv')
seu <- subset(seu, subset = celltype_r2 == 'Doublet', invert = T)
# rerun UMAP
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
myeloids <- c('Mast','pDC','cDC1', 
              'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'mregDC', 'MoDC', 
              'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
              'Macro_IL1B', 'Macro_INHBA', 'Macro_FN1', 'Macro_SPP1', 'Macro-ISG', 
              'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2', 'Cycling')
seu$celltype_r2 <- factor(seu$celltype_r2, levels = myeloids)
p <- DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
             cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
qs_save(p, 'data/umap_myeloids.qs2')
p <- qs_read('data/umap_myeloids.qs2')
p
ggsave('figures/UMAP/UMAP_myeloids.png', height = 4, width = 6, dpi = 300)


seu <- CreateSeuratObject(counts = GetAssayData(seu, layer = 'counts'), meta.data = seu@meta.data) |> NormalizeData()
genes_to_check <- c('TPSAB1','LILRA4','SELL',
                    'CLEC9A','XCR1',
                    'CD1C','FCER1A','IL1B','NLRP3','NFKBIA','ICAM1',
                    'CD1A','CD207','LAMP3','CCR7','IDO1','NFKB1','DUSP4','MARCKSL1','BIRC3',
                    'S100A8','S100A9','FCN1','CD14','CD300E','FCGR3A','LST1','LILRB2',
                    'CXCL2','TREM1','AREG','EREG',
                    'VCAN','VEGFA','INHBA','FN1','SPP1','OLR1',
                    'ISG15','IFIT1','IFIT2','IFI6','CXCL9','CXCL10','TNFRSF14','IL4I1',
                    'CCL2','CCL3','CCL4','DUSP1','TNF','JUNB','ATF3','CCL3L1',
                    'NR4A2','LYVE1','FOLR2','MRC1','MSR1',
                    'C1QC',"MARCKS",'APOE','TREM2','GPNMB','MMP9','TOP2A','MKI67','RPL10A','RPL11') # 'TGFB1','CD80','CD86'
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
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_myeloids.csv', row.names = F)



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
seu <- JoinLayers(seu)
seu <- subset(seu, subset = celltype_main == 'Mono/macro')
DotPlot2(seu, features = c("CD274", "PDCD1LG2","FAS", "CD200","SOCS1","SOCS2"),group.by = 'celltype_r2',color_scheme = "BuRd")


