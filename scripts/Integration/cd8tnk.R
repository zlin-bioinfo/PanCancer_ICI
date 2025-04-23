pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','COSG','BPCells','SeuratExtend','MetBrewer','ggplot2','CytoTRACE2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(max.print = 10000)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)  # Set to 2 GiB
setwd("/home/zlin/workspace/PanCancer_ICI")

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
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters_full', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  RotatedAxis()
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)), label = T)
seu$celltype_r2[seu$seurat_clusters_full == 5] <- 'Cycling'
seu$celltype_r2[seu$seurat_clusters_full %in% c(11,12,14)] <- 'Doublet'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$celltype_bped_main), c('CD4+ T-cells', 'CD8+ T-cells', 'NK cells','unknown', NA))] <- 'Doublet'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$scGate_multi), c('CD8T', 'CD4T', 'NK','unknown', NA))] <- 'Doublet'
seu$cell.id[seu$cohort == 'NSCLC_Liu'] <- str_replace(colnames(seu)[seu$cohort == 'NSCLC_Liu'],'_15','')
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
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
ggsave('figures/UMAP/UMAP_CD8T.png', height = 4, width = 6.5, device = 'png')

genes_to_check <- c('CD3D','CD8A','CD8B','CD4', 'TCF7','LEF1','CCR7','SELL',
                    'IL7R', 'ANXA1', 'CD55', 'CD27','LMNA','ZNF683','CXCR6','ITGAE','CXCR5',
                    'GZMB','GZMH','GZMK','EOMES','CXCR3','IFNG','PRF1','CXCL13',
                    'CTLA4','PDCD1','LAG3','TIGIT','LAYN','ISG15','IFIT1','IFIT3', "HSPA1A","HSPA6","NR4A1",# 'MYL12A','MYL12B',
                    'SLC4A10','KLRB1','RORA',
                    'TRDV2','TRGV9','CX3CR1','TBX21','FGFBP2','KLRF1','FCGR3A',
                    'GNLY','TYROBP','XCL1','XCL2','NCAM1','MKI67','TOP2A','RPL11','RPL10A','CD19','CD79A')
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
DimPlot(seu, group.by = 'celltype_main', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$cohort)), 1)),
        label = T)

seu <- CreateSeuratObject(counts = GetAssayData(seu, layer = 'counts'), meta.data = seu@meta.data) |> NormalizeData()
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
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_cd8tnk.csv', row.names = F)



