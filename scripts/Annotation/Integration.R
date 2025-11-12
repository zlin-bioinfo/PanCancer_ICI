pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','COSG','BPCells','SeuratExtend','MetBrewer','ggplot2','CytoTRACE2','AnnotationDbi','org.Hs.eg.db')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(max.print = 10000)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)  # Set to 2 GiB

# CD4+T cells
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
seu$celltype_r2 <- factor(seu$celltype_r2, levels =  c('CD4_T-naive','CD4_Tcm','CD4_Tstr','CD4_Tctl','CD4_Tfh','CD4_Th17','CD4_Treg','CD4_T-ISG','Cycling'))
p <- DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', 
             cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=1, override.aes = list(size = 3)))
p + theme(legend.position = "none") 
ggsave('figures/UMAP/UMAP_CD4T.png', height = 4, width = 4, dpi = 300)
get_legend(p) |> as_ggplot()
ggsave('figures/UMAP/legend_CD4T.pdf')

seu <- CreateSeuratObject(counts = GetAssayData(seu, layer = 'counts'), meta.data = seu@meta.data) |> NormalizeData()
seu$celltype_r2 <- factor(seu$celltype_r2, levels =  c('CD4_T-naive','CD4_Tcm','CD4_Tfh','CD4_Th17','CD4_Treg','CD4_T-ISG','CD4_Tctl','CD4_Tstr','Cycling'))
genes_to_check <- c('TCF7','LEF1','CCR7','SELL',
                    'IL7R','CXCR4','CD55','GPR183','CD69','LMNA','ANXA1','DUSP2','NR4A2',
                    'CD40LG','STAT4','CXCR5','IL21','BCL6','CXCL13','NMB',"GNG4", "CD200", "IGFL2",'PDCD1','LAG3',
                    'BATF','KLRB1','IL17A','IL26','RORC','RORA','CXCR6','CCR6',
                    'FOXP3','TNFRSF9','HAVCR2','IL2RA','CTLA4','LAYN','RTKN2', 
                    'ISG15',"IFIT1","IFIT3",'IFI44L',
                    'TNF','GZMA','GZMB','GZMK','TBX21','CCL4','CCL5','IFNG','PRF1','CX3CR1','KLRG1','EOMES', #IL2	CD40LG	TBX21	STAT4	IL12RB2
                    "HSPA1B","HSPA6", 'FOS','JUN','MKI67','TOP2A')
DotPlot(seu,genes_to_check, group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() + 
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
ggsave('figures/Dotplot/CD4T.pdf', height = 3.5, width = 12)

gene_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = rownames(seu),
                                          columns = c("ENSEMBL", "GENETYPE", "GENENAME"), # Include other useful columns
                                          keytype = "SYMBOL")
protein_coding_rna <- gene_annotations |> filter(GENETYPE %in% c('protein-coding','other')) |> pull(SYMBOL)
seu <- seu[protein_coding_rna,]
seu <- seu[!str_detect(rownames(seu), '^MT-|^RP'),]
markers_wilcox <- FindAllMarkers(seu, group.by = 'celltype_r2', only.pos = T) |> mutate(pct.diff = pct.1-pct.2)
marker <- markers_wilcox |> group_by(cluster) |> slice_head(n=20) |> ungroup()
write.csv(markers_wilcox, 'tables/markers_cd4t.csv', row.names = F)

# CD8+T/NK cells
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
adjusted_cd8tnk <- seu@meta.data |> 
  filter(celltype_r2 %in% c('Cycling','Doublet')) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(adjusted_cd8tnk, 'tables/adjusted_cd8tnk.csv')
seu <- subset(seu, subset = celltype_r2 == 'Doublet', invert = T)

genes_to_check <- c('CD3D','CD8A','CD8B','CD4', 'TCF7','LEF1','CCR7','SELL',
                    'IL7R', 'ANXA1', 'CD55', 'CD27','LMNA','ZNF683','CXCR6','ITGAE','CXCR5',
                    'GZMB','GZMH','GZMK','EOMES','CXCR3','IFNG','PRF1','CXCL13',
                    'CTLA4','PDCD1','LAG3','TIGIT','LAYN','ISG15','IFIT1','IFIT3', "HSPA1A","HSPA6","NR4A1",# 'MYL12A','MYL12B',
                    'SLC4A10','KLRB1','RORA',
                    'TRDV2','TRGV9','CX3CR1','TBX21','FGFBP2','KLRF1','FCGR3A',
                    'GNLY','TYROBP','XCL1','XCL2','NCAM1','MKI67','TOP2A','RPL11','RPL10A','CD19','CD79A')
DotPlot(seu, genes_to_check, group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
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
         size = guide_legend(direction = 'vertical',
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
seu$celltype_r2 <- as.character(seu$celltype_r2)
seu$celltype_r2[seu$celltype_r2 == 'MAIT'] = 'CD8_Tc17'
seu$celltype_r2 <- factor(seu$celltype_r2, levels = c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
                                                      "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
                                                      "CD8_NK-like","CD8_Tc17",'gdT','NK_CD56loCD16hi','NK_CD56hiCD16lo','Cycling'))
genes_to_check <- c('CD3D','CD8A','CD8B',
                    'TCF7','LEF1','CCR7','SELL',
                    'IL7R', 'ANXA1', 'CD55', 'CD27','LMNA','ZNF683','CXCR6','ITGAE','CXCR5',
                    'GZMB','GZMH','GZMK','EOMES','CXCR3','IFNG','PRF1','CXCL13',
                    'CTLA4','PDCD1','LAG3','TIGIT','LAYN',
                    'GNG4','CD200',
                    'ISG15','IFIT1','IFIT3', "HSPA1A","HSPA6","NR4A1",# 'MYL12A','MYL12B',
                    'SLC4A10','KLRB1','RORA','RORC',"TRAV1-2",'TRDV1','TRDV2','TRGV9','CX3CR1','TBX21','FGFBP2','KLRF1','FCGR3A',
                    'GNLY','TYROBP','XCL1','XCL2','NCAM1','MKI67','TOP2A')
DotPlot(seu, genes_to_check, group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() + 
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
ggsave('figures/Dotplot/CD8TNK.pdf', height = 4, width = 13)

gene_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = rownames(seu),
                                          columns = c("ENSEMBL", "GENETYPE", "GENENAME"), # Include other useful columns
                                          keytype = "SYMBOL")
protein_coding_rna <- gene_annotations |> filter(GENETYPE %in% c('protein-coding','other')) |> pull(SYMBOL)
seu <- seu[protein_coding_rna,]
seu <- seu[!str_detect(rownames(seu), '^MT-|^RP'),]
markers_wilcox <- FindAllMarkers(seu, group.by = 'celltype_r2', only.pos = T) |> mutate(pct.diff = pct.1-pct.2)
write.csv(markers_wilcox, 'tables/markers_cd8tnk.csv', row.names = F)

# Myeloid cells
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
                    'C1QC',"MARCKS",'APOE','TREM2','GPNMB','MMP9','TOP2A','MKI67') # 'TGFB1','CD80','CD86','RPL10A','RPL11'
myeloids <- c('Mast','pDC','cDC1', 'cDC2', 'mregDC',
              'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
              'Macro_IL1B', 'Macro_INHBA', 'Macro_FN1', 'Macro_SPP1', 'Macro-ISG', 
              'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2', 'Cycling')
seu$celltype_r2 <- factor(seu$celltype_r2, levels = myeloids)
DotPlot(seu, genes_to_check, group.by = 'seurat_clusters_full', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() + 
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
ggsave('figures/Dotplot/Myeloids.pdf', height = 4.5, width = 15)

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
                    'CD1C','FCER1A','IL1B','NLRP3','NFKBIA','ICAM1','STAT3','IL6',
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

gene_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = rownames(seu),
                                          columns = c("ENSEMBL", "GENETYPE", "GENENAME"), # Include other useful columns
                                          keytype = "SYMBOL")
protein_coding_rna <- gene_annotations |> filter(GENETYPE %in% c('protein-coding')) |> pull(SYMBOL)
seu <- seu[protein_coding_rna,]
seu <- seu[!str_detect(rownames(seu), '^MT-|^RP'),]
Idents(seu) <- seu$celltype_r2
markers_wilcox <- FindAllMarkers(seu, group.by = 'celltype_r2', only.pos = T) |> mutate(pct.diff = pct.1-pct.2)
marker <- markers_wilcox |> group_by(cluster) |> slice_head(n=30) |> ungroup()
write.csv(markers_wilcox, 'tables/markers_myeloids.csv', row.names = F)

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


gene_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = rownames(seu),
                                          columns = c("ENSEMBL", "GENETYPE", "GENENAME"), # Include other useful columns
                                          keytype = "SYMBOL")
protein_coding_rna <- gene_annotations |> filter(GENETYPE %in% c('protein-coding','other')) |> pull(SYMBOL)
seu <- seu[protein_coding_rna,]
seu <- seu[!str_detect(rownames(seu), '^MT-|^RP'),]
markers_wilcox <- FindAllMarkers(seu, group.by = 'celltype_r2', only.pos = T) |> mutate(pct.diff = pct.1-pct.2)
marker <- markers_wilcox |> group_by(cluster) |> slice_head(n=30) |> ungroup()
write.csv(markers_wilcox, 'tables/markers_bplasma.csv', row.names = F)

# Non-immune cells
non_immune <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein","Pericytes","SMC", "Myofibroblasts", "CAF_SFRP2", 
                "CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap")
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao',
              'HNSC_Franken',  
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 
              'PCa_Hawley','RCC_Bi','HCC_Ma')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) |> 
    subset(subset = celltype_r2 %in% non_immune) 
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu <- seu |> 
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
DotPlot(seu, unlist(genes_to_check), group.by = 'seurat_clusters_full', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') +
  RotatedAxis()
DimPlot(seu, group.by = 'seurat_clusters_full', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$seurat_clusters_full)), 1)), label = T)
seu$celltype_r2[seu$seurat_clusters_full == 8] <- 'Cycling'
seu$celltype_r2[seu$seurat_clusters_full %in% c(7,9,11,13)] <- 'Doublet'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$celltype_bped_main), c("Fibroblasts", "Endothelial cells", "Myocytes",'unknown', NA))] <- 'Doublet'
seu$celltype_r2[seu$scGate_multi %in% setdiff(unique(seu$scGate_multi), c("Fibroblast", "Endothelial",'unknown', NA))] <- 'Doublet'
adjusted_nonimmune <- seu@meta.data |> 
  filter(celltype_r2 %in% c('Cycling','Doublet')) |> 
  select(cell.id, celltype_r2, seurat_clusters_full)
write.csv(adjusted_nonimmune, 'tables/adjusted_nonimmune.csv')
seu <- subset(seu, subset = celltype_r2 == 'Doublet', invert = T)

seu <- CreateSeuratObject(counts = GetAssayData(seu, layer = 'counts'), meta.data = seu@meta.data) |> NormalizeData()
genes_to_check <- c('PROX1', 'LYVE1','FLT4','CCL21','COL9A3',
                    'GJA5','FBLN5','GJA4',
                    'CA4','CD36','RGCC',
                    'COL4A1','KDR','ESM1','CXCR4',
                    'ACKR1','SELP','CLU',
                    'RGS5','ACTA2','MYH11',
                    "COL1A1", 'LUM','DCN', 'TIMP1',
                    'TGFB1','CDH2',
                    'FAP','POSTN','LRRC15','WNT5A','GREM1','SFRP2','SFRP4','IGF1',
                    'PI16','CD34','CD55','MFAP5',
                    'MMP1','MMP3','CXCL13','ISG15','IL7R','IL6','CXCL1','CXCL2','CXCL12','CEBPD','NFKB1',
                    'C7','ADAMDEC1','CD74','B2M',"HLA-A",'HLA-DRB1', 'HLA-DRA','TOP2A','MKI67','RAMP1') 
DotPlot(seu, genes_to_check, group.by = 'celltype_r2', col.min = -1, col.max = 1, dot.scale = 5, cols = 'RdBu') + 
  theme_minimal() + 
  # scale_y_discrete(position = "right") +
  RotatedAxis() +
  # coord_flip() + 
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
ggsave('figures/Dotplot/Non-immune.pdf', height = 3.5, width = 12)
Idents(seu) <- seu$celltype_r2

gene_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = rownames(seu),
                                          columns = c("ENSEMBL", "GENETYPE", "GENENAME"), # Include other useful columns
                                          keytype = "SYMBOL")
protein_coding_rna <- gene_annotations |> filter(GENETYPE == 'protein-coding') |> pull(SYMBOL)
seu <- seu[protein_coding_rna,]
seu <- seu[!str_detect(rownames(seu), '^MT-|^RP-'),]
seu_sub <- subset(seu, subset=celltype_r2 %in% c("Pericytes","SMC", "Myofibroblasts","CAF-desmo",  "CAF_SFRP2", 
                                                 "CAF-prog", "iCAF_MMP1", "iCAF_IL6", "CAF-ap"))
markers_wilcox <- FindAllMarkers(seu_sub, group.by = 'celltype_r2', only.pos = T)
write.csv(markers_wilcox, 'tables/markers_caf.csv', row.names = F)
seu_sub <- subset(seu, subset=celltype_r2 %in% c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein"))
markers_wilcox <- FindAllMarkers(seu_sub, group.by = 'celltype_r2', only.pos = T)
write.csv(markers_wilcox, 'tables/markers_endo.csv', row.names = F)




