rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','ggsci','ggplot2','tibble','qs2','Matrix','janitor','RColorBrewer','SeuratExtend','COSG')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
setwd("/home/zlin/workspace/PanCancer_ICI")
source('scripts/Celltype_classification.R')
metadata <- read.csv('tables/meta_all.csv')
filter_sample <- metadata |> 
  distinct(sample, .keep_all = T) |> 
  pull(sample)
# CD4+T
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 'SCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'NSCLC_Liu',
              'PCa_Hawley')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = celltype_r2 %in% c('CD4_T-naive','CD4_Tcm','CD4_Treg','CD4_T-ISG','CD4_Tfh','CD4_Tstr','CD4_Tctl','CD4_Th17')) |> 
    subset(subset = sample %in% filter_sample)
  seu$celltype_r2[seu$celltype_main %in% c('Cycling', 'Cycling T/NK')] <- 'CD4_T-cycling'
  seu <- NormalizeData(seu)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_cd4t.csv', row.names = F)
genes_to_check <- c('TCF7','LEF1','CCR7','SELL',
                    'IL7R','CXCR4','CD55','GPR183','CD69','LMNA','ANXA1',
                    'DUSP2','NR4A2',
                    'TNF','GZMA','GZMB','GZMK','TBX21','CCL4','CCL5','IFNG','PRF1','CX3CR1','KLRG1','EOMES', #IL2	CD40LG	TBX21	STAT4	IL12RB2
                    'CD40LG','STAT4','CXCR5','IL21','BCL6','CXCL13',"GNG4", "CD200", "IGFL2",
                    'PDCD1','LAG3','HAVCR2','IL2RA','CTLA4','LAYN','RTKN2', 'FOXP3','TNFRSF9',
                    'CXCR6','CCR6','KLRB1','RORA','RORC','IL17A','IL26',
                    'ISG15',"IFIT1","IFIT3",'IFI44L',
                    "HSPA1B","HSPA6", 'FOS','JUN','MKI67','TOP2A')
seu$celltype_r2 <- factor(seu$celltype_r2, levels = c('CD4_T-naive','CD4_Tcm','CD4_Tctl','CD4_Tfh','CD4_Treg','CD4_Th17','CD4_T-ISG','CD4_Tstr','CD4_T-cycling'))
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
# qs_save(seu, 'data/seu_cd4.qs2')

# CD8+T
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 'SCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'NSCLC_Liu',
              'PCa_Hawley')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = celltype_r2 %in% c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
                                       "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
                                       "CD8_NK-like",'MAIT','gdT','NK_CD56loCD16hi','NK_CD56hiCD16lo')) |> 
    subset(subset = sample %in% filter_sample)
  seu$celltype_r2[seu$celltype_r2 %in% c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
                                     "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
                                     "CD8_NK-like",'MAIT','gdT') & seu$celltype_main %in% c('Cycling', 'Cycling T/NK')] <- 'CD8_T-cycling'
  seu$celltype_r2[seu$celltype_r2 %in% c('NK_CD56loCD16hi','NK_CD56hiCD16lo') & seu$celltype_main %in% c('Cycling', 'Cycling T/NK')] <- 'NK-cycling'
  seu <- NormalizeData(seu)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_cd8t.csv', row.names = F)
genes_to_check <- c('CD3D','CD8A','CD8B',
                    'TCF7','LEF1','CCR7','SELL',
                    'IL7R', 'ANXA1', 'CD55', 'CD27','LMNA','ZNF683','CXCR6','ITGAE',
                    'GZMB','GZMH','GZMK','EOMES','CXCR3','IFNG','PRF1','CXCL13',
                    'CTLA4','PDCD1','LAG3','TIGIT','LAYN','ISG15','IFIT1','IFIT3', "HSPA1A","HSPA6","NR4A1",# 'MYL12A','MYL12B',
                    'SLC4A10','KLRB1','RORA',
                    'TRDV2','TRGV9','CX3CR1','TBX21','FGFBP2','KLRF1','FCGR3A',
                    'GNLY','TYROBP','XCL1','XCL2','NCAM1','MKI67','TOP2A')
seu$celltype_r2 <- factor(seu$celltype_r2, levels = c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
                                                      "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
                                                      "CD8_NK-like",'MAIT','CD8_T-cycling','gdT','NK_CD56loCD16hi','NK_CD56hiCD16lo', 'NK-cycling'))
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
ggsave('figures/Dotplot/CD8TNK.pdf', height = 12, width = 6.5)

# B/plasma cells
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 
              'PCa_Hawley')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = celltype_r2 %in% bplasma) |> 
    subset(subset = sample %in% filter_sample)
  seu <- NormalizeData(seu)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_bplasma.csv', row.names = F)
seu$celltype_r2 <- factor(seu$celltype_r2, levels = bplasma)
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

# Myeloid cells
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 
              'PCa_Hawley')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = celltype_r2 %in% mye) |> 
    subset(subset = sample %in% filter_sample)
  seu$celltype_r2[seu$celltype_main %in% 'Cycling'] <- "Cycling"
  seu <- NormalizeData(seu)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_myeloids.csv', row.names = F)
seu$celltype_r2 <- factor(seu$celltype_r2, levels = c('Mast','pDC','cDC1', 
                                                      'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'mregDC', 'MoDC', 
                                                      'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
                                                      'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro-ISG', 
                                                      'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2','Cycling'))
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
                    'C1QC',"MARCKS",'APOE','TREM2','GPNMB','MMP9','CD36','HLA-DQA1','TOP2A','MKI67')
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
ggsave('figures/Dotplot/Myeloids.pdf', height = 12, width = 6.5)

# Non-immune
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao',
              'HNSC_Franken',
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 
              'PCa_Hawley')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = celltype_r2 %in% nonimmune) |> 
    subset(subset = sample %in% filter_sample)
  seu$celltype_r2[seu$celltype_main %in% 'Cycling'] <- "Cycling"
  seu <- NormalizeData(seu)
  return(seu)
})
seu <- merge(x = seu_list[[1]], y=seu_list[2:length(seu_list)])
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_nonimmune.csv', row.names = F)
genes_to_check <- c('PROX1', 'LYVE1','FLT4','CCL21',
                   'GJA5','FBLN5','GJA4',
                   'CA4','CD36','RGCC',
                   'COL4A1','KDR','ESM1','CXCR4',
                   'ACKR1','SELP','CLU',
                   'RGS5','ACTA2','MYH11',
                   'HOPX','TGFB1',
                   'SFRP2','SFRP4','IGF1',
                   'PI16','CD34','CD55','MFAP5',
                   'POSTN','FAP','LRRC15',
                   'WNT5A','GREM1',
                   'MMP1','MMP3','CXCL13','ISG15','IL6','CXCL1','CXCL2','CEBPD','NFKB1',
                   'ADAMDEC1','ADAM28','C7','CD74','B2M','MKI67','TOP2A') 
seu$celltype_r2 <- factor(seu$celltype_r2, levels = c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein",
                                                      "Pericytes","SMC",  "CAF_SFRP2", 
                                                      "CAF-prog","Myofibroblasts", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap", 'Cycling'))
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
ggsave('figures/Dotplot/Non-immune.pdf', height = 10, width = 5.5)


