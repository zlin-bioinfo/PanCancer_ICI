pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','COSG','BPCells','SeuratExtend','MetBrewer','ggplot2','CytoTRACE2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(max.print = 10000)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)  # Set to 2 GiB
setwd("/home/zlin/workspace/PanCancer_ICI")

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
# rerun UMAP
seu <- RunUMAP(seu, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")
non_immune <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein","Pericytes","SMC", "Myofibroblasts","CAF-desmo",  "CAF_SFRP2", 
                "CAF-prog", "iCAF_MMP1", "iCAF_IL6", "CAF-ap","Cycling")
seu$celltype_r2 <- factor(seu$celltype_r2, levels = non_immune)
DimPlot(seu, group.by = 'celltype_r2', reduction = 'umap.full', alpha = 1, 
        cols = rev(color_pro(length(unique(seu$celltype_r2)), 1)), label = F) + 
  ggtitle('') + theme_void() + theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol=2, override.aes = list(size = 3)))
ggsave('figures/UMAP/UMAP_nonimmune.png', height = 4, width = 6, dpi = 300)

seu <- CreateSeuratObject(counts = GetAssayData(seu, layer = 'counts'), meta.data = seu@meta.data) |> NormalizeData()
genes_to_check <- c('PROX1', 'LYVE1','FLT4','CCL21','COL9A3',
                    'GJA5','FBLN5','GJA4',
                    'CA4','CD36','RGCC',
                    'COL4A1','KDR','ESM1','CXCR4',
                    'ACKR1','SELP','CLU',
                    'RGS5','ACTA2','MYH11',
                    'HOPX','TGFB1','HMGA2','HMGA1','CDKN2A','CDH2',
                    'FAP','POSTN','LRRC15','WNT5A','GREM1','SFRP2','SFRP4','IGF1',
                    'PI16','CD34','CD55','MFAP5',
                    'MMP1','MMP3','CXCL13','ISG15','IL7R','IL6','CXCL1','CXCL2','CXCL12','CEBPD','NFKB1',
                    'C7','ADAMDEC1','CD74','B2M','TOP2A','MKI67') 
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
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=50, expressed_pct=0.1)
write.csv(marker_cosg$name, 'tables/marker_nonimmune.csv', row.names = F)
# enrichment analysis
seu <- subset(seu, subset = celltype_r2 == 'Cycling', invert=T)
seu <- subset(seu, subset = celltype_main == 'CAF')
DefaultAssay(seu) <- 'RNA'
seu <- JoinLayers(seu)
seu <- CreateSeuratObject(counts = GetAssayData(seu, layer = 'counts'), meta.data = seu@meta.data) |> NormalizeData()
Idents(seu) <- seu$celltype_r2
marker_cosg <- cosg(seu, groups='all', assay='RNA', slot='data', mu=100, n_genes_user=100, expressed_pct=0.1)
marker_cosg$names |> head()
marker_list <- as.list(marker_cosg$names)
msigdbr(species = "Homo sapiens", category = 'H') 
runGSEA <- function(genes, universe=NULL,
                    category="H", subcategory=NULL,
                    species="Homo sapiens",
                    custom.db=NULL,
                    pval.thr=0.05) {
  
  
  if (!requireNamespace("fgsea", quietly = TRUE) |
      !requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Function 'runGSEA' requires the 'fgsea' and 'msigdbr' packages.
            Please install them.", call. = FALSE)
  }  
  
  if (any(duplicated(genes))) {
    genes <- genes[!duplicated(genes)]
  }
  
  if (!is.null(custom.db)) {
    #check format: should be a named list of signatures
    if (!is.list(custom.db) || is.null(names(custom.db))) 
      stop("custom.db should be a named list", call. = FALSE)
    DB_list <- custom.db
  } else {  # use signatures from mSigDB
    msig_df <- msigdbr::msigdbr(species = species, category = category, subcategory=subcategory)
    DB_list <- split(x=msig_df$gene_symbol, f=msig_df$gs_name)
  }
  
  fgRes <- fgsea::fora(pathways = DB_list,
                       genes = genes,
                       universe = universe)
  
  fgRes <- fgRes[fgRes$pval <= pval.thr,]
  return(fgRes)
}

marker_gene_list <- as.list(marker_cosg$names)

res_list <- list()
for (i in 1:ncol(marker_cosg$names)){
  print(names(marker_cosg$names)[i])
  hallmark <- runGSEA(marker_cosg$names[,i], universe = rownames(seu), category = 'H') |> data.frame() |> filter(padj<0.05) 
  # reactome <- runGSEA(marker_cosg$names[,i], universe = rownames(seu), category = 'C2', subcategory = "CP:REACTOME") |> data.frame() |> filter(padj<0.05)  
  # wiki <- runGSEA(marker_cosg$names[,i], universe = rownames(seu), category = 'C2', subcategory = "CP:WIKIPATHWAYS") |> data.frame() |> filter(padj<0.05)
  kegg <- runGSEA(marker_cosg$names[,i], universe = rownames(seu), category = 'C2', subcategory = "CP:KEGG") |> data.frame() 
  c5 <- runGSEA(marker_cosg$names[,i], universe = rownames(seu), category = 'C5') |> data.frame() |> filter(padj<0.05) 
  res <- rbind(hallmark,kegg,c5)
  res$celltype <- names(marker_cosg$names)[i]
  res_list[[i]] <- res
}

res <- do.call(rbind, res_list) |> filter(overlap > 15)
qs_save(res, 'tables/marker_caf.qs2')

caf_pathways <- list(
  # Core shared pathways with annotations
  Shared = c(
    "GOBP_TISSUE_DEVELOPMENT",                         # Pro-tumor: stromal support
    "GOBP_VASCULATURE_DEVELOPMENT",                    # Pro-tumor: angiogenesis
    "GOBP_CELL_MIGRATION",                             # Pro-tumor: invasion
    "GOBP_REGULATION_OF_CELL_POPULATION_PROLIFERATION", # Pro-tumor: growth
    "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX" # Pro-tumor: ECM; Anti-immune: barrier
  ),
  
  # Distinct pathways with annotations
  iCAF_MMP1 = c(
    "GOBP_CYTOKINE_PRODUCTION",               # Pro-tumor, Pro-immune, Anti-immune
    "GOBP_INNATE_IMMUNE_RESPONSE"            # Pro-tumor, Pro-immune, Anti-immune
  ),
  Myofibroblasts = c(
    "GOBP_BLOOD_VESSEL_MORPHOGENESIS",                # Pro-tumor: vascularization
    "GOBP_REGULATION_OF_CELLULAR_COMPONENT_MOVEMENT"   # Pro-tumor: migration regulation
  ),
  iCAF_IL6 = c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",       # Pro-tumor, Pro-immune, Anti-immune
    "GOBP_RESPONSE_TO_LIPID"                          # Pro-tumor: lipid signaling
  ),
  CAF_desmo = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",      # Pro-tumor: EMT
    "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT" # Pro-tumor, Anti-immune
  ),
  CAF_SFRP2 = c(
    "GOBP_ENZYME_LINKED_RECEPTOR_PROTEIN_SIGNALING_PATHWAY", # Pro-tumor: kinase signaling
    "GOBP_ANIMAL_ORGAN_MORPHOGENESIS"                 # Pro-tumor: structural support
  ),
  CAF_ap = c(
    "GOBP_DEFENSE_RESPONSE",                  # Pro-tumor, Pro-immune, Anti-immune
    "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS" # Pro-tumor, Pro-immune, Anti-immune
  ),
  CAF_prog = c(
    "GOBP_POSITIVE_REGULATION_OF_DEVELOPMENTAL_PROCESS",   # Pro-tumor: developmental support
    "GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_STIMULUS" # Pro-tumor, Anti-immune
  )
)
pathways <- unlist(tme_caf_pathways) |> unique()





