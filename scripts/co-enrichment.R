rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','dittoSeq','RColorBrewer','tibble','pheatmap','qs','MetBrewer','viridis','forcats','grid','corrplot','ComplexHeatmap','colorRamp2','corrr','igraph','ggraph','tidygraph','graphlayouts','ggforce','ggnetwork','ggnewscale')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

meta_int <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_int.csv')
pt_df <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_patient.csv')
head(meta_int)
range(meta_int$freq_r2_comp)
# Immune cells
immune <- c('CD4_Naive','CD4_Tm_CREM-','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM', 'CD4_Tm_CCL5', 
            'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
            'CD4_Treg_Early', 'CD4_Treg_ISG15', 'CD4_Treg_TNFRSF9', 
            'CD4_Th_ISG15', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif',
            'CD8_Prolif', 'CD8_Naive', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
            'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13',
            'CD8_Tex_ISG15', 'CD8_Temra_CX3CR1', 'CD8_NK-like', 
            'MAIT', 'gdT', 'NK_CD56loCD16hi', 'NK_CD56hiCD16lo',
            'B_Naive', 'B_ISG15', 'B_HSP', 'B_MT2A', 'ACB_EGR1', 'ACB_NR4A2', 'ACB_CCR7', 'B_Memory', 'B_AtM',
            'GCB_Pre', 'GCB_SUGCT', 'GCB_LMO2', 'GCB_Prolif', 'Plasmablast', 'Plasma_cell',
            'Mast','pDC','cDC1', 
            'cDC2_CD1C', 'cDC2_IL1B','cDC2_ISG15', 'cDC2_CXCL9', 'DC_LC-like', 'MigrDC', 'MoDC', 
            'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
            'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro_ISG15', 
            'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2')

# component
pdf('/bigdata/zlin/Melanoma_meta/figures/Co-enrichment/immune.pdf', height = 10, width = 10)
mat <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  pivot_longer(cols = c('Pre', 'Post'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  filter(dataset != 'NSCLC_Liu', 
         celltype_r2 %in% immune) |> 
  mutate(sample = paste0(patient, '_', time_point)) |> 
  select(sample, freq_r2_comp, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = sample, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
M <- cor(t(mat))
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
         tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
         mar = c(0,0,0.7,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = paste0('Immune Cells (n=', ncol(mat), ')'), diag = F, method = 'square')
dev.off()


pdf('/bigdata/zlin/Melanoma_meta/figures/Co-enrichment/immune_pre_post.pdf', height = 10, width = 18)
par(mfrow=c(1,2))
mat <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  filter(dataset != 'NSCLC_Liu', 
         celltype_r2 %in% immune) |>
  select(patient, Pre, celltype_r2) |> 
  pivot_wider(values_from = Pre, names_from = patient, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
M <- cor(t(mat))
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
         tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
         mar = c(0,0,0.7,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = 'Immune Cells(Pre-Tx)', diag = F, method = 'square')

mat <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  filter(dataset != 'NSCLC_Liu', 
         celltype_r2 %in% immune) |>
  select(patient, Post, celltype_r2) |> 
  pivot_wider(values_from = Post, names_from = patient, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
M <- cor(t(mat))
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
         tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
         mar = c(0,0,0.7,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = 'Immune Cells(Post-Tx)', diag = F, method = 'square')
dev.off()

pdf('/bigdata/zlin/Melanoma_meta/figures/Co-enrichment/TME.pdf', height = 10, width = 10)
mat <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  pivot_longer(cols = c('Pre', 'Post'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley', 'HNSC_Franken'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014")) |> 
  mutate(sample = paste0(patient, '_', time_point)) |> 
  select(sample, freq_r2_comp, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = sample, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
M <- cor(t(mat))
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
         tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
         mar = c(0,0,0.7,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = paste0('TME (n=', ncol(mat), ')'), diag = F, method = 'square')
dev.off()


pdf('/bigdata/zlin/Melanoma_meta/figures/Co-enrichment/TME_pre_post.pdf', height = 10, width = 18)
par(mfrow=c(1,2))
mat <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley', 'HNSC_Franken'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014"), 
         celltype_r2 != 'B_MT2A') |> 
  select(patient, Pre, celltype_r2) |> 
  pivot_wider(values_from = Pre, names_from = patient, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
M <- cor(t(mat))
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
         tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
         mar = c(0,0,0.7,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = 'TME (Pre-Tx)', diag = F, method = 'square')

mat <- meta_int |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley', 'HNSC_Franken'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014"), 
         celltype_r2 != 'B_MT2A') |> 
  select(patient, Post, celltype_r2) |> 
  pivot_wider(values_from = Post, names_from = patient, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
M <- cor(t(mat))
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cols <- colorRampPalette(c("#336699", "white", "#CC0000")) 
corrplot(M, p.mat = testRes$p, col = cols(100), tl.srt=45,
         tl.cex = 0.5, pch.cex = 0.6, tl.col = 'black', order = 'hclust', hclust.method = "ward.D2",
         mar = c(0,0,0.7,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = 'TME (Post-Tx)', diag = F, method = 'square')
dev.off()



