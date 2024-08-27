rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','dittoSeq','RColorBrewer','tibble','pheatmap','qs','MetBrewer','viridis','ComplexHeatmap','colorRamp2','corrr','ggnewscale','NMF','ggpubr','corrplot')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

meta_int <- read.csv('/bigdata/zlin/PanCancer_ICI/tables/meta_int.csv') 
pt_df <- read.csv('/bigdata/zlin/PanCancer_ICI/tables/meta_patient.csv')

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
nonimmune <- c('Endo_lymphatic','Endo_artery','Endo_capillary','Endo_tip','Endo_vein',
               'EndMT','CAF_inflammatory', 'CAF_adipogenic', 'CAF_PN', 'CAF_AP', 'Myofibroblast')

mtx_immune <- meta_int |> 
  filter(!dataset %in% c('NSCLC_Liu'), 
         !patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su011', 'BCC/SCC_Yost_su012', 'BCC/SCC_Yost_su014'),
         celltype_r2 %in% immune) |> 
  select(sample,  freq_r2_comp, celltype_r2) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = celltype_r2, values_fill = 0) |> 
  column_to_rownames(var = 'sample') |> t()

res.nmf <- nmf(mtx_immune, rank = 2:8, nrun = 40, seed = 10000)
pdf('/bigdata/zlin/PanCancer_ICI/figures/NMF/rank_selection.pdf', height = 6, width =8)
plot(res.nmf)
dev.off()
nrank=3
nmf.res = nmf(mtx_immune, rank=nrank, nrun=40, seed=10000)
w <- basis(nmf.res)
colnames(w) = paste0('NMF', seq(1,nrank))
write.csv(w, file = '/bigdata/zlin/PanCancer_ICI/tables/nmf_w.csv')
w_ = w  |>  t()  |>  scale()
for (i in 1:3) {
  print(rownames(w_)[i])
  print(sort(w_[i,], decreasing = T) |> head(n=15))}
pdf('/bigdata/zlin/PanCancer_ICI/figures/NMF/ht.pdf', height = 3, width = 12)
Heatmap(w_, width = 11, height = 4, column_km = nrank, cluster_rows = F, clustering_method_columns = 'single', row_names_side = "left",
        col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("#154999", "white", "#CF0034")), na_col = 'lightgray',
        heatmap_legend_param = list(title = "Loading \n(z-score)", at = c(-0.5, 0, 1.5), labels = c("Min", "", "Max")),
        rect_gp = gpar(col = "gray", lwd = 0.5),
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 10)) 
dev.off()
h <- coef(nmf.res) |> t() |> as.data.frame()
colnames(h) = paste0('NMF', seq(1,nrank)) 
write.csv(h, file = '/bigdata/zlin/PanCancer_ICI/tables/nmf_h.csv')
h$patient <- str_replace(rownames(h), '_Pre|_On', '')
h$response <- pt_df$response[match(h$patient, pt_df$patient)]
h$dataset <- pt_df$dataset[match(h$patient, pt_df$patient)]
h$matrice <- pt_df$res_metric[match(h$patient, pt_df$patient)]

h$time_point <- 'Pre'
h$time_point[str_detect(rownames(h), 'On')] <- 'On'
h$sample <- paste0(h$patient, '_', h$time_point)
rownames(h) <- NULL

h |> pivot_longer(names_to = 'NMF', values_to = 'value', col = contains('NMF')) |> 
  select(-sample) |> 
  pivot_wider(names_from = 'time_point', values_from = 'value') |> 
  mutate(dynamic=On-Pre) |> 
  ggplot(aes(x=response, y=dynamic)) +
  geom_violin() +
  geom_boxplot(aes(col=response),alpha=.5, width=.6, lwd=.8, outlier.color = NA) +
  geom_point(aes(col=response), position = position_jitter(width = .15)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  facet_wrap(~NMF, nrow=1) +
  stat_compare_means(method = 't.test', label.y = 0.05) +
  scale_color_manual(values = c('RE'='#A2001F', 'NR'='#005D89'), name='Response') +
  theme_pubr() + ylab('Dynamic \n(NMF score )') + xlab('') +
  theme(plot.title = element_text(face = "bold", size=16),
        strip.text = element_text(size = 12, color='black'),
        strip.background =element_rect(fill=NA, color = NA, size = .5), 
        axis.text.x = element_text(size=12, angle = 0, hjust = .5, vjust=.5),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 14),
        axis.line = element_line(size = .5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = 'right',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 12, face = 'bold'),
        plot.margin = margin(t = 0, l=0)) 
ggsave('/bigdata/zlin/PanCancer_ICI/figures/NMF/dyn_score.pdf', width = 8, height = 3)

h_tme <- h |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'CRC_Li', 'PCa_Hawley', 'HNSC_Franken')) |> 
  select(NMF1, NMF2, NMF3, sample) |> 
  column_to_rownames(var = 'sample')
mtx_nonimmune <- meta_int |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC_Yost', 'CRC_Li', 'CRC_Chen', 'PCa_Hawley', 'HNSC_Franken'), 
         !patient %in% c('BCC/SCC_Yost_su009', 'BCC/SCC_Yost_su011', 'BCC/SCC_Yost_su012', 'BCC/SCC_Yost_su014'),
         celltype_r2 %in% nonimmune) |> 
  select(sample,  freq_r2_comp, celltype_r2) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = celltype_r2, values_fill = 0) |> 
  column_to_rownames(var = 'sample') 
common_sample <- intersect(rownames(h_tme), rownames(mtx_nonimmune))
df <- cbind(h_tme[common_sample,], mtx_nonimmune[common_sample,])
cor.res <- cor(df) |> as.data.frame()
cor.res <- cor.res[1:3, -c(1:3)] 

cor.test.res <- cor.mtest(df, conf.level = 0.95) 
cor.test.res <- cor.test.res$p[1:3, -c(1:3)]

pdf('/bigdata/zlin/PanCancer_ICI/figures/NMF/cor_nonimmune.pdf', width = 4, height = 2)
Heatmap(cor.res, name = "mat", col = circlize::colorRamp2(c(-0.3, 0, 0.3), c("#154999", "white", "#CF0034")),
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "left", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Pearson's \nCorrelation", 
                                    at = c(-0.8, 0, 0.8), 
                                    labels = c("0.8", "0", "0.8"),
                                    legend_direction = "horizontal", 
                                    legend_width = unit(1.5, "cm"), 
                                    legend_side = 'bottom',
                                    title_position = "topcenter"),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(cor.res)*unit(6, "mm"), 
        height = nrow(cor.res)*unit(6, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(cor.test.res[i, j]) & cor.test.res[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(cor.test.res[i, j]) & cor.test.res[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(cor.test.res[i, j]) & cor.test.res[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
dev.off()


