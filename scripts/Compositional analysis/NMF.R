rm(list=ls())
pkgs <- c('qs2','tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','tibble','pheatmap','MetBrewer','viridis','ComplexHeatmap','colorRamp2','corrr','ggnewscale','NMF','ggpubr','corrplot','rstatix')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
setwd("/home/zlin/workspace/PanCancer_ICI")

metadata <- read.csv('tables/meta_all.csv') 
metadata$freq_r2_comp[metadata$celltype_r2 == 'Malignant(CNA+)'] <- metadata$freq_r2[metadata$celltype_r2 == 'Malignant(CNA+)']
patient_info <- metadata |> distinct(patient, .keep_all = T)
pt <- metadata |>
  distinct(sample, .keep_all = TRUE) |>
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |>
  pull(patient)
source('scripts/Celltype_classification.R')
unwanted_celltypes <- c('Cycling T','Cycling NK', "GCB-cycling", "PC-cycling", 'Cycling myeloids')
mtx_immune <- metadata |> 
  dplyr::filter(celltype_r2 %in% c(t_nk, mye, bplasma),
         subset %in% c('All TME', 'CD45+sorted'),
         cohort != 'SKCM_this study') |> 
  dplyr::filter(!celltype_r2 %in% unwanted_celltypes) |>
  select(sample,  freq_r2_comp, celltype_r2) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = celltype_r2, values_fill = 0) |> 
  column_to_rownames(var = 'sample') |> t()
res.nmf <- nmf(mtx_immune, rank = 2:8, nrun = 30, seed = 1234)
qs_save(res.nmf,'data/res_nmf_seed_1234.qs2')
res.nmf <- qs_read('data/res_nmf_seed_1234.qs2')
pdf('figures/NMF/rank_survey.pdf', height = 6, width =8)
plot(res.nmf)
dev.off()

nrank=3
nmf.res = nmf(mtx_immune, rank=nrank, nrun=30, seed=1234)
qs_save(nmf.res,'data/res_nmf_k3_seed_1234.qs2')
nmf.res <- qs_read('data/res_nmf_k3_seed_1234.qs2')
w <- basis(nmf.res)
colnames(w) = paste0('NMF', seq(1,nrank))
w_ = w |> 
  t() |>  
  scale()
w_df <- w_ |> t() |> data.frame() |> rownames_to_column('CellType')
celltypes <- list(
  NMF1 = w_df |>
    mutate(NMF1_diff = NMF1 - (NMF2 + NMF3)) |>
    filter(NMF1_diff > 0) |>
    arrange(desc(NMF1_diff)) |>
    pull(CellType) |> head(n = 10),
  
  NMF2 = w_df |>
    mutate(NMF2_diff = NMF2 - (NMF1 + NMF3)) |>
    filter(NMF2_diff > 0) |>
    arrange(desc(NMF2_diff)) |>
    pull(CellType) |> head(n = 10),
  
  NMF3 = w_df |>
    mutate(NMF3_diff = NMF3 - (NMF1 + NMF2)) |>
    filter(NMF3_diff > 0) |>
    arrange(desc(NMF3_diff)) |>
    pull(CellType) |> head(n = 10)
)
lapply(celltypes, function(x) { print(x) })

pdf('figures/NMF/ht_rank3.pdf', height = 2.5, width = 12)
Heatmap(w_, width = 11, height = 4,  cluster_rows = F, 
        row_names_side = "left",
        column_names_rot = 45,
        col = circlize::colorRamp2(c(-1.154, 0, 1.154), c("#154999", "white", "#CF0034")), na_col = 'lightgray',
        heatmap_legend_param = list(title = "Loading \n(scaled)", at = c(-1.154, 0, 1.154), labels = c("Min", "", "Max")),
        rect_gp = gpar(col = "lightgray", lwd = 0.5),
        column_names_gp = gpar(fontsize = 9),
        row_names_gp = gpar(fontsize = 10)) 
dev.off()

h <- coef(nmf.res) |> t() |> as.data.frame()
colnames(h) = paste0('NMF', seq(1,nrank))
h$patient <- str_replace(rownames(h), '_Pre|_On|_Post', '')
h$time_point <- str_split(rownames(h), '_', simplify = T)[,4]
h$time_point[h$time_point == 'Post'] <- 'On'
h$response <- patient_info$response[match(h$patient, patient_info$patient)]
h$response <- factor(h$response, levels = c('R','NR'))
h <- h |> drop_na()
h$cohort <- patient_info$cohort[match(h$patient, patient_info$patient)]
h$subtype <- patient_info$subtype[match(h$patient, patient_info$patient)]
h$time_point[h$cohort == 'RCC_Bi'] <- 'ICI_exposed'
h$tx_status <- 'Baseline'
h$tx_status[h$time_point %in% c('On','ICI_exposed')] <- 'Treated'
h$sample <- paste0(h$patient, '_', h$tx_status)
h$status_response <- paste0(h$tx_status, '.', h$response)
h$status_response <- factor(h$status_response, levels = c('Baseline.R', 'Baseline.NR', 'Treated.R', 'Treated.NR' ))
rownames(h) <- NULL

stat.test <- h |> select(-sample) |> 
  pivot_longer(cols = contains('NMF'), names_to = 'NMF', values_to = 'score') |> 
  drop_na() |> 
  group_by(NMF) |> 
  wilcox_test(score ~ status_response) |> 
  add_significance() |> 
  filter(!((group1 == 'Baseline.R'& group2=='Treated.NR')|
           (group1 == 'Baseline.NR'& group2=='Treated.R'))) |>
  mutate(
    p.label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ NA_character_
    )
  ) |> 
  filter(p<0.1)
stat.test <- stat.test |> add_y_position(step.increase = 0.2, scales = 'free_y')
h_long <- h |> select(-sample) |> 
  pivot_longer(cols = contains('NMF'), names_to = 'NMF', values_to = 'score') |> 
  drop_na() 
h_long |> ggplot(aes(x=status_response, y=score)) +
  geom_boxplot(aes(fill=status_response), width=.6, lwd=.8, outlier.color = NA) +
  geom_point(position = position_jitter(width = .15), size = 0.3) +
  geom_hline(yintercept = 0, linetype='dashed') +
  facet_wrap(~NMF, nrow=1, scales = 'free_y') +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size = 5) +
  scale_fill_manual(values = met.brewer(name="Egypt",n=4)) +
  theme_pubr() + ylab('NMF score') + xlab('') +
  theme(plot.title = element_text(face = "bold", size=16),
        strip.text = element_text(size = 12, color='black'),
        strip.background =element_rect(fill=NA, color = NA, size = .5), 
        axis.text.x = element_text(size=10, angle = 30, hjust = .5, vjust=.5),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 14),
        axis.line = element_line(size = .5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = 'none') +
  coord_cartesian(ylim = ) 
ggsave('figures/NMF/response.pdf', height = 4, width = 10)

# validation
mtx_immune_becker <- metadata |> 
  dplyr::filter(celltype_r2 %in% c(t_nk, mye, bplasma),
                subset %in% c('All TME', 'CD45+sorted'),
                cohort == 'SKCM_this study') |> 
  dplyr::filter(!celltype_r2 %in% unwanted_celltypes) |>
  select(sample,  freq_r2_comp, celltype_r2) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = celltype_r2, values_fill = 0) |> 
  column_to_rownames(var = 'sample')
common_ct <- intersect(rownames(w), colnames(mtx_immune_becker))
score <- as.matrix(mtx_immune_becker[,common_ct]) %*% as.matrix(w[common_ct,])
score <- data.frame(score)
score$patient <- str_replace(rownames(score), '_Pre|_On|_Post', '')
score$time_point <- str_split(rownames(score), '_', simplify = T)[,4]
score$response <- patient_info$response[match(score$patient, patient_info$patient)]
score$response <- factor(score$response, levels = c('R','NR'))
score$cohort <- patient_info$cohort[match(score$patient, patient_info$patient)]
score$subtype <- patient_info$subtype[match(score$patient, patient_info$patient)]
score$tx_status <- 'Baseline'
score$tx_status[score$time_point %in% c('On','ICI_exposed')] <- 'Treated'
score$sample <- paste0(score$patient, '_', score$tx_status)
score$status_response <- paste0(score$tx_status, '.', score$response)
score$status_response <- factor(score$status_response, levels = c('Baseline.R', 'Baseline.NR', 'Treated.R', 'Treated.NR' ))
rownames(score) <- NULL
stat.test <- score |> select(-sample) |> 
  pivot_longer(cols = contains('NMF'), names_to = 'NMF', values_to = 'score') |> 
  drop_na() |> 
  group_by(NMF) |> 
  wilcox_test(score ~ status_response) |> 
  add_significance() |> 
  filter(!((group1 == 'Baseline.R'& group2=='Treated.NR')|
             (group1 == 'Baseline.NR'& group2=='Treated.R'))) |>
  mutate(
    p.label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ NA_character_
    )
  )

stat.test <- stat.test |> add_y_position(step.increase = 0.2, scales = 'free_y')
score_long <- score |> select(-sample) |> 
  pivot_longer(cols = contains('NMF'), names_to = 'NMF', values_to = 'score') |> 
  drop_na() 
score_long |> ggplot(aes(x=status_response, y=score)) +
  # geom_violin() +
  geom_boxplot(aes(fill=status_response), width=.6, lwd=.8, outlier.color = NA) +
  geom_point(position = position_jitter(width = .15), size = 0.8) +
  geom_hline(yintercept = 0, linetype='dashed') +
  facet_wrap(~NMF, nrow=1, scales = 'free_y') +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0.01, size = 5) +
  scale_fill_manual(values = met.brewer(name="Egypt",n=4)) +
  theme_pubr() + ylab('NMF score') + xlab('') +
  theme(plot.title = element_text(face = "bold", size=16),
        strip.text = element_text(size = 12, color='black'),
        strip.background =element_rect(fill=NA, color = NA, size = .5), 
        axis.text.x = element_text(size=10, angle = 30, hjust = .5, vjust=.5),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 14),
        axis.line = element_line(size = .5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = 'none') 
ggsave('figures/NMF/response_validation.pdf', height = 4, width = 10)

# dynamics
h |> filter(patient %in% pt) |> 
  select(-c(status_response, sample, time_point)) |> 
  pivot_longer(cols = contains('NMF'), names_to = 'NMF', values_to = 'score') |> 
  pivot_wider(names_from = 'tx_status', values_from = 'score') |> 
  drop_na() |> 
  mutate(dynamic=Treated-Baseline) |> 
  ggplot(aes(x=response, y=dynamic)) +
  geom_boxplot(aes(fill=response), width=.6, lwd=.8, outlier.color = NA) +
  geom_point(position = position_jitter(width = .15), size = 0.3) +
  geom_hline(yintercept = 0, linetype='dashed') +
  facet_wrap(~NMF, nrow=1, scales = 'free_y') +
  stat_compare_means(method = 'wilcox', label.x = 1) +
  scale_fill_manual(values = c('R'='#A2001F', 'NR'='#005D89'), name='Response') +
  theme_pubr() + ylab('delta NMF score \n(On-Pre)') + xlab('') +
  theme(plot.title = element_text(face = "bold", size=16),
        strip.text = element_text(size = 12, color='black'),
        strip.background =element_rect(fill=NA, color = NA, size = .5), 
        axis.text.x = element_text(size=12, angle = 0, hjust = .5, vjust=.5),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 14),
        axis.line = element_line(size = .5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = 'none',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 12, face = 'bold')) 
ggsave('figures/NMF/delta.pdf', height = 3.5, width = 12)
# validation
score |> filter(patient %in% pt) |> 
  select(-c(status_response, sample, time_point)) |> 
  pivot_longer(cols = contains('NMF'), names_to = 'NMF', values_to = 'score') |> 
  pivot_wider(names_from = 'tx_status', values_from = 'score') |> 
  drop_na() |> 
  mutate(dynamic=Treated-Baseline) |> 
  ggplot(aes(x=response, y=dynamic)) +
  geom_boxplot(aes(fill=response), width=.6, lwd=.8, outlier.color = NA) +
  geom_point(position = position_jitter(width = .15), size = 0.8) +
  geom_hline(yintercept = 0, linetype='dashed') +
  facet_wrap(~NMF, nrow=1, scales = 'free_y') +
  stat_compare_means(method = 'wilcox', label.x = 1) +
  scale_fill_manual(values = c('R'='#A2001F', 'NR'='#005D89'), name='Response') +
  theme_pubr() + ylab('delta NMF score \n(On-Pre)') + xlab('') +
  theme(plot.title = element_text(face = "bold", size=16),
        strip.text = element_text(size = 12, color='black'),
        strip.background =element_rect(fill=NA, color = NA, size = .5), 
        axis.text.x = element_text(size=12, angle = 0, hjust = .5, vjust=.5),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 14),
        axis.line = element_line(size = .5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = 'none',
        legend.key.size = unit(1.4,'line'),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 12, face = 'bold')) 
ggsave('figures/NMF/delta_validation.pdf', height = 3.5, width = 12)


h <- coef(nmf.res) |> t() |> as.data.frame()
colnames(h) = paste0('NMF', seq(1,nrank)) 
h <- h[!str_detect(rownames(h), 'RCC_Bi'),]
h$patient <- str_replace(rownames(h), '_Pre|_On|_Post', '')
h$time_point <- str_split(rownames(h), '_', simplify = T)[,4]
h$response <- patient_info$response[match(h$patient, patient_info$patient)]
h$response <- factor(h$response, levels = c('R','NR'))
h <- h |> filter(response != 'NE')
h$cohort <- patient_info$cohort[match(h$patient, patient_info$patient)]
h$subtype <- patient_info$subtype[match(h$patient, patient_info$patient)]
h$sample <- paste0(h$patient, '_', h$time_point)
h$time_point_response <- paste0(h$time_point, '.', h$response)
h$time_point_response <- factor(h$time_point_response, levels = c('Pre.R','Pre.NR','On.R','On.NR'))
h <- h |> filter(subtype != 'PCa')
rownames(h) <- NULL

stat.test <- h |> select(-sample) |> 
  pivot_longer(cols = contains('NMF'), names_to = 'NMF', values_to = 'value') |> 
  drop_na() |> 
  group_by(NMF, subtype) |> 
  wilcox_test(value ~ time_point_response) |> 
  filter(!((group1 == 'Pre.R'&group2=='On.NR')|
             (group1 == 'Pre.NR'&group2=='On.R'))) |> 
  mutate(p.signif = ifelse(p<0.05, '*', 'ns'))
stat.test <- stat.test |> add_y_position(step.increase = 0.3) 
h |> select(-sample) |> 
  pivot_longer(cols = contains('NMF'), names_to = 'NMF', values_to = 'value') |> 
  drop_na() |> 
  ggplot(aes(x=time_point_response, y=value)) +
  geom_violin() +
  geom_boxplot(aes(col=time_point_response),alpha=.5, width=.6, lwd=.8, outlier.color = NA) +
  geom_point(aes(col=time_point_response), position = position_jitter(width = .15)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  facet_wrap(~ NMF + subtype, nrow=3) +
  stat_pvalue_manual(stat.test, label = "p.signif", tip.length = 0.01, size = 5) +
  scale_color_manual(values = met.brewer(name="Egypt",n=4)) +
  theme_pubr() + ylab('NMF score') + xlab('') +
  theme(plot.title = element_text(face = "bold", size=16),
        strip.text = element_text(size = 12, color='black'),
        strip.background =element_rect(fill=NA, color = NA, size = .5), 
        axis.text.x = element_text(size=10, angle = 30, hjust = .5, vjust=.5),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 14),
        axis.line = element_line(size = .5),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = 'none') 
ggsave('figures/NMF/response_sub.pdf', height = 8, width = 12)

h <- coef(nmf.res) |> t() |> as.data.frame()
colnames(h) = paste0('NMF', seq(1,nrank)) 
sample_tme <- metadata |> 
  distinct(sample, .keep_all = T) |> 
  filter(subset == 'All TME') |> 
  select(sample) |> pull()
common <- intersect(rownames(h), sample_tme)
h <- h[common,]
mtx <- metadata |> 
  filter(sample %in% common) |> 
  distinct(sample, celltype_r2, .keep_all = T)
# unique(mtx$celltype_r2) |> setdiff(c(t_nk, mye, bplasma, nonimmune))
mtx$freq_r2_comp[mtx$celltype_r2 == "Malignant(CNA+)"] <- mtx$freq_cna[mtx$celltype_r2 == "Malignant(CNA+)"]
mtx <- mtx |> 
  filter(celltype_r2 %in% c(nonimmune, "Malignant(CNA+)")) |> 
  select(sample, celltype_r2, freq_r2_comp) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = celltype_r2, values_fill = 0)
til <- metadata |> 
  filter(sample %in% mtx$sample) |> 
  distinct(sample, .keep_all = T) |> 
  select(freq_t, sample)
names(til)[1] <- 'T cell infiltration'
mtx <- mtx |> 
  merge(til, by='sample') 
mtx <- cbind(h[mtx$sample,], mtx)
mtx$sample <- NULL
# add mean expression of PDL12
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'RCC_Bi','HCC_Guo','HCC_Ma')
list_df <- lapply(datasets, function(dataset){ 
  print(dataset)
  df <- read.csv(paste0('data/', dataset, '/avg_expr_pdl12.csv'), row.names = 'X')
  rownames(df) <- str_replace(rownames(df), 'Post', 'On')
  return(df)
})
df <- do.call(rbind, list_df)
rownames(mtx) <- str_replace(rownames(mtx), '_ICI_exposed|_NoICI', '')
common <- intersect(rownames(df), rownames(mtx))
df <- df[common,]
mtx <- mtx[common,]
mtx <- cbind(mtx, 'CD274(myeloids)'=df$CD274)
mtx <- mtx[,-which(names(mtx) == 'Cycling non-immune')]

cor.mat <- cor_mat(mtx)
cor.mat.p <- cor.mat |>  cor_get_pval()
cor.mat <- cor.mat |> column_to_rownames(var = 'rowname')
cor.mat <- cor.mat[ -c(1:3),1:3] 
cor.mat.p <- cor.mat.p |> column_to_rownames(var = 'rowname')
cor.mat.p <- cor.mat.p[ -c(1:3),1:3] 

# ct_order <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein",
#               "Pericytes","SMC", "Myofibroblasts", "CAF_SFRP2", 
#               "CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap", "Malignant(CNA+)", "T cell infiltration")
pdf('figures/NMF/ht_nonimmune.pdf', height =6, width = 4)
Heatmap(cor.mat, name = "mat", col = circlize::colorRamp2(c(-0.6, 0, 0.6), c("#154999", "white", "#CF0034")),
        cluster_rows = T, cluster_columns = T, column_names_rot =45, row_names_side = "left", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Pearson's \nCorrelation", 
                                    at = c(-0.6, 0, 0.6), 
                                    labels = c("0.6", "0", "0.6"),
                                    legend_direction = "horizontal", 
                                    title_position = "topcenter"),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(cor.mat)*unit(6, "mm"), 
        height = nrow(cor.mat)*unit(6, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          gb = textGrob("*", gp = gpar(fontsize=14))
          gb_h = convertHeight(grobHeight(gb), "mm")
          gb_w = convertWidth(grobWidth(gb), "mm")
          if(!is.na(cor.mat.p[i, j]) & cor.mat.p[i, j] < 0.001) {
            grid.text('***', x, y - gb_h*0.5 + gb_w*0.4,gp=gpar(fontsize=15))
          }
          else if(!is.na(cor.mat.p[i, j]) & cor.mat.p[i, j] < 0.01) {
            grid.text('**', x, y - gb_h*0.5 + gb_w*0.4,gp=gpar(fontsize=15))
          }
          else if(!is.na(cor.mat.p[i, j]) & cor.mat.p[i, j] < 0.05) {
            grid.text('*', x, y - gb_h*0.5 + gb_w*0.4,gp=gpar(fontsize=15))
          }
        })
dev.off()



