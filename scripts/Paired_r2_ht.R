rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggplot2','gtools','ComplexHeatmap','RColorBrewer','tibble','rstatix','purrr','effsize','ggpubr')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

meta_int <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_int.csv')

paired_t <- function(meta_combi, cohort, condition, values){
  test_list <- list()
  df_check <- meta_combi |>
    filter(dataset %in% cohort) |>
    select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior, modality) |>
    distinct(celltype_r2, patient, time_point, .keep_all = T) |>
    pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0)
  get_subtype <- function(value){
    df_check |>
      filter(eval(parse(text = sprintf("%s == '%s'", condition, value)))) |>
      group_by(celltype_r2) |>
      summarise(count = n()) |>
      filter(count >= 3) |>
      pull(celltype_r2)
  }
  subtypes <- values |> map(get_subtype) |> reduce(intersect)
  check_and_add_columns <- function(data, columns) {
    for (column in columns) {
      if (!column %in% names(data)) {
        data[[column]] <- 0
      }
    }
    return(data)
  }
  # paired t test with multiple hypothosis testing
  for (value in values){
    print(value)
    # calculate for each value in the condition
    df <- filter(meta_combi, dataset %in% cohort, eval(parse(text = condition)) == value) |> 
      select(freq_r2_comp, time_point, patient, celltype_r2) |> 
      distinct(celltype_r2, patient, time_point, .keep_all = T)
    dep_comp <- lapply(subtypes, function(subtype){
      print(subtype)
      df_sub <- df |> 
        filter(celltype_r2 == subtype) |>
        pivot_wider(names_from = time_point, values_from = freq_r2_comp, values_fill = 0)
      # Check for columns 'Pre' and 'Post'
      df_sub <- check_and_add_columns(df_sub, c("Pre", "Post"))
      t <- t.test(df_sub$Post, df_sub$Pre, paired = T, alternative = 'two.sided')
      cohen_d <- cohen.d(df_sub$Post, df_sub$Pre, paired = TRUE, hedges.correction = TRUE)
      comp <- c(subtype, t$statistic, t$p.value, cohen_d$estimate)
      return(comp)
    })
    results <- do.call(rbind, dep_comp) |> data.frame()
    rownames(results) <- subtypes
    colnames(results) <- c('celltype','t_score', 'pvalue','effectsize')
    results$fdr <- p.adjust(results$pvalue, method = 'fdr', n = nrow(results))
    test_list[[which(values == value)]] <- results 
  }
  names(test_list) <- values
  df <- reduce(test_list, full_join, by = "celltype", suffix = paste0('_', names(test_list)))
  df <- df |> column_to_rownames(var = 'celltype')
  return(df)
}


res_comp <- paired_t(meta_int, cohort = unique(meta_int$dataset), condition = 'modality', values = c('Mono','Dual')) 
celltypes <- rownames(res_comp)
pdf('/bigdata/zlin/Melanoma_meta/figures/Change/scatter_modality.pdf', height = 5, width = 6.5)
res_comp |> 
  apply(2, as.numeric) |> 
  as.data.frame() |> 
  mutate(celltype_label = ifelse(abs(effectsize_Mono - effectsize_Dual) > 0.2 & 
                                   (abs(effectsize_Mono) > 0.2 & abs(effectsize_Dual) > 0.2), celltypes,
                                 ifelse((pvalue_Mono < 0.05 & pvalue_Dual>= 0.05)|(pvalue_Mono >= 0.05 & pvalue_Dual < 0.05), celltypes,
                                        ifelse(pvalue_Mono < 0.05 & pvalue_Dual<0.05, celltypes, ''))),
         Significance = ifelse(pvalue_Mono < 0.05 & pvalue_Dual<0.05, 'Significant in Both', 
                               ifelse(pvalue_Mono < 0.05 & pvalue_Dual>= 0.05, 'Significant in Mono', 
                                      ifelse(pvalue_Mono >= 0.05 & pvalue_Dual < 0.05, 'Significant in Dual', 'Non-significant')))) |>
  mutate(Significance = factor(Significance, levels = c('Significant in Mono', 'Significant in Dual', 'Significant in Both', 'Non-significant'))) |> 
  ggplot(aes(x = effectsize_Mono, y = effectsize_Dual)) + 
  geom_point(aes(color = Significance), size = 2) + 
  scale_color_manual(values = c('#E41A1C', '#377EB8', '#4DAF4A', '#999999'), name = 'Significance\n(Paired t-test)') +
  geom_abline(intercept = 0.2, slope = 1, linetype = "dashed", size = 0.1) +
  geom_abline(intercept = -0.2, slope = 1, linetype = "dashed", size = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 0.1) +
  geom_hline(yintercept = -0.5, linetype = "dashed", size = 0.1) +
  geom_vline(xintercept = -0.5, linetype = "dashed", size = 0.1) +
  theme_classic() + xlim(c(-0.9, 0.9)) + ylim(c(-0.9, 0.9)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "lightgray", size = 0.05, linetype = "dashed"),
        panel.grid.minor = element_line(color = "lightgray", size = 0.05, linetype = "dashed")) +
  geom_text_repel(aes(label = celltype_label), box.padding = 0.8, size = 2.5, max.overlaps = 50) + 
  xlab("Mono\n(Cohen's d)") +
  ylab("Dual\n(Cohen's d)") + ggtitle('Modality')
y_text = 0.033
x1 = 0.15
x2 = 0.7
grid.segments(
  x = unit(x1 + 0.05, "npc"),
  x1 = unit(x1 - 0.05, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x2 - 0.05, "npc"),
  x1 = unit(x2 + 0.05, "npc"),
  y = unit(y_text+0.03, "npc"),
  y1 = unit(y_text+0.03, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Down", x = unit(x1, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
grid.text("Up", x = unit(x2, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))

x_text = 0.05
y1 = 0.17
y2 = 0.87
grid.segments(
  x = unit(x_text + 0.02, "npc"),
  x1 = unit(x_text + 0.02, "npc"),
  y = unit(y1 + 0.05, "npc"),
  y1 = unit(y1 - 0.05, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x_text + 0.02, "npc"),
  x1 = unit(x_text + 0.02, "npc"),
  y = unit(y2 - 0.05, "npc"),
  y1 = unit(y2 + 0.05, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Down", y = unit(y1, "npc"), x = unit(x_text, "npc"), gp = gpar(fontsize = 8), rot = 90)
grid.text("Up", y = unit(y2, "npc"), x = unit(x_text, "npc"), gp = gpar(fontsize = 8), rot = 90)
dev.off()




res <- res_comp |> 
  filter(pvalue_Mono<0.05|pvalue_Dual<0.05) |> 
  arrange(desc(effectsize_Mono))
mat_es <- res[,str_detect(colnames(res), 'effectsize')] |> sapply(as.numeric)
rownames(mat_es) <- rownames(res)
colnames(mat_es) <- c('Mono','Dual')
mat_sig <- res[,str_detect(colnames(res), 'pvalue')] |> sapply(as.numeric)
pdf('/bigdata/zlin/Melanoma_meta/figures/Change/ht_monovsdual.pdf', height = 8, width = 3)
Heatmap(mat_es, name = "mat", col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), column_title = 'Treatment Modality',
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "left", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Effect Size \n(Cohen's d)", at = c(-0.8, 0, 0.8), labels = c("0.8", "0", "0.8")),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_es)*unit(6, "mm"), 
        height = nrow(mat_es)*unit(6, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
dev.off()

# meta_int |> 
#   filter(celltype_r2 == 'Mono_CD16') |> 
#   distinct(sample, .keep_all = T) |> 
#   group_by(patient) |> 
#   mutate(pt_count = n()) |> 
#   filter(pt_count == 2) |> 
#   ungroup() |> 
#   mutate(time_point = factor(time_point, levels = c('Pre', 'Post')),
#          modality = factor(modality, levels = c('Mono', 'Dual'))) |> 
#   ggplot(aes(x = time_point, y = freq_r2_comp)) +
#   geom_boxplot(aes(fill = time_point)) +
#   # geom_point(aes(color = response)) +
#   geom_point() +
#   geom_line(aes(group = patient), color = "gray", linetype = 'dotted') +
#   facet_wrap(.~modality) + 
#   stat_compare_means(method = 't.test', paired = T, label.x = 1.3, label.y = 0.15) +
#   scale_fill_manual(values = c('#CF0034','#154999')) +
#   xlab('Time Point') + ylab('Freqency') + ggtitle('Mono_CD16') +
#   theme_classic() + theme(legend.position='none')
# ggsave('/bigdata/zlin/Melanoma_meta/figures/Change//boxplot_Mono_CD16.pdf', height = 4, width = 6)

res_comp <- paired_t(meta_int, cohort = 'HNSC_Franken' , condition = 'treatment', values = c('aPDL1','aPDL1+CTLA4')) 
res <- res_comp |> 
  filter(pvalue_aPDL1<0.05|`pvalue_aPDL1+CTLA4`<0.05) |> 
  arrange(desc(effectsize_aPDL1))
mat_es <- res[,str_detect(colnames(res), 'effectsize')] |> sapply(as.numeric)
rownames(mat_es) <- rownames(res)
colnames(mat_es) <- c('aPDL1','aPDL1+CTLA4')
mat_sig <- res[,str_detect(colnames(res), 'pvalue')] |> sapply(as.numeric)
pdf('/bigdata/zlin/Melanoma_meta/figures/Change/ht_hnsc_pdl1.pdf', height = 6, width = 3)
Heatmap(mat_es, name = "mat", col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), column_title = 'aPDL1',
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "left", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Effect Size \n(Cohen's d)", at = c(-0.8, 0, 0.8), labels = c("0.8", "0", "0.8")),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_es)*unit(6, "mm"), 
        height = nrow(mat_es)*unit(6, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
dev.off()

pdf('/bigdata/zlin/Melanoma_meta/figures/Change/ht_hnsc_pdl1_ls.pdf', height = 3, width = 10)
Heatmap(t(mat_es), name = "mat", col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), column_title = '',
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "left", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Effect Size \n(Cohen's d)", at = c(-0.8, 0, 0.8), labels = c("0.8", "0", "0.8"),
                                    legend_direction = "horizontal", legend_side = 'bottom', title_position = "topcenter"),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_es)*unit(6, "cm"), 
        height = nrow(mat_es)*unit(1, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(t(mat_sig)[i, j]) & t(mat_sig)[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(t(mat_sig)[i, j]) & t(mat_sig)[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(t(mat_sig)[i, j]) & t(mat_sig)[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
dev.off()



res_comp <- paired_t(meta_int, cohort = c('HNSC_Luoma','HNSC_IMCISION') , condition = 'treatment', values = c('aPD1','aPD1+CTLA4')) 
res <- res_comp |> 
  filter(pvalue_aPD1 < 0.05|`pvalue_aPD1+CTLA4` < 0.05) |> 
  arrange(desc(effectsize_aPD1))
mat_es <- res[,str_detect(colnames(res), 'effectsize')] |> sapply(as.numeric)
rownames(mat_es) <- rownames(res)
colnames(mat_es) <- c('aPD1','aPD1+CTLA4')
mat_sig <- res[,str_detect(colnames(res), 'pvalue')] |> sapply(as.numeric)
pdf('/bigdata/zlin/Melanoma_meta/figures/Change/ht_hnsc_pd1.pdf', height = 6, width = 3)
Heatmap(mat_es, name = "mat", col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), column_title = 'aPD1',
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "left", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Effect Size \n(Cohen's d)", at = c(-0.8, 0, 0.8), labels = c("0.8", "0", "0.8")),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_es)*unit(6, "mm"), 
        height = nrow(mat_es)*unit(6, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
dev.off()

pdf('/bigdata/zlin/Melanoma_meta/figures/Change/ht_hnsc_pd1_ls.pdf', height = 3, width = 10)
Heatmap(t(mat_es), name = "mat", col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), column_title = '',
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "left", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Effect Size \n(Cohen's d)", at = c(-0.8, 0, 0.8), labels = c("0.8", "0", "0.8"),
                                    legend_direction = "horizontal", legend_side = 'bottom', title_position = "topcenter"),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_es)*unit(6, "cm"), 
        height = nrow(mat_es)*unit(0.8, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(t(mat_sig)[i, j]) & t(mat_sig)[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(t(mat_sig)[i, j]) & t(mat_sig)[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(t(mat_sig)[i, j]) & t(mat_sig)[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })
dev.off()

res_comp <- paired_t(meta_int, cohort = 'SKCM_Becker' , condition = 'treatment', values = c('aPD1','aPD1+CTLA4')) 
res <- res_comp |> filter(pvalue_aPD1 < 0.05|`pvalue_aPD1+CTLA4` < 0.05) 
mat_es <- res[,str_detect(colnames(res), 'effectsize')] |> sapply(as.numeric)
rownames(mat_es) <- rownames(res)
colnames(mat_es) <- c('aPD1','aPD1+CTLA4')
mat_sig <- res[,str_detect(colnames(res), 'pvalue')] |> sapply(as.numeric)

Heatmap(mat_es, name = "mat", col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), column_title = 'SKCM_Becker',
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "right", column_names_side = 'bottom',
        heatmap_legend_param = list(title = "Effect Size \n(Cohen's d)", at = c(-0.8, 0.8), labels = c("Low","High")),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_es)*unit(7, "mm"), 
        height = nrow(mat_es)*unit(7, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=15))
          }
          else if(!is.na(mat_sig[i, j]) & mat_sig[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=15))
          }
        })



