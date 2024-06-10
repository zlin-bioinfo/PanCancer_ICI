rm(list=ls())
pkgs <- c('tidyr','dplyr','stringr','ggsci','patchwork','ggplot2','gtools','dittoSeq','RColorBrewer','ggpubr','tibble','epitools', 'cowplot','forcats','rstatix','pheatmap','janitor','ggmosaic','gridExtra','ggalluvial','forcats')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Frequency by composition
# main level
meta_combi <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_int.csv') |> select(!X)
head(meta_combi)
freq_filt <- meta_combi |> 
  select(patient, time_point, celltype_main, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
  distinct(patient, time_point, celltype_main, .keep_all = T) |> 
  pivot_wider(values_from = count_r2, names_from = time_point, values_fill = 0) |> 
  filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5))
freq_filt$pt_main <- paste0(freq_filt$patient, '_', freq_filt$celltype_main)
freq_wide <- meta_combi |> 
  select(patient, time_point, celltype_main, interval, cancertype, response, res_metric, treatment, prior, modality, freq_main, dataset, component) |> 
  distinct(patient, time_point, celltype_main, .keep_all = T) |> 
  pivot_wider(values_from = freq_main, names_from = time_point, values_fill = 0)
freq_wide$pt_main <- paste0(freq_wide$patient, '_', freq_wide$celltype_main)
freq_wide <- filter(freq_wide, pt_main %in% freq_filt$pt_main)

ggpaired(freq_wide, cond1 = "Pre", cond2 = "Post",
         fill = "condition", palette = "jco", ylab = "Freqency", xlab = "Time Point") + 
  facet_wrap(.~ celltype_main, scales = 'free') + 
  stat_compare_means(method = 't.test', paired = T) + theme_pubclean()
  
# deep level
freq_filt <- meta_combi |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = count_r2, names_from = time_point) |> 
  filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5))
freq_filt$pt_r2 <- paste0(freq_filt$patient, '_', freq_filt$celltype_r2)
freq_long <- meta_combi |> 
  select(patient, time_point, celltype_r2, celltype_main, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) 
freq_long$pt_r2 <- paste0(freq_long$patient, '_', freq_long$celltype_r2)
freq_long <- filter(freq_long, pt_r2 %in% freq_filt$pt_r2)
freq_long$celltype_cate <- freq_long$celltype_main
freq_long$celltype_cate <- ifelse(freq_long$celltype_cate %in% c("Mast", "pDC", "cDC", "Mono", "Macro"), "Myeloids", 
                                  ifelse(freq_long$celltype_cate %in% c("CD8+T", "NK"), "CD8+T/NK",
                                         ifelse(freq_long$celltype_cate %in% c("Plasma", "B"), "B/Plasma",
                                                ifelse(freq_long$celltype_cate == "CD4+T", "CD4+T", "Non-immune"))))

# freq_long$celltype_cate <- factor(freq_long$celltype_cate, levels = c("CD4+T", "CD8+T/NK", "Myeloids", "B/Plasma", "Non-immune"))
# freq_long$time_point <- factor(freq_long$time_point, levels = c('Pre', 'Post'))

# list_main <- lapply(unique(freq_long$celltype_cate), function(major){
#   print(major)
#   cellstates <- freq_long |> filter(celltype_cate == major) |> pull(celltype_r2) |> unique()
#   res_comp <- list()
#   for (i in 1:length(cellstates)){
#     print(cellstates[[i]])
#     pt_filt <- meta_combi |> 
#       filter(celltype_r2 == cellstates[[i]]) |> 
#       select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
#       distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
#       pivot_wider(values_from = count_r2, names_from = time_point) |> 
#       filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5)) |> 
#       mutate(int_cat = case_when(interval < 21 ~ "<3w",
#                                  interval >=21 ~ ">=3w")) |> 
#       group_by(int_cat) |> 
#       mutate(pt_count = n()) |> 
#       filter(pt_count >= 3) |> pull(patient)
#     if (length(pt_filt) < 3){
#       res_comp[[i]] <- c(cellstates[[i]], NA, NA)
#     } else {
#       effectsize <- freq_long |> 
#         filter(celltype_r2 == cellstates[[i]], patient %in% pt_filt) |> 
#         cohens_d(freq_r2_comp ~ time_point, paired = TRUE, ref.group = 'Post')
#       ttest <- freq_long |> 
#         filter(celltype_r2 == cellstates[[i]], patient %in% pt_filt) |> 
#         t_test(freq_r2_comp ~ time_point, paired = TRUE, ref.group = 'Post')
#         fdr <- 
#       res_comp[[i]] <- c(cellstates[[i]], effectsize$effsize, ttest$p)
#     }
#   }
#   res_comp <- do.call(rbind, res_comp) |> data.frame()
#   colnames(res_comp) <- c("Cell States", "Cohen's d", "Pvalue")
#   res_comp <- res_comp |> 
#     adjust_pvalue(p.col = 'Pvalue', method = 'fdr') |> 
#     dplyr::filter(Pvalue.adj < 0.05)
#   return(res_comp)
# })
# names(list_main) <- unique(freq_long$celltype_cate)
# cs_sig <- unname(unlist(lapply(list_main, function(df) {df$`Cell States`})))

list_cs <- lapply(unique(freq_long$celltype_r2), function(cellstate){
  print(cellstate)
  pt_filt <- meta_combi |> 
    filter(celltype_r2 == cellstate) |> 
    select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
    distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
    pivot_wider(values_from = count_r2, names_from = time_point) |> 
    filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5)) |> 
    mutate(int_cat = case_when(interval < 21 ~ "<3w",
                               interval >=21 ~ ">=3w")) |> 
    group_by(int_cat) |> 
    mutate(pt_count = n()) |> 
    filter(pt_count >= 3) |> pull(patient)
  if (length(pt_filt) > 0){
    sub_long <- freq_long |> 
      filter(celltype_r2 == cellstate, patient %in% pt_filt, response != 'NE') |> 
      mutate(time_point2 = case_when(time_point == "Pre" ~ "Pre",
                                     time_point == "Post" & interval < 21 ~ "<3w",
                                     time_point == "Post" & interval >= 21 ~ ">=3w"),
             time_point2 = fct_relevel(time_point2, "Pre", "<3w", ">=3w"),
             pct_r2_comp_log1p = log1p(freq_r2_comp*100))
    if (length(unique(sub_long$time_point2)) > 2){
      stat.test <- rbind(sub_long |> filter(interval < 21) |> t_test(pct_r2_comp_log1p ~ time_point2, ref.group = '<3w'), 
                         sub_long |> filter(interval > 21) |> t_test(pct_r2_comp_log1p ~ time_point2, ref.group = '>=3w'))
    } else {
      stat.test <- sub_long |> t_test(pct_r2_comp_log1p ~ time_point2)
      print(unique(sub_long$time_point2))
    }
    if (stat.test |> filter(p < 0.05) |> nrow()>0) {
      print(cellstate)
      return(cellstate)
    }
  }
})
cs_sig <- unlist(list_cs) |> na.omit()

palette_cancertype <- setNames(brewer.pal(length(unique(freq_long$cancertype)),"Set1"), unique(freq_long$cancertype))
list_boxplot <- lapply(cs_sig, function(cellstate){
  pt_filt <- meta_combi |> 
    filter(celltype_r2 == cellstate) |> 
    select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
    distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
    pivot_wider(values_from = count_r2, names_from = time_point) |> 
    filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5)) |> pull(patient)
  sub_long <- freq_long |> 
    filter(celltype_r2 == cellstate, patient %in% pt_filt, response != 'NE') |> 
    mutate(time_point2 = case_when(time_point == "Pre" ~ "Pre",
                                   time_point == "Post" & interval < 21 ~ "<3w",
                                   time_point == "Post" & interval >= 21 ~ ">=3w"),
           time_point2 = fct_relevel(time_point2, "Pre", "<3w", ">=3w"),
           pct_r2_comp_log1p = log1p(freq_r2_comp*100))
  p <- sub_long |> ggplot(aes(x = time_point2, y = pct_r2_comp_log1p)) +
    geom_boxplot(color = 'black', outlier.shape = NA) +
    geom_line(aes(group = patient), color = "gray", linetype = 'dotted') +
    geom_point(aes(color = cancertype, shape = modality), size=2, alpha = 0.6) +
    scale_color_manual(values = palette_cancertype) + 
    # facet_wrap(.~ response, scales = 'free') + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, size = 12, color = "black")) +
    xlab("") + ylab("Log(pct+1)") + ggtitle(cellstate) +
    labs(color = "Cancer Type", shape = "Modality") 
  stat.test <- rbind(sub_long |> filter(interval < 21) |> t_test(pct_r2_comp_log1p ~ time_point2, ref.group = '<3w') |> add_xy_position(x = "time_point2", dodge = 0.8), 
                     sub_long |> filter(interval > 21) |> t_test(pct_r2_comp_log1p ~ time_point2, ref.group = '>=3w') |> add_xy_position(x = "time_point2", dodge = 0.8))  
  p <- p + stat_pvalue_manual(stat.test) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  return(p)
})
wrap_plots(list_boxplot, ncol = 4) + plot_layout(guides = 'collect')

# One-by-one
cellstate <- 'CD8_Tm_IL7R'
pt_filt <- meta_combi |> 
  filter(celltype_r2 == cellstate) |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = count_r2, names_from = time_point) |> 
  filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5)) |> pull(patient)
sub_long <- freq_long |> 
  filter(celltype_r2 == cellstate, patient %in% pt_filt, response != 'NE') |> 
  mutate(time_point2 = case_when(time_point == "Pre" ~ "Pre",
                                 time_point == "Post" & interval < 21 ~ "<3w",
                                 time_point == "Post" & interval >= 21 ~ ">=3w"),
         time_point2 = fct_relevel(time_point2, "Pre", "<3w", ">=3w"),
         pct_r2_comp_log1p = log1p(freq_r2_comp*100))
p <- sub_long |> ggplot(aes(x = time_point2, y = pct_r2_comp_log1p)) +
  geom_boxplot(color = 'black', outlier.shape = NA) +
  geom_line(aes(group = patient), color = "gray", linetype = 'dotted') +
  geom_point(aes(color = cancertype, shape = modality), size=2, alpha = 0.6) +
  scale_color_manual(values = palette_cancertype) + 
  # facet_wrap(.~ response, scales = 'free') + 
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 0, size = 12, color = "black")) +
  xlab("") + ylab("Log(pct+1)") + ggtitle(cellstate) +
  labs(color = "Cancer Type", shape = "Modality") 

stat.test <- rbind(sub_long |> filter(interval < 21) |> t_test(pct_r2_comp_log1p ~ time_point2, ref.group = '<3w') |> add_xy_position(x = "time_point2", dodge = 0.8), 
      sub_long |> filter(interval > 21) |> t_test(pct_r2_comp_log1p ~ time_point2, ref.group = '>=3w') |> add_xy_position(x = "time_point2", dodge = 0.8))  
p + stat_pvalue_manual(stat.test)+scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# by response
list_main_res <- lapply(unique(freq_long$celltype_cate), function(major){
  print(major)
  cellstates <- freq_long |> filter(celltype_cate == major) |> pull(celltype_r2) |> unique()
  list_res <- lapply(c('RE', 'NR'), function(res){
    res_comp <- list()
    for (i in 1:length(cellstates)){
      print(cellstates[[i]])
      pt_filt <- meta_combi |> 
        filter(celltype_r2 == cellstates[[i]], response == res) |> 
        select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
        distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
        pivot_wider(values_from = count_r2, names_from = time_point) |> 
        filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5)) |> pull(patient)
      if (length(pt_filt) <= 3){
        res_comp[[i]] <- c(cellstates[[i]], NA, NA)
      } else {
        effectsize <- freq_long |> 
          filter(celltype_r2 == cellstates[[i]], patient %in% pt_filt) |> 
          cohens_d(freq_r2_comp ~ time_point, ref.group = 'Post', paired = TRUE)
        ttest <- freq_long |> 
          filter(celltype_r2 == cellstates[[i]], patient %in% pt_filt) |> 
          t_test(freq_r2_comp ~ time_point, ref.group = 'Post', paired = TRUE)
        fdr <- 
          res_comp[[i]] <- c(cellstates[[i]], effectsize$effsize, ttest$p)
      }
    }
    res_comp <- do.call(rbind, res_comp) |> data.frame()
    colnames(res_comp) <- c("Cell States", "Cohen's d", "Pvalue")
    res_comp$Response <- res
    res_comp <- res_comp |> 
      adjust_pvalue(p.col = 'Pvalue', method = 'fdr') |> 
      dplyr::filter(Pvalue < 0.05)
    return(res_comp)
  })
  df_res <- do.call(rbind, list_res)
  return(df_res)
})
names(list_main_res) <- unique(freq_long$celltype_cate)






