rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','gtools','ComplexHeatmap','dittoSeq','RColorBrewer','ggpubr','tibble','epitools','dior','cowplot','ggnewscale','export','forcats','qs','purrr','effsize','corrplot','lmerTest','rstatix')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Frequency by composition
meta_combi <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_int.csv') |> 
  select(!X) |> 
  filter(proliferating == 'yes') |> 
  distinct(sample, celltype_r2, .keep_all = T)

freq_wide <- meta_combi |> 
  filter(count_r2 >= 10) |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_prol, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_prol, names_from = time_point) |> 
  replace_na(list(Post = 0, Pre = 0)) |> 
  mutate(change = log2((Post + 0.01)/(Pre + 0.01)), diff = (Post - Pre))
pal <- colorRampPalette(brewer.pal(10, "RdBu"))
# T/NK cell
celltype <- c('CD4_Naive','CD4_Tn_ADSL','CD4_Tm_TNF','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM',
              'CD4_Tm_CCL5', 'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
              'CD4_Treg_TNFRSF9-', 'CD4_Treg_S1PR1', 'CD4_Treg_TNFRSF9', 'CD4_Treg_ISG', 'CD4_Th_ISG', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_STR',
              'CD8_Naive', 'CD8_MAIT_SLC4A10', 'CD8_Tm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
              'CD8_Tex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13', 'CD8_Tex_OXPHOS-', 'CD8_STR', 'CD8_ISG', 'CD8_Temra_CX3CR1', 'CD8_NK-like')
mat_change <- freq_wide |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014")) |> 
  select(patient, diff, celltype_r2) |> 
  pivot_wider(values_from = diff, names_from = patient) |> 
  column_to_rownames(var = 'celltype_r2')
mat_change[is.na(mat_change)] <- 0
M <- cor(t(mat_change))
testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
corrplot(M, p.mat = testRes$p, method = 'color', col = rev(pal(100)), 
         tl.cex = 0.6, pch.cex = 1, tl.col = 'black', order = 'hclust', mar = c(0,1,1,0),
         sig.level = 0.05, insig = 'label_sig', title = 'TME')

uni_fe <- function(meta_combi, n.min = 10, pt.min = 15){
  df_check <- meta_combi |>
    select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
    distinct(celltype_r2, patient, time_point, .keep_all = T) |>
    pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
    filter(Pre >= n.min | Post >= n.min)
  subtypes <- 
    df_check |> 
    group_by(celltype_r2) |>
    summarise(count = n()) |> 
    filter(count >= pt.min) |>
    pull(celltype_r2)
  freq_mat <- distinct(meta_combi, sample, celltype_r2, .keep_all = T)
  uni_models <- lapply(subtypes, function(subtype){
    print(subtype)
    freq_mat$timepoint <- ifelse(freq_mat$time_point == 'Pre', 0, 1)
    formula <- as.formula(paste("timepoint ~ freq_r2_prol + dataset + response + (1 | patient)"))
    model <- lmer(formula, data = freq_mat[freq_mat$celltype_r2 == subtype,], REML = FALSE)
    model_summary <- summary(model)
    coef_table <- coef(summary(model))
    confint_table <- confint(model, level = 0.95)
    subtype_results <- coef_table['freq_r2_prol', , drop = FALSE]
    subtype_confint <- confint_table['freq_r2_prol', , drop = FALSE]
    combined_results <- data.frame(
      Celltypes = subtype,
      Estimate = subtype_results[1],
      StdError = subtype_results[2],
      tValue = subtype_results[4],
      pValue = coef_table['freq_r2_prol', 'Pr(>|t|)'],
      CI_lower = subtype_confint[1],
      CI_upper = subtype_confint[2]
    )
    return(combined_results)
  })
  df <- do.call(rbind, uni_models) |> data.frame()
  df$fdr <- p.adjust(df$pValue, method = 'fdr', n = nrow(df))
  # celltype_keep <- filter(df, pValue < 0.05) |> pull(Celltypes) |> unique()
  # df <- filter(df, Celltypes %in% celltype_keep)
  return(df)
}
# Overall
res_fe <- uni_fe(meta_combi); 
res_fe <- res_fe |> arrange(desc(Estimate))
res_fe |> filter(pValue < 0.05)

meta_combi$int_cat <- ifelse(meta_combi$interval < 21, 'Short','Long')
uni_logi <- function(meta_combi, timepoint = 'Pre', n.min = 10, pt.min = 20, vars_to_include = c("freq_r2_prol", "int_cat", "datset", "modality")){
  meta_df <- meta_combi |>
    filter(!response == 'NE', time_point == timepoint, celltype_r2 >= n.min) 
  subtypes <- meta_df |>
    group_by(celltype_r2, response) |>
    mutate(n_res = n()) |>
    filter(response == 'RE' & n_res >= pt.min | response == 'NR' & n_res >= pt.min) |>
    distinct(celltype_r2, .keep_all = T) |>
    select(celltype_r2, response) |>
    group_by(celltype_r2) |>
    summarize(n = n()) |>
    filter(n == 2) |>
    pull(celltype_r2)
  uni_models <- lapply(subtypes, function(subtype){
    print(subtype)
    df <- filter(meta_df, celltype_r2 == subtype)
    df$response <- ifelse(df$response == 'RE', 1, 0)
    vars_to_include <- vars_to_include[sapply(vars_to_include, function(x) length(unique(df[[x]])) > 1)]
    formula <- paste("response ~", paste(vars_to_include, collapse = " + "))
    model <- glm(formula, data = df, family = binomial)
    coef_table <- coef(summary(model))
    confint_table <- confint(model, level = 0.95)
    subtype_results <- coef_table['freq_r2_prol', , drop = FALSE]
    subtype_confint <- confint_table['freq_r2_prol', , drop = FALSE]
    combined_results <- data.frame(
      Celltypes = subtype,
      Estimate = subtype_results[1],
      StdError = subtype_results[2],
      pValue = coef_table['freq_r2_prol', 'Pr(>|z|)'],
      CI_lower = subtype_confint[1],
      CI_upper = subtype_confint[2]
    )
    return(combined_results)
  })
  results <- do.call(rbind, uni_models) |> data.frame()
  results$fdr <- p.adjust(results$pValue, method = 'fdr', n = nrow(results))
  # celltype_keep <- filter(results, fdr < 0.05) |> pull(Celltypes)
  # results <- filter(results, Celltypes %in% celltype_keep)
  return(results)
}
res_logi_pre <- uni_logi(meta_combi, pt.min = 15, timepoint = 'Pre'); res_logi_pre
res_logi_post <- uni_logi(meta_combi, pt.min = 15, timepoint = 'Post'); res_logi_post

res_logi_pre |> filter(pValue < 0.05)
res_logi_post |> filter(pValue < 0.05)




