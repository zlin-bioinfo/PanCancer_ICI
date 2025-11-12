pkgs <- c('qs2','tidyr','dplyr','plyr','stringr','tibble','janitor','ggplot2','RColorBrewer','ComplexHeatmap','MetBrewer','survminer','survival','Seurat','CellChat','clusterProfiler','enrichplot','ggrepel','GSVA','forestplot','Seurat','Startrac','colorRamp2','rstatix','org.Hs.eg.db')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1) 
setwd("/home/zlin/PanCancer_ICI")
sample_info <- read.csv('tables/sample_info_time_updated.csv') 
source('scripts/Celltype_classification.R')
rownames(sample_info) <- NULL

deg.list <- qs_read('data/bulk_datasets/deg_list_inflamed.qs2')
deg.list$I2M <- rev(deg.list$I2M)
genelist <- unlist(deg.list)
seu$response <- sample_info$response[match(seu$sample, str_replace_all(sample_info$sample, '_', '-'))]
seu <- seu |> subset(subset = sample_group %in% c('i2bq','i2m'))
scores <- score(GetAssayData(seu, layer = 'data'), deg.list)
scores <- data.frame(scores, check.names = F)
seu$score <- scores$`I2B/Q`-scores$I2M
expr_mat_sub <- FetchData(seu, vars = genelist)
metadata <- seu@meta.data[names(sort(seu$score)),]
metadata$sample_group <- ifelse(metadata$sample_group == 'i2bq', 'I->B/Q', 'I->M')
metadata$subtype <- str_split(metadata$sample, '-', simplify = T)[,1]
col_fun = colorRamp2(c(-1, 1), c("white", "#CC0C00FF"))
col_ha = HeatmapAnnotation(
  Score = metadata$score,
  Response = metadata$response,
  `Group` = metadata$sample_group,
  col = list(`Group` = c('I->B/Q'="#009E73", 'I->M'="#E69F00"), 
             Response = c('R'='#CC0C00FF','NR'='#5C88DAFF'),
             Score = col_fun),
  annotation_legend_param = list(
    Score = list(
      title_position = "topcenter",
      at = c(-1, 0, 1),
      labels = c("Min", "", "Max")
    )
  ))
# Genes to highlight
genes_to_mark <- c("PRF1", "GNLY", "GZMH", 'GZMM', "NKG7",'KLRD1', "PDCD1", "LAG3", "ZAP70", 'CD8A',"CD247",'HLA-A',
                   "LYZ", 'MSR1', 'CD163L1', 'OLR1','GPNMB', 'PTGS2', 'PTAFR',"FOLR2",'FCGR2A',
                   'CCL3L1', 'CCL3', 'CCL2', 'CXCL8', 'CXCL11' 
)
ha = rowAnnotation(foo = anno_mark(at = which(genelist %in% genes_to_mark), 
                                   labels = genes_to_mark,
                                   labels_gp = gpar(fontsize = 7)))
pdf('figures/Dynamics/Inflamed/ht_expr.pdf', width = 5, height = 6)
Heatmap(expr_mat_sub[names(sort(seu$score)), genelist] |> scale() |> t(), 
        cluster_rows = F, cluster_columns = F, border = T,
        show_column_names = F,
        show_row_names = F,
        top_annotation = col_ha, 
        right_annotation = ha,
        col = circlize::colorRamp2(c(-3, 0, 3), c("#154999", "white", "#CF0034")),
        heatmap_legend_param = list(
          title = 'Expression',
          title_position = "topcenter",
          at = c(-3, 0, 3), labels = c("Min", "", "Max"), border = F),
)
dev.off()

metadata <- read.csv('tables/meta_all.csv')
paired_pt_info <- sample_info |>
  dplyr::filter(paired == 'Yes') |>
  dplyr::select(patient, tx_status, group, response, subtype, interval, res_metric, cohort, neoadjuvant, subtype) |> 
  pivot_wider(values_from = 'group', names_from = 'tx_status') |> 
  mutate(dynamics = ifelse(Baseline == Treated, 'Stable', 'Shifted'))
pt_i2q <- paired_pt_info |> 
  filter(Baseline == "Immune Inflamed", Treated == "Immune Quiescent") |> 
  pull(patient)
pt_i2b <- paired_pt_info |> 
  filter(Baseline == "Immune Inflamed", Treated == "B Cell-enriched") |> 
  pull(patient)
pt_i2m <- paired_pt_info |> 
  filter(Baseline == "Immune Inflamed", Treated == "Myeloid-enriched") |> 
  pull(patient)
pt_i2i <- paired_pt_info |>
  filter(Baseline == "Immune Inflamed", Treated == "Immune Inflamed", response == 'R') |>
  pull(patient)

# Boxplot 
seu <- qs_read('data/inflamed_seu_TME_bt.qs2')
gene_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = rownames(seu),
                                          columns = c("ENSEMBL", "GENETYPE", "GENENAME"), # Include other useful columns
                                          keytype = "SYMBOL")
protein_coding_rna <- gene_annotations |> filter(GENETYPE %in% c('protein-coding')) |> pull(SYMBOL)
seu <- seu[protein_coding_rna,]
seu <- seu |> 
  subset(subset = celltype_r2 %in% immune) |> 
  AverageExpression(group.by = 'sample', layer = 'data', return.seurat = T)
seu$sample_group <- 'i2bq'
seu$sample_group[seu$sample %in% str_replace_all(sample_i2i,'_','-')] <- 'i2i'
seu$sample_group[seu$sample %in% str_replace_all(sample_i2m,'_','-')] <- 'i2m'
scores <- score(GetAssayData(seu, layer = 'data'), deg.list)
scores <- data.frame(scores, check.names = F)
seu$score <- scores$`I2B/Q`-scores$I2M
metadata <- seu@meta.data 
metadata <- metadata |> mutate(sample_group = case_when(sample_group == 'i2bq' ~ 'I->B/Q',
                                                        sample_group == 'i2i' ~ 'I->I',
                                                        sample_group == 'i2m' ~ 'I->M'))
stat.test <- metadata |> 
  wilcox_test(score ~ sample_group) |> 
  add_significance(p.col = 'p.adj') |> 
  add_xy_position() |> 
  mutate(p.label = case_when(p < 0.001 ~"***", p < 0.01 ~"**", p < 0.05 ~"*", TRUE ~ NA_character_)) %>% filter(!is.na(p.label))
metadata |>
  ggplot() +
  geom_boxplot(aes(x = sample_group, y = score, fill = sample_group), width = 0.5, outliers = FALSE) +
  scale_fill_manual(values = c("#009E73","#56B4E9","#D55E00"), name = '') +
  stat_pvalue_manual(stat.test, hide.ns = 'p', label = 'p.label', tip.length = 0.001, bracket.size = 0.5, label.size = 5, position = position_dodge(0.8), inherit.aes = FALSE) +
  theme_classic2() +
  ylab('Score') +
  xlab('') + # Set x-axis label to blank
  theme(axis.text.x = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave('figures/Dynamics/Inflamed/boxplot_score_group.pdf', height = 3, width = 2.5)

# Correlation
# loading all baseline Inflamed samples
samples_included <- sample_info |> filter(group == 'Immune Inflamed', tx_status == 'Baseline') |> pull(sample)
datasets <- sample_info |> filter(group == 'Immune Inflamed', tx_status == 'Baseline') |> pull(cohort) |> unique()
datasets[which(datasets == "BCC&SCC_Yost")] <- 'BCC_Yost'
list_seu <- lapply(datasets, function(dataset) {
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = sample %in% samples_included & celltype_r2 %in% immune) |>
    NormalizeData() 
  return(seu)
})
seu <- merge(x = list_seu[[1]], y=list_seu[2:length(list_seu)])
seu <- JoinLayers(seu)
gene_annotations <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = rownames(seu),
                                          columns = c("ENSEMBL", "GENETYPE", "GENENAME"), # Include other useful columns
                                          keytype = "SYMBOL")
protein_coding_rna <- gene_annotations |> filter(GENETYPE %in% c('protein-coding')) |> pull(SYMBOL) |> unique()
seu <- seu[protein_coding_rna,]
seu <- seu[!str_detect(rownames(seu), '^MT-|^RP'),]
# pseudobulk for immune
expr_mat <- AverageExpression(seu, group.by = 'sample', layer = 'data')
# scoring for all samples
scores <- score(expr_mat$RNA, deg.list)
scores <- scores |> 
  data.frame(check.names = F) |> 
  mutate(score = `I2B/Q` - I2M)
sample_info_sub <- sample_info |> filter(group == 'Immune Inflamed', tx_status == 'Baseline', subset == 'TME')
sample_info_sub$score <- scores$score[match(str_replace_all(sample_info_sub$sample, '_', '-'), rownames(scores))]
# immune
mtx_immune <- metadata |> 
  filter(sample %in% sample_info_sub$sample, 
         count_immune >= 50, non_malignant_count >= 100, 
         celltype_r2 %in% immune,
         # subset == 'TME',
         subset %in% c('TME', 'CD45+sorted'),
         cohort != 'SKCM_this study') |> 
  mutate(celltype_r2 = case_when(celltype_r2 %in% c('iCAF_MMP1', 'iCAF_IL6') ~ 'iCAF',
                                 .default = celltype_r2)) |> 
  dplyr::select(sample,  freq_r2_comp, celltype_r2) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = celltype_r2, values_fill = 0) |> 
  column_to_rownames(var = 'sample') 
mtx <- cbind(mtx_immune[sample_info_sub$sample,], sample_info_sub$score)
names(mtx)[63] <- 'Score'
cor.res <- cor(mtx)
testRes <- corrplot::cor.mtest(mtx, conf.level = 0.95)
cor.test.res <- testRes$p
diag(cor.test.res) <- NA
df.res <- data.frame(cor_eff = cor.res[,'Score'], pvalue = cor.test.res[,'Score']) 
df.res |> filter(pvalue < 0.05) |> dplyr::arrange(desc(cor_eff))

df.res <- df.res %>%
  mutate(
    color_group = case_when(
      pvalue < 0.05 & cor_eff > 0 ~ "Positive",
      pvalue < 0.05 & cor_eff < 0 ~ "Negative",
      TRUE ~ "Non-significant" # For p-values between 0.05 and 0.1
    ), 
    cell_type = rownames(df.res)
  )

# Reorder the cell_type factor by cor_eff for plotting
df.res$cell_type <- factor(df.res$cell_type, levels = df.res$cell_type[order(df.res$cor_eff)])
df.res <- df.res[-which(rownames(df.res)=='Score'),]
# Define the custom color palette
custom_colors <- c("Positive" = "#CC0C00FF", "Negative" = "#5C88DAFF", "Non-significant" = "gray")

# Create the bar plot
df.res |> filter(abs(cor_eff)> 0.15) |> 
  mutate(color_group = factor(color_group, levels = c("Positive", "Non-significant", "Negative"))) |> 
  ggplot(aes(x = cell_type, y = cor_eff, fill = color_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = custom_colors, name = "") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add a line at y=0
  labs(
    title = "",
    x = "",
    y = "Correlation Coefficient"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(angle = 0,size = 14),
    axis.text.y = element_text(size = 13, color = 'black',angle = 0, hjust = 1, vjust = 0.5),
    axis.text.x = element_text(size = 11, color = 'black')
  )
ggsave('figures/Dynamics/Inflamed//barplot_cor.pdf', height = 7, width = 6)

# Prognostic Analysis
# 
deg.list <- qs_read('data/bulk_datasets/deg_list_inflamed.qs2')
expr_list <- qs_read('data/bulk_datasets/expr_list.qs2')
score_list <- lapply(expr_list, score, gs_list = deg.list)
scores <- do.call(rbind, score_list) 
scores <- data.frame(scores, check.names = F)
scores$sample <- rownames(scores)
scores$score <- scores$`I2B/Q`-scores$I2M

subtype <- "Immune Inflamed"
clin_info <- read.csv('tables/sample_info_time_bulk.csv')
clin_info$treatment[clin_info$cohort == 'NSCLC_SU2C-MARK'] <- 'aPD(L)1'
clin_info$treatment[clin_info$cohort == 'Melanoma_Gide'] <- 'aPD1/aPD1+CTLA4'
# Response
clin_info_filter <- clin_info |> filter(group == subtype, response %in% c('R','NR'), tx_status == 'Baseline')
n_pt <- nrow(clin_info_filter)
clin_info_filter$score <- scores$score[match(clin_info_filter$sample, scores$sample)] |> scale()
clin_info_filter$response <- ifelse(clin_info_filter$response == 'R', 1, 0) 
res <- glm(response ~ score + cohort, data = clin_info_filter, family = binomial(link = "logit"))
s_res <- summary(res)
coef_table <- s_res$coefficients
OR <- coef_table["score", "Estimate"] |> exp() |> round(2)
conf_interval_odds <- confint(res)
lower_CI <- conf_interval_odds["score", "2.5 %"] |> exp() |> round(2)
upper_CI <- conf_interval_odds["score", "97.5 %"] |> exp() |> round(2)
p_value <- coef_table["score", "Pr(>|z|)"]
p_value <- ifelse(p_value<0.001, '<0.001', round(p_value, 3))
cohorts <- clin_info_filter |> tabyl(cohort) |> 
  filter(n > 20) |>
  pull(cohort)
cohorts <- cohorts[-which(cohorts == "HNSC_CLB-IHN")]
list_res <- lapply(cohorts, function(x){
  print(x)
  clin_info_filter <- clin_info_filter |> filter(cohort == x)
  n_pt <- nrow(clin_info_filter)
  res <- glm(response ~ score, data = clin_info_filter, family = binomial(link = "logit"))
  s_res <- summary(res); s_res
  coef_table <- s_res$coefficients
  p_value <- coef_table["score", "Pr(>|z|)"]
  OR <- coef_table["score", "Estimate"] |> exp()
  conf_interval_odds <- confint(res)
  lower_CI <- conf_interval_odds["score", "2.5 %"] |> exp()
  upper_CI <- conf_interval_odds["score", "97.5 %"] |> exp()
  # Create the res_df for the current cohort
  res_df <- data.frame(
    cancertype = unique(clin_info_filter$cancertype),
    cohort = x,
    n = n_pt,
    tx = unique(clin_info_filter$treatment),
    OR = round(OR,2),
    mean = round(OR,2),
    lower = lower_CI,
    upper = upper_CI,
    p_value = round(p_value,3)
  )
  return(res_df)
})
df <- do.call(rbind, list_res)
p <- df |> forestplot(labeltext = c(cohort, n, tx, OR, p_value),
                      clip = c(0.1, 6),
                      vertices = TRUE,
                      xlog = T) |> 
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_lines("black") |>
  fp_add_header(cohort =  "Cohort",
                n = "No.patients",
                tx = "Treatment",
                OR = "OR",
                p_value = 'P value') |> 
  # fp_add_header("Mye-Inflamed vs Inflamed stable") |>
  fp_append_row(mean  = OR,
                lower = lower_CI,
                upper = upper_CI,
                cohort = "Overall (adjusted)",
                OR = OR,
                n = n_pt,
                p_value = p_value,
                is.summary = TRUE) |>
  fp_set_zebra_style("#EFEFEF") |> 
  fp_decorate_graph(box = gpar(lty = 2, col = "lightgray"),
                    graph.pos = 4) ;p
pdf('figures/Dynamics/Inflamed/fp_res.pdf', height = 3, width = 10)
p
dev.off()

# PFS
clin_info_filter <- clin_info |> filter(group == subtype, !is.na(pfs), tx_status == 'Baseline')
n_pt <- nrow(clin_info_filter)
clin_info_filter$score <- scores$score[match(clin_info_filter$sample, scores$sample)] |> scale()
clin_info_filter$os <- ifelse(clin_info_filter$os %in% c(0, 'Alive'), 0, 1)
res <- coxph(Surv(time_pfs, pfs) ~ score + cohort, data = clin_info_filter)
s_res <- summary(res)
HR <- s_res$coefficients["score", "exp(coef)"] |> round(2)
lower_CI <- s_res$conf.int["score", "lower .95"] |> round(2)
upper_CI <- s_res$conf.int["score", "upper .95"] |> round(2)
p_value <- s_res$coefficients["score", "Pr(>|z|)"] |> round(3)
cohorts <- clin_info_filter |> tabyl(cohort) |> 
  filter(n>20) |>
  pull(cohort)
list_res <- lapply(cohorts, function(x){
  print(x)
  df_filter <- clin_info_filter |> filter(cohort == x)
  res <- coxph(Surv(time_pfs, pfs) ~ score, data = df_filter)
  s_res <- summary(res)
  HR <- s_res$coefficients["score", "exp(coef)"]
  lower_CI <- s_res$conf.int["score", "lower .95"]
  upper_CI <- s_res$conf.int["score", "upper .95"]
  p_value <- s_res$coefficients["score", "Pr(>|z|)"]
  res_df <- data.frame(
    cancertype = unique(df_filter$cancertype),
    cohort = x,
    n = nrow(df_filter),
    tx = unique(df_filter$treatment),
    HR = round(HR,2),
    mean = round(HR,2),
    lower = lower_CI,
    upper = upper_CI,
    p_value = round(p_value,3)
  )
  return(res_df)
})
df <- do.call(rbind, list_res)
p <- df |> forestplot(labeltext = c(cohort, n, tx, HR, p_value),
                      clip = c(0.001, 2),
                      vertices = TRUE,
                      xlog = TRUE) |> 
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_lines("steelblue") |>
  fp_add_header(cohort = "Cohort",
                n = "No.patients",
                tx = "Treatment",
                HR = "HR(PFS)",
                p_value = 'P value') |> 
  # fp_add_header("") |>
  fp_append_row(mean  = HR,
                lower = lower_CI,
                upper = upper_CI,
                cohort = "Overall (adjusted)",
                HR = HR,
                n = n_pt,
                p_value = p_value,
                is.summary = TRUE) |>
  fp_set_zebra_style("#EFEFEF") |> 
  fp_decorate_graph(box = gpar(lty = 2, col = "lightgray"),
                    graph.pos = 4) ;p
pdf('figures/Dynamics/Inflamed/fp_pfs.pdf', height = 2.5, width = 10)
p
dev.off()
# OS
clin_info_filter <- clin_info |> filter(group == subtype, !is.na(os), tx_status == 'Baseline')
n_pt <- nrow(clin_info_filter)
clin_info_filter$score <- scores$score[match(clin_info_filter$sample, scores$sample)] |> scale()
clin_info_filter$os <- ifelse(clin_info_filter$os %in% c(0, 'Alive'), 0, 1)
res <- coxph(Surv(time_os, os) ~ score + cohort, data = clin_info_filter)
s_res <- summary(res)
s_res <- summary(res)
HR <- s_res$coefficients["score", "exp(coef)"] |> round(2)
lower_CI <- s_res$conf.int["score", "lower .95"] |> round(2)
upper_CI <- s_res$conf.int["score", "upper .95"] |> round(2)
p_value <- s_res$coefficients["score", "Pr(>|z|)"] |> round(3)
cohorts <- clin_info_filter |> tabyl(cohort) |> 
  filter(n>20) |>
  pull(cohort)
list_res <- lapply(cohorts, function(x){
  print(x)
  df_filter <- clin_info_filter |> filter(cohort == x)
  # df_filter$score <- scores$score[match(df_filter$sample, scores$sample)] |> scale()
  res <- coxph(Surv(time_os, os) ~ score, data = df_filter)
  s_res <- summary(res)
  HR <- s_res$coefficients["score", "exp(coef)"]
  lower_CI <- s_res$conf.int["score", "lower .95"]
  upper_CI <- s_res$conf.int["score", "upper .95"]
  p_value <- s_res$coefficients["score", "Pr(>|z|)"]
  res_df <- data.frame(
    cancertype = unique(df_filter$cancertype),
    cohort = x,
    n = nrow(df_filter),
    tx = unique(df_filter$treatment),
    HR = round(HR,2),
    mean = round(HR,2),
    lower = lower_CI,
    upper = upper_CI,
    p_value = round(p_value,3)
  )
  return(res_df)
})
df <- do.call(rbind, list_res)
p <- df |> forestplot(labeltext = c(cohort, n, tx, HR, p_value),
                      clip = c(0.001, 2),
                      vertices = T,
                      xlog = TRUE) |> 
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_lines("black") |>
  fp_add_header(cohort =  "Cohort",
                n = "No.patients",
                tx = "Treatment",
                HR = "HR(OS)",
                p_value = 'P value') |> 
  # fp_add_header("") |>
  fp_append_row(mean  = HR,
                lower = lower_CI,
                upper = upper_CI,
                cohort = "Overall (adjusted)",
                HR = HR,
                n = n_pt,
                p_value = p_value,
                is.summary = TRUE) |>
  fp_set_zebra_style("#EFEFEF") |> 
  fp_decorate_graph(box = gpar(lty = 2, col = "lightgray"),
                    graph.pos = 4) ;p
pdf('figures/Dynamics/Inflamed/fp_os.pdf', height = 3, width = 10)
p
dev.off()

# K-M plot
subtype <- "Immune Inflamed"
clin_info <- read.csv('tables/sample_info_time_bulk.csv')
clin_info_filter <- clin_info |> filter(group == subtype, !is.na(os), tx_status == 'Baseline')
clin_info_filter$score <- scores$score[match(clin_info_filter$sample, scores$sample)] 
clin_info_filter$os <- ifelse(clin_info_filter$os %in% c(0, 'Alive'), 0, 1)
res.cut <- surv_cutpoint(clin_info_filter, time = "time_os", event = "os", variables = 'score', minprop = 0.35)
res.cat <- surv_categorize(res.cut)
clin_info_filter$group_bt <- res.cat[, 'score']
fit <- survfit(Surv(time_os, os) ~ group_bt, data = clin_info_filter)
p <- ggsurvplot(fit, data = clin_info_filter,conf.int = F, 
                pval = T, 
                pval.coord = c(40, 0.75),
                legend.title = '',
                legend.labs = c('High','Low'),
                # surv.median.line = 'hv',
                risk.table = F,        
                risk.table.col = "strata",
                risk.table.height = 0.25,
                ggtheme = theme_classic(),
                palette = c('High' = '#CC0C00FF','Low' = '#5C88DAFF'),
                xlab = 'Time(OS)', title=''); p
pdf('figures/Dynamics/Inflamed/KM_os.pdf', height = 4, width = 4)
print(p, newpage = FALSE)
dev.off()

clin_info_filter <- clin_info |> filter(group == subtype, !is.na(pfs), tx_status == 'Baseline')
clin_info_filter$score <- scores$score[match(clin_info_filter$sample, scores$sample)] 
clin_info_filter$pfs <- ifelse(clin_info_filter$pfs %in% c(0, 'Alive'), 0, 1)
res.cut <- surv_cutpoint(clin_info_filter, time = "time_pfs", event = "pfs", variables = 'score', minprop = 0.35)
res.cat <- surv_categorize(res.cut)
clin_info_filter$group_bt <- res.cat[, 'score']
fit <- survfit(Surv(time_pfs, pfs) ~ group_bt, data = clin_info_filter)
p <- ggsurvplot(fit, data = clin_info_filter,conf.int = F, 
                pval = T, 
                pval.coord = c(30, 0.75),
                legend.title = '',
                legend.labs = c('High','Low'),
                risk.table = F,        
                risk.table.col = "strata",
                risk.table.height = 0.25,
                ggtheme = theme_classic(),
                palette = c('High' = '#CC0C00FF','Low' = '#5C88DAFF'),
                xlab = 'Time(PFS)', title=''); p
pdf('figures/Dynamics/Inflamed/KM_pfs.pdf', height = 4, width = 4)
print(p, newpage = FALSE)
dev.off()







