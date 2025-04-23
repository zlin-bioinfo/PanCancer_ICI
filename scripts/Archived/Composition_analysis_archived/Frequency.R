rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','gtools','RColorBrewer','tibble','forcats','qs','effsize','corrplot','lmerTest','rstatix','MetBrewer')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# Frequency by composition
meta_combi <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_int.csv') |> select(!X)
head(meta_combi)
range(meta_combi$freq_r2_comp)
meta_combi$int_cat <- ifelse(meta_combi$interval < 21, '< 21d', '>= 21d')
freq_filt <- meta_combi |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, dataset, component, count_r2) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = count_r2, names_from = time_point)
freq_filt$pt_r2 <- paste0(freq_filt$patient, '_', freq_filt$celltype_r2)
freq_wide <- meta_combi |> 
  select(patient, time_point, celltype_r2, interval, cancertype, response, res_metric, treatment, prior, modality, freq_r2_comp, dataset, component) |> 
  distinct(patient, time_point, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(change = log2((Post + 0.01)/(Pre + 0.01)), diff = (Post - Pre))
freq_wide$pt_r2 <- paste0(freq_wide$patient, '_', freq_wide$celltype_r2)
freq_wide <- filter(freq_wide, pt_r2 %in% freq_filt$pt_r2)

# Correlation (cell states)
# Immune cells
celltype <- c('CD4_Naive','CD4_Tm_TNF','CD4_Tm_AREG','CD4_Tm_TIMP1','CD4_Tm_CAPG','CD4_Tm_CREM',
              'CD4_Tm_CCL5', 'CD4_Tem_GZMK', 'CD4_Temra_CX3CR1', 'CD4_pre-Tfh_CXCR5','CD4_Tfh_CXCR5','CD4_TfhTh1_IFNG', 
              'CD4_Treg_Early', 'CD4_Treg_ISG', 'CD4_Treg_TNFRSF9', 'CD4_Th_ISG', 'CD4_Th17_IL26','CD4_Th17_CCR6','CD4_Prolif',
              'CD8_Prolif', 'CD8_Naive', 'CD8_Tcm_IL7R', 'CD8_Trm_ZNF683', 'CD8_Tem_Early', 'CD8_Tem_GZMK', 
              'CD8_Tpex_TCF7', 'CD8_Tex_GZMK', 'CD8_Tex_CXCL13', 'CD8_Tex_OXPHOS-', 
              'CD8_ISG', 'CD8_Temra_CX3CR1', 'CD8_NK-like', 'MAIT', 'gdT', 'NK_CD56loCD16hi', 'NK_CD56hiCD16lo',
              'Naive B cells', 'Non-switched memory B cells', 'Switched memory B cells', 'Exhausted B cells', 'Plasma',
              'Mast','pDC','cDC1', 'cDC2_CD1C', 'cDC2_IL1B', 'cDC2_ISG15', 'cDC2_CXCL9', 'DC_LC-like', 'cCD2_MoDC', 'DC_Migr',
              'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
              'Macro_NLRP3','Macro_ISG15', 'Macro_INHBA', 'Macro_FN1', 'Macro_SPP1',
              'Macro_LYVE1','Macro_IL1B','Macro_C1QC','Macro_TREM2')
mat_change <- freq_wide |> 
  filter(dataset != 'NSCLC_Liu', 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014"), 
         celltype_r2 %in% celltype) |> 
  select(patient, diff, celltype_r2) |> 
  pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
  column_to_rownames(var = 'celltype_r2')
M <- cor(t(mat_change))
testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
pdf('/bigdata/zlin/Melanoma_meta/figures/Co-regulating/freq_immune.pdf')
pal <- colorRampPalette(met.brewer("Benedictus"))
corrplot(M, p.mat = testRes$p, method = 'color', col = rev(pal(100)), 
         tl.cex = 0.5, pch.cex = 0.7, tl.col = 'black', order = 'hclust', mar = c(0,1,1,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = 'Immune')
dev.off()
# TME
mat_change <- freq_wide |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014")) |> 
  select(patient, diff, celltype_r2) |> 
  pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
  column_to_rownames(var = 'celltype_r2') 
corrplot(M, order = 'hclust', addrect = 2)
corrplot(M, method = 'color', col = rev(pal(100)), 
         tl.cex = 0.4, pch.cex = 0.4, tl.col = 'black', order = 'hclust', mar = c(0,1,1,0), diag = FALSE, title = '', addrect =7)
M <- cor(t(mat_change))
testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
pdf('/bigdata/zlin/Melanoma_meta/figures/Co-regulating/freq_TME.pdf')
corrplot(M, p.mat = testRes$p, method = 'color', col = rev(pal(100)), 
         tl.cex = 0.4, pch.cex = 0.4, tl.col = 'black', order = 'AOE', mar = c(0,1,1,0),
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = 'TME')
dev.off()

mat_change <- freq_wide |> 
  filter(dataset %in% c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'BCC/SCC_Yost', 'BCC/SCC_Yost', 'PCa_Hawley'), 
         !patient %in% c("BCC/SCC_Yost_su009", "BCC/SCC_Yost_su011", "BCC/SCC_Yost_su012", "BCC/SCC_Yost_su014"),
         response == 'NR') |> 
  select(patient, diff, celltype_r2) |> 
  pivot_wider(values_from = diff, names_from = patient, values_fill = 0) |> 
  column_to_rownames(var = 'celltype_r2') 
M <- cor(t(mat_change))
testRes <- cor.mtest(t(mat_change), conf.level = 0.95)
pdf('/bigdata/zlin/Melanoma_meta/figures/Co-regulating/freq_NR_TME.pdf')
corrplot(M, p.mat = testRes$p, method = 'color', col = rev(brewer.pal(8, 'RdBu')), 
         tl.cex = 0.4, pch.cex = 0.6, tl.col = 'black', order = 'hclust', mar = c(0,1,1,0), 
         sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', title = 'NR_TME', diag = F)
dev.off()


# Mixed-Effects Models
uni_lmer <- function(meta_combi, n.sample = 10){
  df_check <- meta_combi |>
    select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
    distinct(celltype_r2, patient, time_point, .keep_all = T) |>
    pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
    filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5))
  subtypes <- 
    df_check |> 
    group_by(celltype_r2) |>
    summarise(count = n()) |> 
    filter(count >= n.sample) |>
    pull(celltype_r2)
  freq_mat <- distinct(meta_combi, sample, celltype_r2, .keep_all = T)
  uni_models <- lapply(subtypes, function(subtype){
    print(subtype)
    freq_mat$timepoint <- ifelse(freq_mat$time_point == 'Pre', 0, 1)
    formula <- as.formula("timepoint ~ freq_r2_comp + dataset + response + modality + int_cat + (1 | patient)")
    model <- lmer(formula, data = freq_mat[freq_mat$celltype_r2 == subtype,], REML = FALSE)
    model_summary <- summary(model)
    coef_table <- coef(summary(model))
    confint_table <- confint(model, level = 0.95)
    subtype_results <- coef_table['freq_r2_comp', , drop = FALSE]
    subtype_confint <- confint_table['freq_r2_comp', , drop = FALSE]
    combined_results <- data.frame(
      Celltypes = subtype,
      Estimate = subtype_results[1],
      StdError = subtype_results[2],
      tValue = subtype_results[4],
      pValue = coef_table['freq_r2_comp', 'Pr(>|t|)'],
      CI_lower = subtype_confint[1],
      CI_upper = subtype_confint[2]
    )
    return(combined_results)
  })
  df <- do.call(rbind, uni_models) |> data.frame()
  df$fdr <- p.adjust(df$pValue, method = 'fdr', n = nrow(df)) 
  celltype_keep <- filter(df, pValue < 0.05) |> pull(Celltypes) |> unique()
  df <- filter(df, Celltypes %in% celltype_keep)
  return(df)
}
# Overall
res_lmer <- uni_lmer(meta_combi); res_lmer
res_lmer <- res_lmer |> arrange(desc(Estimate))
res_lmer$fdr_cat <- '< 0.1'
res_lmer$fdr_cat[res_lmer$fdr>0.1] <- '> 0.1'
# Making plot
pdf('/bigdata/zlin/Melanoma_meta/figures/Change/uni_lmer_all.pdf', height = 5, width = 5)
p <- ggplot(res_lmer, aes(x= factor(Celltypes, levels = rev(res_lmer$Celltypes)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(aes(color = fdr_cat), size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.001), -log10(0.005), -log10(0.01)), labels = c('<0.001', '<0.005', '<0.01'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-10, 13)) + xlab("") +
  scale_color_manual(values = c("black", "grey")) +
  coord_flip() + theme_minimal() +
  labs(size = "P value", color = "FDR") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.35, "npc"),
  x1 = unit(0.25, "npc"),
  y = unit(0.17, "npc"),
  y1 = unit(0.17, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.7, "npc"),
  x1 = unit(0.8, "npc"),
  y = unit(0.17, "npc"),
  y1 = unit(0.17, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.3, "npc"), y = unit(0.15, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.75, "npc"), y = unit(0.15, "npc"), gp = gpar(fontsize = 8))
dev.off()

# specific cell states
subtype <- 'CD8_Tem_GZMK'
res_lmer |> filter(Celltypes == subtype)
freq_mat <- meta_combi |> 
  filter(celltype_r2 == subtype) |>
  distinct(sample, celltype_r2, .keep_all = T) 

uni_ds <- lapply(setdiff(unique(freq_mat$dataset), 'PCa_Hawley'), function(ds){
  mat <- freq_mat[freq_mat$dataset == ds,]
  print(ds)
  mat$timepoint <- ifelse(mat$time_point == 'Pre', 0, 1)
  variables_to_use <- c('freq_r2_comp', 'response'[length(unique(mat$response)) > 1])
  variables_to_use <- c(variables_to_use, "(1 | patient)")
  formula <- as.formula(paste("timepoint ~", paste((variables_to_use), collapse = " + ")))
  model <- lmer(formula, data = mat, REML = FALSE)
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  subtype_results <- coef_table['freq_r2_comp', , drop = FALSE]
  subtype_confint <- confint_table['freq_r2_comp', , drop = FALSE]
  combined_results <- data.frame(
    Celltypes = subtype,
    Estimate = subtype_results[1],
    StdError = subtype_results[2],
    tValue = subtype_results[4],
    pValue = coef_table['freq_r2_comp', 'Pr(>|t|)'],
    CI_lower = subtype_confint[1],
    CI_upper = subtype_confint[2],
    dataset = ds
  )
  return(combined_results)
})

res_df <- do.call(rbind, uni_ds)
res_subtype_all <- res_lmer |> filter(Celltypes == subtype) |> select(!fdr_cat)
res_subtype_all$fdr <- 'Overall'
names(res_subtype_all)[[8]] <- 'dataset'
res_df <- rbind(res_subtype_all, res_df) 
res_df <- res_df |> replace_na(list(CI_lower = -20, CI_upper = 40))
res_df$CI_lower[res_df$CI_lower < -20] <- -20
res_df$CI_upper[res_df$CI_upper > 40] <- 40
order <- res_df$dataset
res_df <- res_df |> filter(!dataset == 'SCC_Yost')
pdf(paste0('/bigdata/zlin/Melanoma_meta/figures/Change/', subtype, '.pdf'), height = 4, width = 5)
p <- ggplot(res_df, aes(x= factor(dataset, levels = rev(order)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, color = '#A4DBFF', position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.05), -log10(0.1), -log10(0.5)), labels = c('<0.05', '<0.1', '<0.5'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-20, 40), breaks = c(-20, -10, 0, 10, 20, 30, 40), labels = c("<-20", "-10", "0", "10", "20", "30", ">40")) + 
  coord_flip() + theme_minimal() + ggtitle(unique(res_df$Celltypes)) +
  labs(size = "P value") + xlab("") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.35, "npc"),
  x1 = unit(0.25, "npc"),
  y = unit(0.2, "npc"),
  y1 = unit(0.2, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.7, "npc"),
  x1 = unit(0.8, "npc"),
  y = unit(0.2, "npc"),
  y1 = unit(0.2, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.3, "npc"), y = unit(0.18, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.75, "npc"), y = unit(0.18, "npc"), gp = gpar(fontsize = 8))
dev.off()

subtype <- 'CD8_Tcm_IL7R'
res_lmer |> filter(Celltypes == subtype)
freq_mat <- meta_combi |> 
  filter(celltype_r2 == subtype) |>
  distinct(sample, celltype_r2, .keep_all = T) 

uni_ds <- lapply(setdiff(unique(freq_mat$dataset), 'PCa_Hawley'), function(ds){
  mat <- freq_mat[freq_mat$dataset == ds,]
  print(ds)
  mat$timepoint <- ifelse(mat$time_point == 'Pre', 0, 1)
  variables_to_use <- c('freq_r2_comp', 'response'[length(unique(mat$response)) > 1])
  variables_to_use <- c(variables_to_use, "(1 | patient)")
  formula <- as.formula(paste("timepoint ~", paste((variables_to_use), collapse = " + ")))
  model <- lmer(formula, data = mat, REML = FALSE)
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  subtype_results <- coef_table['freq_r2_comp', , drop = FALSE]
  subtype_confint <- confint_table['freq_r2_comp', , drop = FALSE]
  combined_results <- data.frame(
    Celltypes = subtype,
    Estimate = subtype_results[1],
    StdError = subtype_results[2],
    tValue = subtype_results[4],
    pValue = coef_table['freq_r2_comp', 'Pr(>|t|)'],
    CI_lower = subtype_confint[1],
    CI_upper = subtype_confint[2],
    dataset = ds
  )
  return(combined_results)
})

res_df <- do.call(rbind, uni_ds)
res_subtype_all <- res_lmer |> filter(Celltypes == subtype) |> select(!fdr_cat)
res_subtype_all$fdr <- 'Overall'
names(res_subtype_all)[[8]] <- 'dataset'
res_df <- rbind(res_subtype_all, res_df) 
res_df <- res_df |> replace_na(list(CI_lower = -10, CI_upper = 60))
res_df$CI_lower[res_df$CI_lower < -10] <- -10
res_df$CI_upper[res_df$CI_upper > 60] <- 60
order <- res_df$dataset
pdf(paste0('/bigdata/zlin/Melanoma_meta/figures/Change/', subtype, '.pdf'), height = 4, width = 5)
p <- ggplot(res_df, aes(x= factor(dataset, levels = rev(order)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, color = '#A4DBFF', position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.05), -log10(0.1), -log10(0.5)), labels = c('<0.05', '<0.1', '<0.5'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-10, 60), breaks = c(-10, 0, 10, 20, 30, 40, 50, 60), labels = c("<-10", "0", "10", "20", "30", "40", "50", ">60")) + 
  coord_flip() + theme_minimal() + ggtitle(unique(res_df$Celltypes)) +
  labs(size = "P value") + xlab("") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.35, "npc"),
  x1 = unit(0.25, "npc"),
  y = unit(0.21, "npc"),
  y1 = unit(0.21, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.7, "npc"),
  x1 = unit(0.8, "npc"),
  y = unit(0.21, "npc"),
  y1 = unit(0.21, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.3, "npc"), y = unit(0.19, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.75, "npc"), y = unit(0.19, "npc"), gp = gpar(fontsize = 8))
dev.off()

# CD4_Th_ISG
subtype <- 'CD4_Th_ISG'
res_lmer |> filter(Celltypes == subtype)
freq_mat <- meta_combi |> 
  filter(celltype_r2 == subtype) |>
  distinct(sample, celltype_r2, .keep_all = T) 

uni_ds <- lapply(setdiff(unique(freq_mat$dataset), c('SKCM_Becker', 'PCa_Hawley', 'CRC_Li', 'SCC_Yost')), function(ds){
  mat <- freq_mat[freq_mat$dataset == ds,]
  print(ds)
  mat$timepoint <- ifelse(mat$time_point == 'Pre', 0, 1)
  variables_to_use <- c('freq_r2_comp', 'response'[length(unique(mat$response)) > 1])
  variables_to_use <- c(variables_to_use, "(1 | patient)")
  formula <- as.formula(paste("timepoint ~", paste((variables_to_use), collapse = " + ")))
  model <- lmer(formula, data = mat, REML = FALSE)
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  subtype_results <- coef_table['freq_r2_comp', , drop = FALSE]
  subtype_confint <- confint_table['freq_r2_comp', , drop = FALSE]
  combined_results <- data.frame(
    Celltypes = subtype,
    Estimate = subtype_results[1],
    StdError = subtype_results[2],
    tValue = subtype_results[4],
    pValue = coef_table['freq_r2_comp', 'Pr(>|t|)'],
    CI_lower = subtype_confint[1],
    CI_upper = subtype_confint[2],
    dataset = ds
  )
  return(combined_results)
})

res_df <- do.call(rbind, uni_ds)
res_subtype_all <- res_lmer |> filter(Celltypes == subtype) |> select(!fdr_cat)
res_subtype_all$fdr <- 'Overall'
names(res_subtype_all)[[8]] <- 'dataset'
res_df <- res_df |> filter(!dataset == 'BCC/SCC_Yost')
res_df <- rbind(res_subtype_all, res_df) 
res_df <- res_df |> replace_na(list(CI_lower = -30, CI_upper = 15))
res_df$CI_lower[res_df$CI_lower < -30] <- -30
res_df$CI_upper[res_df$CI_upper > 15] <- 15
order <- res_df$dataset
pdf(paste0('/bigdata/zlin/Melanoma_meta/figures/Change/', subtype, '.pdf'), height = 4, width = 5)
p <- ggplot(res_df, aes(x= factor(dataset, levels = rev(order)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, color = '#A4DBFF', position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.05), -log10(0.1), -log10(0.5)), labels = c('<0.05', '<0.1', '<0.5'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-30, 15), breaks = c(-30, -20,-10, 0, 10, 15), labels = c("<-30", "-20", "-10", "0", "10", ">15")) + 
  coord_flip() + theme_minimal() + ggtitle(unique(res_df$Celltypes)) +
  labs(size = "P value") + xlab("") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.35, "npc"),
  x1 = unit(0.25, "npc"),
  y = unit(0.2, "npc"),
  y1 = unit(0.2, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.7, "npc"),
  x1 = unit(0.8, "npc"),
  y = unit(0.2, "npc"),
  y1 = unit(0.2, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.3, "npc"), y = unit(0.18, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.75, "npc"), y = unit(0.18, "npc"), gp = gpar(fontsize = 8))
dev.off()

# cDC2_ISG15
subtype <- 'cDC2_ISG15'
res_lmer |> filter(Celltypes == subtype)
freq_mat <- meta_combi |> 
  filter(celltype_r2 == subtype) |>
  distinct(sample, celltype_r2, .keep_all = T) 

uni_ds <- lapply(setdiff(unique(freq_mat$dataset), c('PCa_Hawley', 'BCC/SCC_Yost','CRC_Li','TNBC_Zhang')), function(ds){
  mat <- freq_mat[freq_mat$dataset == ds,]
  print(ds)
  mat$timepoint <- ifelse(mat$time_point == 'Pre', 0, 1)
  variables_to_use <- c('freq_r2_comp', 'response'[length(unique(mat$response)) > 1])
  variables_to_use <- c(variables_to_use, "(1 | patient)")
  formula <- as.formula(paste("timepoint ~", paste((variables_to_use), collapse = " + ")))
  model <- lmer(formula, data = mat, REML = FALSE)
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  subtype_results <- coef_table['freq_r2_comp', , drop = FALSE]
  subtype_confint <- confint_table['freq_r2_comp', , drop = FALSE]
  combined_results <- data.frame(
    Celltypes = subtype,
    Estimate = subtype_results[1],
    StdError = subtype_results[2],
    tValue = subtype_results[4],
    pValue = coef_table['freq_r2_comp', 'Pr(>|t|)'],
    CI_lower = subtype_confint[1],
    CI_upper = subtype_confint[2],
    dataset = ds
  )
  return(combined_results)
})

res_df <- do.call(rbind, uni_ds)
res_subtype_all <- res_lmer |> filter(Celltypes == subtype) |> select(!fdr_cat)
res_subtype_all$fdr <- 'Overall'
names(res_subtype_all)[[8]] <- 'dataset'
res_df <- res_df |> filter(!dataset %in% c('HNSC_Luoma','SKCM_Becker'))
res_df <- rbind(res_subtype_all, res_df)
order <- res_df$dataset
res_df <- res_df |> replace_na(list(CI_lower = -20, CI_upper = 10)) 
res_df$CI_lower[res_df$CI_lower < -20] <- -20
res_df$CI_upper[res_df$CI_upper > 10] <- 10

pdf(paste0('/bigdata/zlin/Melanoma_meta/figures/Change/', subtype, '.pdf'), height = 3.5, width = 4.5)
p <- ggplot(res_df, aes(x= factor(dataset, levels = rev(order)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, color = '#A4DBFF', position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.05), -log10(0.1), -log10(0.5)), labels = c('<0.05', '<0.1', '<0.5'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-20, 5), breaks = c(-15, -10, -5, 0, 5), labels = c("-15", "-10", "-5", "0", "5")) + 
  coord_flip() + theme_minimal() + ggtitle(unique(res_df$Celltypes)) +
  labs(size = "P value") + xlab("") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.37, "npc"),
  x1 = unit(0.27, "npc"),
  y = unit(0.24, "npc"),
  y1 = unit(0.24, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.68, "npc"),
  x1 = unit(0.78, "npc"),
  y = unit(0.24, "npc"),
  y1 = unit(0.24, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.32, "npc"), y = unit(0.21, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.73, "npc"), y = unit(0.21, "npc"), gp = gpar(fontsize = 8))
dev.off()

subtype <- 'CD4_TfhTh1_IFNG'
res_lmer |> filter(Celltypes == subtype)
freq_mat <- meta_combi |> 
  filter(celltype_r2 == subtype) |>
  distinct(sample, celltype_r2, .keep_all = T) 

uni_ds <- lapply(setdiff(unique(freq_mat$dataset), 'PCa_Hawley'), function(ds){
  mat <- freq_mat[freq_mat$dataset == ds,]
  print(ds)
  mat$timepoint <- ifelse(mat$time_point == 'Pre', 0, 1)
  variables_to_use <- c('freq_r2_comp', 'response'[length(unique(mat$response)) > 1])
  variables_to_use <- c(variables_to_use, "(1 | patient)")
  formula <- as.formula(paste("timepoint ~", paste((variables_to_use), collapse = " + ")))
  model <- lmer(formula, data = mat, REML = FALSE)
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  subtype_results <- coef_table['freq_r2_comp', , drop = FALSE]
  subtype_confint <- confint_table['freq_r2_comp', , drop = FALSE]
  combined_results <- data.frame(
    Celltypes = subtype,
    Estimate = subtype_results[1],
    StdError = subtype_results[2],
    tValue = subtype_results[4],
    pValue = coef_table['freq_r2_comp', 'Pr(>|t|)'],
    CI_lower = subtype_confint[1],
    CI_upper = subtype_confint[2],
    dataset = ds
  )
  return(combined_results)
})

res_df <- do.call(rbind, uni_ds)
res_subtype_all <- res_lmer |> filter(Celltypes == subtype) |> select(!fdr_cat)
res_subtype_all$fdr <- 'Overall'
names(res_subtype_all)[[8]] <- 'dataset'
res_df <- rbind(res_subtype_all, res_df)
order <- res_df$dataset
res_df$CI_lower[res_df$CI_lower < -15] <- -15
res_df$CI_upper[res_df$CI_upper > 20] <- 20
pdf(paste0('/bigdata/zlin/Melanoma_meta/figures/Change/', subtype, '.pdf'), height = 4.5, width = 4.5)
p <- ggplot(res_df, aes(x= factor(dataset, levels = rev(order)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, color = '#A4DBFF', position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.05), -log10(0.1), -log10(0.5)), labels = c('<0.05', '<0.1', '<0.5'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-15.1, 20.1), breaks = c(-15, -10, -5, 0, 5, 10, 15, 20), labels = c("<-15", "-10", "-5", "0", "5", "10", "15", ">20")) + 
  coord_flip() + theme_minimal() + ggtitle(unique(res_df$Celltypes)) +
  labs(size = "P value") + xlab("") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.35, "npc"),
  x1 = unit(0.25, "npc"),
  y = unit(0.18, "npc"),
  y1 = unit(0.18, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.7, "npc"),
  x1 = unit(0.8, "npc"),
  y = unit(0.18, "npc"),
  y1 = unit(0.18, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.32, "npc"), y = unit(0.16, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.73, "npc"), y = unit(0.16, "npc"), gp = gpar(fontsize = 8))
dev.off()

subtype <- 'gdT'
res_lmer |> filter(Celltypes == subtype)
freq_mat <- meta_combi |> 
  filter(celltype_r2 == subtype) |>
  distinct(sample, celltype_r2, .keep_all = T) 

uni_ds <- lapply(setdiff(unique(freq_mat$dataset), c('PCa_Hawley','SCC_Yost')), function(ds){
  mat <- freq_mat[freq_mat$dataset == ds,]
  print(ds)
  mat$timepoint <- ifelse(mat$time_point == 'Pre', 0, 1)
  variables_to_use <- c('freq_r2_comp', 'response'[length(unique(mat$response)) > 1])
  variables_to_use <- c(variables_to_use, "(1 | patient)")
  formula <- as.formula(paste("timepoint ~", paste((variables_to_use), collapse = " + ")))
  model <- lmer(formula, data = mat, REML = FALSE)
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  subtype_results <- coef_table['freq_r2_comp', , drop = FALSE]
  subtype_confint <- confint_table['freq_r2_comp', , drop = FALSE]
  combined_results <- data.frame(
    Celltypes = subtype,
    Estimate = subtype_results[1],
    StdError = subtype_results[2],
    tValue = subtype_results[4],
    pValue = coef_table['freq_r2_comp', 'Pr(>|t|)'],
    CI_lower = subtype_confint[1],
    CI_upper = subtype_confint[2],
    dataset = ds
  )
  return(combined_results)
})

res_df <- do.call(rbind, uni_ds)
res_subtype_all <- res_lmer |> filter(Celltypes == subtype) |> select(!fdr_cat)
res_subtype_all$fdr <- 'Overall'
names(res_subtype_all)[[8]] <- 'dataset'
res_df <- rbind(res_subtype_all, res_df)
order <- res_df$dataset
res_df <- res_df |> filter(!dataset %in% c('TNBC_Zhang','BCC/SCC_Yost','TNBC_Shiao'))
res_df$CI_lower[res_df$CI_lower < -15] <- -15
res_df$CI_upper[res_df$CI_upper > 20] <- 20

pdf(paste0('/bigdata/zlin/Melanoma_meta/figures/Change/', subtype, '.pdf'), height = 4, width = 4.5)
p <- ggplot(res_df, aes(x= factor(dataset, levels = rev(order)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, color = '#A4DBFF', position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.05), -log10(0.1), -log10(0.5)), labels = c('<0.05', '<0.1', '<0.5'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-15.1, 20.1), breaks = c(-15, -10, -5, 0, 5, 10, 15, 20), labels = c("<-15", "-10", "-5", "0", "5", "10", "15", ">20")) + 
  coord_flip() + theme_minimal() + ggtitle(unique(res_df$Celltypes)) +
  labs(size = "P value") + xlab("") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.35, "npc"),
  x1 = unit(0.25, "npc"),
  y = unit(0.18, "npc"),
  y1 = unit(0.18, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.7, "npc"),
  x1 = unit(0.8, "npc"),
  y = unit(0.18, "npc"),
  y1 = unit(0.18, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.32, "npc"), y = unit(0.16, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.73, "npc"), y = unit(0.16, "npc"), gp = gpar(fontsize = 8))
dev.off()












# Pairwise comparison (frequency)
meta_combi$int_cat <- ifelse(meta_combi$interval < 21, '< 21d', '>= 21d')
# meta_combi$time_point <- factor(meta_combi$time_point, levels = c('Pre','Post'))
meta_combi$response <- factor(meta_combi$response, levels = c('RE','NR'))
meta_combi$time_point2[meta_combi$time_point != 'Pre' & meta_combi$interval < 21] <- 'On(<3w)'
meta_combi$time_point2[meta_combi$time_point != 'Pre' & meta_combi$interval >= 21] <- 'Post(>=3w)'
meta_combi$time_point2[!meta_combi$time_point2 %in% c('On(<3w)', 'Post(>=3w)')] <- 'Pre'
meta_combi$time_point2 <- factor(meta_combi$time_point2, levels = c('Pre', 'On(<3w)', 'Post(>=3w)'))

# Main level
subtype <- 'T'
pw_boxplot <- function(meta_combi, subtype){
  pt <- meta_combi |>
    filter(celltype_main == subtype) |> 
    select(celltype_main, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
    distinct(celltype_main, patient, time_point, .keep_all = T) |>
    pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
    filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5)) |> 
    pull(patient) |> 
    unique()
  
  df <- meta_combi |> 
    filter(celltype_main == subtype, count_r2 >= 5, patient %in% pt) |> 
    distinct(sample, .keep_all = T) |> 
    group_by(patient) |> 
    mutate(pt_count = n()) |> 
    ungroup() |> 
    filter(pt_count == 2) |> 
    mutate(freq_rela = freq_r2/freq_main) 
  
  p <-  ggplot(df, aes(x = time_point2, y = freq_main)) +
    geom_boxplot(color = 'black', outlier.shape = NA) +
    geom_line(aes(group = patient), color = "gray", linetype = 'dashed') +
    geom_point(aes(color = cancertype, shape = modality), size=2, alpha = 0.6) +
    scale_color_manual(values = brewer.pal(length(unique(meta_combi$cancertype)),"Paired")) +
    # facet_wrap(~response, scales = "free") +
    theme_classic2() + xlab("") + ylab("Frequency") + ggtitle(subtype) +
    labs(color = "Cancer Type", shape = "Modality", linetype = "Interval") 
  return(p)
}
ps <- lapply(unique(meta_combi$celltype_main), function(x) pw_boxplot(meta_combi, subtype = x))







subtype <- 'CD4_TfhTh1_IFNG'
pw_boxplot <- function(meta_combi, subtype){
  pt <- meta_combi |>
    filter(celltype_r2 == subtype) |> 
    select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
    distinct(celltype_r2, patient, time_point, .keep_all = T) |>
    pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
    filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5)) |> 
    pull(patient) |> 
    unique()
  meta_combi$time_point[meta_combi$time_point != 'Pre' & meta_combi$interval < 21] <- 'On(<3w)'
  meta_combi$time_point[meta_combi$time_point != 'Pre' & meta_combi$interval >= 21] <- 'On(>=3w)'
  meta_combi$time_point <- factor(meta_combi$time_point, levels = c('Pre', 'On(<3w)', 'Post(>=3w)'))
  
  df <- meta_combi |> 
    filter(celltype_r2 == subtype, response != 'NE', count_r2 >= 5, patient %in% pt) |> 
    distinct(sample, .keep_all = T) |> 
    group_by(patient) |> 
    mutate(pt_count = n()) |> 
    ungroup() |> 
    filter(pt_count == 2) |> 
    mutate(freq_rela = freq_r2/freq_main) 
  
  p <-  ggplot(df, aes(x = time_point, y = freq_rela)) +
    geom_boxplot(color = 'black', outlier.shape = NA) +
    geom_line(aes(group = patient), color = "gray", linetype = 'dashed') +
    geom_point(aes(color = cancertype, shape = modality), size=3, alpha = 0.6) +
    scale_color_manual(values = brewer.pal(length(unique(meta_combi$dataset)),"Set1")) +
    facet_wrap(~response, scales = "free") +
    theme_classic2() + xlab("") + ylab("Frequency") + ggtitle(subtype) +
    labs(color = "Cancer Type", shape = "Modality", linetype = "Interval") 
  return(p)
}

ps <- lapply(c('cDC1', 'CD4_Treg_TNFRSF9', 'CD8_Tm_IL7R', 'CD4_TfhTh_IFNG', 'Macro_C1QC'), function(x) pw_boxplot(meta_combi, subtype = x))
pw_boxplot(meta_combi, subtype = 'Macro_C1QC')

# Univariate logistic
meta_combi$int_cat <- ifelse(meta_combi$interval < 21, 'Short','Long')
freq_mat <- distinct(meta_combi, sample, celltype_r2, .keep_all = T)
uni_logi <- function(freq_mat , meta_combi, n = 20, covariates = c('freq_r2_comp', 'interval', 'cancertype', 'modality', 'prior')){
  df_check <- meta_combi |>
    filter(!response == 'NE') |> 
    select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
    distinct(celltype_r2, patient, time_point, .keep_all = T) |>
    pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
    filter(abs(Pre-Post) >= 3, (Pre >= 5 | Post >= 5))
  subtypes <- df_check |>
    group_by(celltype_r2, response) |> 
    mutate(n_res = n()) |> 
    filter(response == 'RE' & n_res >= n | response == 'NR' & n_res >= n) |>
    distinct(celltype_r2, .keep_all = T) |> 
    select(celltype_r2, response) |> 
    group_by(celltype_r2) |> 
    summarize(n = n()) |> 
    filter(n == 2) |> 
    pull(celltype_r2)
  uni_models <- lapply(subtypes, function(subtype){
    print(subtype)
    pt <- filter(df_check, celltype_r2 == subtype) |> pull(patient)
    df <- filter(freq_mat, celltype_r2 == subtype, patient %in% pt, time_point == 'Pre', !response == 'NE') 
    df$response <- ifelse(df$response == 'RE', 1, 0)
    # formula <- as.formula(paste("response ~ freq_r2_comp + int_cat + cancertype + modality"))
    formula <- as.formula(paste("response ~ freq_r2_comp + int_cat + modality"))
    model <- glm(formula, data = df, family = binomial)
    coef_table <- coef(summary(model))
    confint_table <- confint(model, level = 0.95)
    subtype_results <- coef_table['freq_r2_comp', , drop = FALSE]
    subtype_confint <- confint_table['freq_r2_comp', , drop = FALSE]
    combined_results <- data.frame(
      Celltypes = subtype,
      Estimate = subtype_results[1],
      StdError = subtype_results[2],
      pValue = coef_table['freq_r2_comp', 'Pr(>|z|)'],
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
res_logi <- uni_logi(freq_mat, meta_combi); res_logi
res_filt <- filter(res_logi, pValue < 0.1) |> arrange(desc(Estimate))
order_row <- res_filt$Celltypes

pdf('/bigdata/zlin/Melanoma_meta/figures/uni_logi_freq_pre.pdf', height = 6, width = 5)
p <- ggplot(res_filt, aes(x= factor(Celltypes, levels = rev(order_row)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  #specify position here
  geom_linerange(size=1,position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  #specify position here too
  geom_point(shape=21, colour="white", fill = '#824362', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.01), -log10(0.05), -log10(0.1)), labels = c('<0.01', '<0.05', '<0.1'), range = c(0,6)) +
  scale_y_continuous(name= "Log(Odd Ratio)", limits = c(-105, 105)) + xlab("") +
  coord_flip() + 
  theme_minimal() +
  labs(size = "P value") +
  guides(size = guide_legend(override.aes = list(fill = "black"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))

grid.draw(p)
grid.segments(
  x = unit(0.4, "npc"), 
  x1 = unit(0.3, "npc"), 
  y = unit(0.13, "npc"), 
  y1 = unit(0.13, "npc"), 
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
# Add right arrow
grid.segments(
  x = unit(0.7, "npc"), 
  x1 = unit(0.8, "npc"), 
  y = unit(0.13, "npc"), 
  y1 = unit(0.13, "npc"), 
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Non-responder", x = unit(0.35, "npc"), y = unit(0.11, "npc"), gp = gpar(fontsize = 8))
grid.text("Responder", x = unit(0.75, "npc"), y = unit(0.11, "npc"), gp = gpar(fontsize = 8))
dev.off()

# # test
# subtype <- 'EC_lymphatic'
# pt <- df_check %>% filter(celltype_r2 == subtype) %>% pull(patient)
# df <- meta_r2 %>% filter(celltype_r2 == subtype, patient %in% pt, response %in% c('RE','NR'))
# model <- summary(glm(response ~  Pre + lfc + interval + cancertype + prior + modality, data = df, family = binomial))
# coef(model)

uni_logi <- function(meta_combi, change, cutoff){
  #
  freq_r2 <- pivot_wider(meta_combi, names_from = time_point, values_from = fraction_r2, id_cols = c('patient', 'celltype_r2')) %>% 
    mutate(
      Pre = map_dbl(Pre, ~ (.x[[1]] %||% 0) + 1e-4),
      Post = map_dbl(Post, ~ (.x[[1]] %||% 0) + 1e-4),
      lfc = log2(Post/Pre),
      diff = Post - Pre)
  #
  meta_r2 <- meta_combi %>% 
    distinct(patient, celltype_r2, .keep_all = T) %>% select(-c(time_point, fraction_r2, count_sample:fraction_r2)) %>% 
    left_join(freq_r2, by = c("patient", "celltype_r2")) 
  meta_r2$modality <- as.factor(meta_r2$modality)
  meta_r2$interval <- as.numeric(meta_r2$interval)
  meta_r2$cancertype <- as.factor(meta_r2$cancertype)
  meta_r2$res_metric <- as.factor(meta_r2$res_metric)
  meta_r2$prior <- as.factor(meta_r2$prior)
  meta_r2$response <- as.factor(meta_r2$response)
  #
  df_check <- meta_combi %>%
    select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior, modality) %>%
    distinct(celltype_r2, patient, time_point, .keep_all = T) %>%
    pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) %>% 
    filter(Pre >= 5| Post >= 5, response != 'NE') 
  #
  valid_subtypes <- df_check %>%
    group_by(celltype_r2) %>%
    filter(length(unique(response)) == 2) %>%
    pull(celltype_r2) %>%
    unique()
  # A data frame to store the results
  results <- data.frame(
    Subtype = character(),
    Pre_coef = numeric(),
    Pre_pval = numeric(),
    change_coef = numeric(),
    change_pval = numeric(),
    stringsAsFactors = FALSE
  )
  
  variables <- c('interval', 'cancertype', 'prior', 'modality')
  for (subtype in valid_subtypes) {
    print(subtype)
    # Filter the patient data based on subtype
    pt <- df_check %>% filter(celltype_r2 == subtype) %>% pull(patient)
    
    # Filter the meta data based on subtype and patient
    df <- meta_r2 %>% filter(celltype_r2 == subtype, patient %in% pt, response %in% c('RE', 'NR'))
    
    # Check if there are any variables in df with only one unique value and remove them
    variables_to_use <- c('Pre', change, variables[sapply(df[variables], function(x) length(unique(x))) > 1])
    
    # If no variables are left after the check, print a message and skip to the next iteration
    if (length(variables_to_use) == 0) {
      warning("No variables left to use in glm for subtype: ", subtype)
      next
    }
    
    # Create the formula for glm
    formula_glm <- as.formula(paste("response ~", paste((variables_to_use), collapse = " + ")))
    
    model <- glm(formula_glm, data = df, family = binomial)
    summary_model <- summary(model)
    
    # Extract the coefficients
    coefs <- coef(summary_model)
    
    # Extract the coefficients and initialize default values
    pre_coef <- NA
    pre_pval <- NA
    change_coef <- NA
    change_pval <- NA
    
    # Check if 'Pre' is in the coefficients and assign values
    if ("Pre" %in% rownames(coefs)) {
      pre_coef <- coefs["Pre", "Estimate"]
      pre_pval <- coefs["Pre", "Pr(>|z|)"]
    }
    
    # Check if 'lfc' is in the coefficients and assign values
    if (change %in% rownames(coefs)) {
      change_coef <- coefs[change, "Estimate"]
      change_pval <- coefs[change, "Pr(>|z|)"]
    }
    
    # Save the results in the data frame
    results <- rbind(results, data.frame(
      Subtype = subtype,
      Pre_coef = pre_coef,
      Pre_pval = pre_pval,
      change_coef = change_coef,
      change_pval = change_pval
    ))
  }
  results <- results %>% filter(Pre_pval < cutoff | change_pval < cutoff) %>% 
    mutate(Category = case_when(
    Pre_pval < cutoff & change_pval < cutoff ~ "Both < 0.1",
    Pre_pval < cutoff ~ "Pre < 0.1",
    change_pval < cutoff  ~ "Change < 0.1"
  ))
  return(results)
}
res_lfc <- uni_logi(meta_combi, 'lfc', 0.1)
res_diff <- uni_logi(meta_combi, 'diff', 0.1)



# A data frame to store the results
results <- data.frame(
  Subtype = character(),
  Pre_coef = numeric(),
  Pre_pval = numeric(),
  change_coef = numeric(),
  change_pval = numeric(),
  stringsAsFactors = FALSE
)

variables <- c('interval', 'cancertype', 'prior', 'modality')
for (subtype in valid_subtypes) {
  print(subtype)
  # Filter the patient data based on subtype
  pt <- df_check %>% filter(celltype_r2 == subtype) %>% pull(patient)
  
  # Filter the meta data based on subtype and patient
  df <- meta_r2 %>% filter(celltype_r2 == subtype, patient %in% pt, response %in% c('RE', 'NR'))
  
  # Check if there are any variables in df with only one unique value and remove them
  variables_to_use <- c('Pre', 'lfc', variables[sapply(df[variables], function(x) length(unique(x))) > 1])
  
  # If no variables are left after the check, print a message and skip to the next iteration
  if (length(variables_to_use) == 0) {
    warning("No variables left to use in glm for subtype: ", subtype)
    next
  }
  
  # Create the formula for glm
  formula_glm <- as.formula(paste("response ~", paste((variables_to_use), collapse = " + ")))
  
  model <- glm(formula_glm, data = df, family = binomial)
  summary_model <- summary(model)
  
  # Extract the coefficients
  coefs <- coef(summary_model)
  
  # Extract the coefficients and initialize default values
  pre_coef <- NA
  pre_pval <- NA
  lfc_coef <- NA
  lfc_pval <- NA
  
  # Check if 'Pre' is in the coefficients and assign values
  if ("Pre" %in% rownames(coefs)) {
    pre_coef <- coefs["Pre", "Estimate"]
    pre_pval <- coefs["Pre", "Pr(>|z|)"]
  }
  
  # Check if 'lfc' is in the coefficients and assign values
  if ("lfc" %in% rownames(coefs)) {
    lfc_coef <- coefs["lfc", "Estimate"]
    lfc_pval <- coefs["lfc", "Pr(>|z|)"]
  }
  
  # Save the results in the data frame
  results <- rbind(results, data.frame(
    Subtype = subtype,
    Pre_coef = pre_coef,
    Pre_pval = pre_pval,
    lfc_coef = lfc_coef,
    lfc_pval = lfc_pval
  ))
}

# Filter the results based on p-value
filtered_results <- results %>% filter(Pre_pval < 0.1 | lfc_pval < 0.1)
filtered_results


























df_obs <- read.csv('/bigdata/zlin/Melanoma_meta/data/obs_integration.csv')
df_obs$patient <- paste0(df_obs$dataset, '_', df_obs$patient)
df_clin <- read.csv('/bigdata/zlin/Melanoma_meta/tables/int.csv')
df_obs$interval <- df_clin$interval[match(df_obs$patient, df_clin$patient)]
df_obs$time_point <- as.character(df_obs$time_point)
df_obs <- filter(df_obs, !patient %in% c('BCC_Yost_su012','BCC_Yost_su009'))
df_obs$interval <- as.numeric(df_obs$interval)
df_obs$time_point[df_obs$time_point %in% c('On','Post') & df_obs$interval <= 21] <- 'On'
df_obs$int_cat <- ifelse(df_obs$interval <= 21, 'Early','Late/Post')
df_obs$int_cat[df_obs$dataset =='HNSC_Franken'] <- 'Early'
df_obs$sample <- paste0(df_obs$dataset, '_', df_obs$sample)

# modality
df_obs$modality <- df_obs$treatment
df_obs$modality <- ifelse(df_obs$modality %in% c("aPD1","aPDL1","Anti-PD-L1+ Chemo","aPD1(pre-Chemo)"), 'Mono-immune', 'Comb-immune')
df_obs$response <- as.character(df_obs$response)
df_obs$response[df_obs$dataset %in% c('SKCM_Becker') & df_obs$response =='NE'] <- 'Not available'
df_obs$response[df_obs$dataset %in% c('HNSC_Franken')] <- 'Not available'
df_obs$response_cat <- ifelse(df_obs$response %in% c('PD','SD','NE','Non-response'), 'Non-responder',
                              ifelse(df_obs$response %in% c('Not available', 'non-pCR'), 'Other', 'Responder'))

df_obs$subset <- paste0(df_obs$dataset, '(', df_obs$treatment, ')')
df_obs$celltype_major <- as.character(df_obs$celltype_major)
df_obs$celltype_major[str_split(df_obs$celltype_r2, '\\.', simplify = T)[,1] =='CD8'] <- 'CD8'
df_obs$celltype_major[str_split(df_obs$celltype_r2, '\\.', simplify = T)[,1] =='CD4'] <- 'CD4'

# Single time point comparison (Univariate Logistic Models)
# Here samples were grouped by response
# Fraction for subsets
# Pre-treatment
timepoint <- 'Post'
df_fraction <- df_obs |> 
  dplyr::filter(time_point == timepoint, response_cat %in% c("Non-responder", "Responder")) |>
  group_by(sample) |> 
  mutate(count_sample = n()) |>
  group_by(celltype_r2, .add = TRUE)  |>  
  mutate(count_subset = n())  |>  
  mutate(fraction = count_subset/count_sample)  |>  
  distinct(patient, celltype_r2, .keep_all = T)

df_fraction$response_cat <- factor(df_fraction$response_cat, levels = c("Non-responder", "Responder"))
data <- df_fraction |> filter(celltype_r2  == 'CAF_myo', response_cat != 'Other')
data$subtype <- as.factor(data$subtype)
data$modality <- as.factor(data$modality)
df_fraction$subtype <- as.factor(df_fraction$subtype)
df_fraction$modality <- as.factor(df_fraction$modality)
# df_fraction$treatment <- as.factor(df_fraction$treatment)

perform_glm <- function(data) {
  print(unique(data$celltype_r2))
  # Find factors with less than 2 levels and exclude them from the formula
  factors_with_enough_levels <- sapply(data, function(x) is.factor(x) && length(unique(x)) > 1)
  factors_with_enough_levels <- names(factors_with_enough_levels[factors_with_enough_levels]) |> setdiff('response_cat')
  # Build the formula dynamically
  formula <- as.formula(paste("response_cat ~ fraction + interval +", paste(factors_with_enough_levels, collapse = " + ")))
  model <- glm(formula, data = data, family = binomial)
  coefficients <- summary(model)$coefficients
  if (coefficients["fraction", "Pr(>|z|)"] < 0.05) {
    res <- data.frame(Estimate = coefficients["fraction", "Estimate"],
                          Std_Error = coefficients["fraction", "Std. Error"],
                          Z_value = coefficients["fraction", "z value"],
                          p_value = coefficients["fraction", "Pr(>|z|)"])
    ci <- confint(model, level = 0.95)
    res$CI_lower = ci["fraction", 1] 
    res$CI_upper = ci["fraction", 2] 
    return(res)
  } else {
    res <- data.frame(Estimate = NA,
                      Std_Error = NA,
                      Z_value = NA,
                      p_value = NA,
                      CI_lower = NA,
                      CI_upper = NA)
    return(res)
  }
}

res <- df_fraction |> 
  filter(response_cat != 'Other') |> 
  group_by(celltype_r2) |> 
  do(perform_glm(.)) |> 
  ungroup() |> 
  drop_na() |> 
  arrange(-as.numeric(Estimate)) 
# Forest Plot
p <- 
  res |>
  ggplot(aes(y = reorder(celltype_r2, Estimate))) + 
  theme_classic() +
  geom_point(aes(x = Estimate), shape=15, size=3) +
  geom_linerange(aes(xmin=CI_lower, xmax=CI_upper)) +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Coefficient(Log Odds Ratio)", y="") + coord_cartesian(ylim=c(1,14), xlim=c(-100, 800)) + 
  theme(axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.title.y= element_blank())
p

pdf(paste0('/bigdata/zlin/Melanoma_meta/figures/lolliplot_', timepoint,'.pdf'), width = 5.5, height = 5)
# Create the ggplot object
plot <- ggplot(filtered_results, aes(x = reorder(celltype_r2, t_value), y = t_value)) +
  geom_point(aes(color = significance, size = size)) +
  geom_segment(aes(xend = reorder(celltype_r2, -t_value), yend = 0), color = "grey") +
  geom_hline(yintercept = 0) +
  scale_color_identity() +
  scale_color_manual(values = c("<0.05" = "red", ">=0.05" = "black"), name = "P adjusted") +
  scale_size_continuous(name = "P value",
                        breaks = breaks_size,
                        labels = labels_p_value,) +
  coord_flip() +
  labs(title = "", x = "", y = "") +
  theme_classic() + 
  theme(axis.title.y = element_text(vjust = -8),
        plot.margin = margin(10, 10, 40, 10)) # Adjust the bottom margin
# Draw the ggplot
grid.draw(plot)
# Add left arrow
grid.segments(
  x = unit(0.45, "npc"), 
  x1 = unit(0.35, "npc"), 
  y = unit(0.14, "npc"), 
  y1 = unit(0.14, "npc"), 
  arrow = arrow(type = "open", length = unit(0.1, "inches"))
)

# Add right arrow
grid.segments(
  x = unit(0.7, "npc"), 
  x1 = unit(0.8, "npc"), 
  y = unit(0.14, "npc"), 
  y1 = unit(0.14, "npc"), 
  arrow = arrow(type = "open", length = unit(0.1, "inches"))
)
grid.text("RE", x = unit(0.75, "npc"), y = unit(0.12, "npc"), gp = gpar(fontsize = 8))
grid.text("NR", x = unit(0.4, "npc"), y = unit(0.12, "npc"), gp = gpar(fontsize = 8))
dev.off()

# Boxplot
df_fraction |> 
  filter(time_point == timepoint,
         celltype_r2 %in% (filtered_results |> filter(significance == '<0.05') |> pull(celltype_r2)),
         response_cat != 'Other') |> 
  ggplot(aes(x=celltype_r2, y=fraction, fill=response_cat)) + 
  geom_boxplot() +
  ylim(0, 0.2) + theme_classic() +
  # theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8, color = "black")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8, color = "black")) +
  xlab(NULL) + ylab('Fraction') +
  scale_fill_manual(name = "Response", 
                      values = brewer.pal(1,'Set1')[1:2],
                      labels = c("RE", "NR"))
ggsave(filename = paste0('/bigdata/zlin/Melanoma_meta/figures/boxplot_', timepoint,'.pdf'), width = 3, height = 3)

# Heatmap
df_fraction <- df_obs |> 
  dplyr::filter(response_cat %in% c("Non-responder", "Responder")) |>
  group_by(sample) |> 
  mutate(count_sample = n()) %>% 
  group_by(celltype_r2, .add = TRUE)  |>  
  mutate(count_subset = n())  |>  
  mutate(fraction = count_subset/count_sample)|>  
  distinct(sample, celltype_r2, .keep_all = T) |> 
  select(treatment, response_cat, celltype_r2, patient, fraction, time_point, sample, dataset, ) |> 
  pivot_wider(values_from = fraction, names_from = celltype_r2) 
df_fraction$time_point <- factor(df_fraction$time_point, levels = c('Pre', 'On', 'Post'))
df_fraction <- df_fraction[with(df_fraction, order(df_fraction$time_point, df_fraction$dataset)),]
df_fraction[is.na(df_fraction)] <- 0
mat <- df_fraction[,-c(1:6)] |> t()
mat[mat>0.05] <- 0.05
# heatmap
# t_vec <- brewer.pal(9, "Set1")[1:length(unique(df_fraction$time_point))]
# names(t_vec) <- unique(df_fraction$time_point)
d_vec <- dittoColors()[1:length(unique(df_fraction$dataset))]
names(d_vec) <- unique(df_fraction$dataset)
# ha_col <- HeatmapAnnotation(`Dataset` = df_fraction$dataset, `Time Point` = df_fraction$time_point, col = list(Dataset = d_vec, `Time Point` = t_vec))
ha_col <- HeatmapAnnotation(`Dataset` = df_fraction$dataset, col = list(Dataset = d_vec))
Heatmap(mat, name = 'Fraction', top_annotation = ha_col, cluster_columns = F, column_split = df_fraction$time_point, row_names_gp = grid::gpar(fontsize = 8))




