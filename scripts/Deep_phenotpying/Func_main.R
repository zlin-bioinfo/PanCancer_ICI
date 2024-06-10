#!/usr/bin/env Rscript
rm(list=ls())
pkgs <- c('Seurat','dplyr','tidyr','stringr','ggsci','qs','readxl','RColorBrewer','pheatmap','ggplotify','ggplot2','patchwork','rstatix', 'tibble','lmerTest','forcast','grid','msigdbr','UCell')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

# Patient info
pt_df <- read.csv('/bigdata/zlin/Melanoma_meta/tables/patient_combi.csv')
# Dataset names
datasets <- c("SKCM_Becker", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "BCC_Yost", "SCC_Yost", "HNSC_IMCISION", "HNSC_Luoma", "NSCLC_Liu", "CRC_Li", "PCa_Hawley")
# CD4T
gs_list <- read_xlsx('/bigdata/zlin/Melanoma_meta/tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 7, skip = 1) |> 
  lapply(function(x){return(x[!is.na(x)])})
gs_list$Exhaustion <- NULL
names(gs_list)[names(gs_list  ) == 'OXPHOS'] = 'Oxidative phosphorylation'
names(gs_list)[names(gs_list) == 'Lipid metabolism'] = 'Fatty acid metabolism'
gs_list$`Fatty acid metabolism`[gs_list$`Fatty acid metabolism` == 'PLA2G16'] <- 'PLAAT3'
gs_list$`Fatty acid metabolism`[gs_list$`Fatty acid metabolism` == 'ARNTL'] <- 'BMAL1'
gs_list$`Oxidative phosphorylation` <- gsub('\\.', '-', gs_list$`Oxidative phosphorylation`)
  
list_meta <- lapply(datasets, function(dataset) {
  print(dataset)
  if (!dataset %in% c("HNSC_Luoma", "PCa_Hawley", "NSCLC_Liu")){
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
      subset(subset = celltype_main == 'CD4+T') |> 
      AggregateExpression(group.by = c('sample','patient','time_point','response'), return.seurat = T)  
    seu <- AddModuleScore_UCell(seu, features = gs_list)
  } else if (dataset == "NSCLC_Liu") {
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
      subset(subset = celltype_main == 'CD4+T') |> 
      AggregateExpression(group.by = c('sample','patient','time_point'), return.seurat = T)  
    seu <- AddModuleScore_UCell(seu, features = gs_list)
    seu$response <- "RE"
  } else {
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
      subset(subset = celltype_main == 'CD4+T') |> 
      AggregateExpression(group.by = c('sample','patient','time_point'), return.seurat = T)  
    seu <- AddModuleScore_UCell(seu, features = gs_list)
    seu$response <- "NR"
  }
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  return(seu@meta.data)
})

df_score <- do.call(rbind, list_meta)
rownames(df_score) <- NULL
df_score <- df_score |> filter(!patient %in% c("CRC-Li-P28", "SKCM-Becker-P10", "SKCM-Becker-P6"))
df_score$patient <- gsub('-', '_', df_score$patient)
df_score <- merge(df_score, pt_df[,c('patient', 'interval', 'treatment')], by = 'patient') |> distinct(sample, .keep_all = T)
df_score$int_cat <- ifelse(df_score$interval < 21, '< 21d', '>= 21d')
df_score$modality <- ifelse(df_score$treatment %in% c("aPD1","aPDL1","Anti-PD-L1+ Chemo","aPD1(pre-Chemo)"), 'Mono-immune', 'Comb-immune')
df_score$time_point <- factor(df_score$time_point, levels = c('Pre', 'Post'))
uni_models <- lapply(names(gs_list), function(geneset){
  print(geneset)
  df_sub <- df_score[c('patient', 'time_point', 'response', 'dataset', 'int_cat', 'modality', geneset)]
  df_sub$timepoint <- ifelse(df_sub$time_point == 'Pre', 0, 1)
  colnames(df_sub)[7] <- 'score'
  formula <- as.formula("timepoint ~ score + dataset + response + int_cat + modality + (1 | patient)")
  model <- lmer(formula, data = df_sub, REML = FALSE)
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  geneset_results <- coef_table['score', , drop = FALSE]
  geneset_confint <- confint_table['score', , drop = FALSE]
  combined_results <- data.frame(
    Signature = geneset,
    Estimate = geneset_results[1],
    StdError = geneset_results[2],
    tValue = geneset_results[4],
    pValue = coef_table['score', 'Pr(>|t|)'],
    CI_lower = geneset_confint[1],
    CI_upper = geneset_confint[2]
  )
  })
res <- do.call(rbind,uni_models) |> arrange(desc(Estimate))
res$p_cat <- '< 0.05'
res$p_cat[res$pValue>0.05] <- '> 0.05'
pdf('/bigdata/zlin/Melanoma_meta/figures/Change/uni_lmer_sig_cd4t.pdf', height = 5, width = 5)
p <- ggplot(res, aes(x= factor(Signature, levels = rev(res$Signature)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(color = "#3E5600",size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(fill = p_cat), shape=21, color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.5), -log10(0.05)), labels = c('<0.5', '<0.05'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-2.2, 4.3)) + xlab("") +
  scale_fill_manual(values = c("#009000", "grey")) +
  coord_flip() + theme_minimal() +
  labs(title = 'CD4+T', size = "P value", fill = "P value") +
  guides(fill = guide_legend(override.aes = list(size = 5)),
         size = "none") +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.47, "npc"),
  x1 = unit(0.37, "npc"),
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
grid.text("Pre", x = unit(0.42, "npc"), y = unit(0.15, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.75, "npc"), y = unit(0.15, "npc"), gp = gpar(fontsize = 8))
dev.off()

# CD8T
gs_list <- read_xlsx('/bigdata/zlin/Melanoma_meta/tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 5, skip = 1) |> 
  lapply(function(x){return(x[!is.na(x)])})
names(gs_list)[names(gs_list) == 'Activation:Effector function'] <- 'Activation/Effector function'
names(gs_list)[names(gs_list) == 'TCR Signaling'] = 'TCR signaling'
names(gs_list)[names(gs_list) == 'IFN Response'] = 'IFN response'
gs_list$`IFN response`[gs_list$`IFN response` == 'DDX58'] <- 'RIGI'
gs_list$`Oxidative phosphorylation` <- gsub('\\.', '-', gs_list$`Oxidative phosphorylation`)

list_meta <- lapply(datasets, function(dataset) {
  print(dataset)
  if (!dataset %in% c("HNSC_Luoma", "PCa_Hawley", "NSCLC_Liu")){
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
      subset(subset = (celltype_main == 'CD8+T')) |> 
      subset(subset =  (celltype_r2 %in% c('gdT', 'MAIT')), invert = T) |> 
      AggregateExpression(group.by = c('sample','patient','time_point','response'), return.seurat = T)  
    seu <- AddModuleScore(seu, features = gs_list, name = "FunctionScore")
  } else if (dataset == "NSCLC_Liu") {
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
      subset(subset = (celltype_main == 'CD8+T')) |> 
      subset(subset =  (celltype_r2 %in% c('gdT', 'MAIT')), invert = T) |> 
      AggregateExpression(group.by = c('sample','patient','time_point'), return.seurat = T)  
    seu <- AddModuleScore(seu, features = gs_list, name = "FunctionScore")
    seu$response <- "RE"
  } else {
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
      subset(subset = (celltype_main == 'CD8+T')) |> 
      subset(subset =  (celltype_r2 %in% c('gdT', 'MAIT')), invert = T) |> 
      AggregateExpression(group.by = c('sample','patient','time_point'), return.seurat = T)  
    seu <- AddModuleScore(seu, features = gs_list, name = "FunctionScore")
    seu$response <- "NR"
  }
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0("FunctionScore", i)] <- names(gs_list)[i]
  }
  return(seu@meta.data)
})

df_score <- do.call(rbind, list_meta)
rownames(df_score) <- NULL
df_score <- df_score |> filter(!patient %in% c("CRC-Li-P28", "SKCM-Becker-P10", "SKCM-Becker-P6"))
df_score$patient <- gsub('-', '_', df_score$patient)
df_score <- merge(df_score, pt_df[,c('patient', 'interval', 'treatment')], by = 'patient') |> distinct(sample, .keep_all = T)
df_score$int_cat <- ifelse(df_score$interval < 21, '< 21d', '>= 21d')
df_score$modality <- ifelse(df_score$treatment %in% c("aPD1","aPDL1","Anti-PD-L1+ Chemo","aPD1(pre-Chemo)"), 'Mono-immune', 'Comb-immune')
df_score$time_point <- factor(df_score$time_point, levels = c('Pre', 'Post'))
uni_models <- lapply(names(gs_list), function(geneset){
  print(geneset)
  df_sub <- df_score[c('patient', 'time_point', 'response', 'dataset', 'int_cat', 'modality', geneset)]
  df_sub$timepoint <- ifelse(df_sub$time_point == 'Pre', 0, 1)
  colnames(df_sub)[7] <- 'score'
  formula <- as.formula("timepoint ~ score + dataset + response + int_cat + modality + (1 | patient)")
  model <- lmer(formula, data = df_sub, REML = FALSE)
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  geneset_results <- coef_table['score', , drop = FALSE]
  geneset_confint <- confint_table['score', , drop = FALSE]
  combined_results <- data.frame(
    Signature = geneset,
    Estimate = geneset_results[1],
    StdError = geneset_results[2],
    tValue = geneset_results[4],
    pValue = coef_table['score', 'Pr(>|t|)'],
    CI_lower = geneset_confint[1],
    CI_upper = geneset_confint[2]
  )
})
res <- do.call(rbind,uni_models) |> arrange(desc(Estimate))
res$p_cat <- '< 0.05'
res$p_cat[res$pValue>0.05] <- '> 0.05'
pdf('/bigdata/zlin/Melanoma_meta/figures/Change/uni_lmer_sig_cd8t.pdf', height = 5, width = 5)
p <- ggplot(res, aes(x= factor(Signature, levels = rev(res$Signature)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(color = "#83000A", size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(fill = p_cat), shape=21, color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.5), -log10(0.05)), labels = c('<0.5', '<0.05'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-2.2, 5.8)) + xlab("") +
  scale_fill_manual(values = c("#E22D4B", "grey")) +
  coord_flip() + theme_minimal() +
  labs(title = 'CD8+T', size = "P value", fill = "P value") +
  guides(fill = guide_legend(override.aes = list(size = 5)),
         size = "none") +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.5, "npc"),
  x1 = unit(0.4, "npc"),
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
grid.text("Pre", x = unit(0.45, "npc"), y = unit(0.15, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.75, "npc"), y = unit(0.15, "npc"), gp = gpar(fontsize = 8))
dev.off()

# Macrophages
kegg_antigen_proc_pres <- msigdbr(species = "Homo sapiens", category = 'C2',subcategory = 'CP:KEGG') |> 
  filter(gs_name  == 'KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION') |> 
  pull(gene_symbol) |> 
  unique()
ifng <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_INTERFERON_GAMMA_RESPONSE') |> 
  pull(gene_symbol) |> 
  unique()
ifng[ifng == 'DDX58'] <- 'RIGI'
complement <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_COMPLEMENT') |> 
  pull(gene_symbol) |> 
  unique()
tnfa <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB') |> 
  pull(gene_symbol) |> 
  unique()
tnfa[tnfa == 'DDX58'] <- 'RIGI'
sig_mac <- read_xlsx('/bigdata/zlin/Melanoma_meta/tables/1-s2.0-S0092867421000106-mmc5.xlsx', skip = 1)
sig_mac$M1[sig_mac$M1 == 'IL23'] <- 'IL23A'
sig_mac$M2[sig_mac$M2 == 'FASL'] <- 'FASLG'
gs_list <- list(`M1 Signature` = sig_mac$M1[!is.na(sig_mac$M1)],
                `M2 Signature` = sig_mac$M2[!is.na(sig_mac$M2)],
                Angiogenesis = sig_mac$Angiogenesis[!is.na(sig_mac$Angiogenesis)],
                Phagocytosis = sig_mac$Phagocytosis[!is.na(sig_mac$Phagocytosis)],
                `Antigen Proc.&Pres.` = kegg_antigen_proc_pres,
                `TNF-α Signaling` = tnfa,
                `IFN-γ Response` = ifng,
                Complement = complement)
# Dataset names
datasets <- c("SKCM_Becker", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "BCC_Yost", "HNSC_IMCISION", "HNSC_Luoma", "CRC_Li", "PCa_Hawley")
list_meta <- lapply(datasets, function(dataset) {
  print(dataset)
  if (!dataset %in% c("HNSC_Luoma", "PCa_Hawley", "NSCLC_Liu")){
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
      subset(subset = celltype_main == 'Macro') |> 
      AggregateExpression(group.by = c('sample','patient','time_point','response'), return.seurat = T)  
    seu <- AddModuleScore(seu, features = gs_list, name = "FunctionScore")
  } else if (dataset == "NSCLC_Liu") {
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
      subset(subset = celltype_main == 'Macro') |> 
      AggregateExpression(group.by = c('sample','patient','time_point'), return.seurat = T)  
    seu <- AddModuleScore(seu, features = gs_list, name = "FunctionScore")
    seu$response <- "RE"
  } else {
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) |> 
      subset(subset = celltype_main == 'Macro') |> 
      AggregateExpression(group.by = c('sample','patient','time_point'), return.seurat = T)  
    seu <- AddModuleScore(seu, features = gs_list, name = "FunctionScore")
    seu$response <- "NR"
  }
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0("FunctionScore", i)] <- names(gs_list)[i]
  }
  return(seu@meta.data)
})
df_score <- do.call(rbind, list_meta)
rownames(df_score) <- NULL
df_score <- df_score |> filter(!patient %in% c("CRC-Li-P28", "SKCM-Becker-P10", "SKCM-Becker-P6"))
df_score$patient <- gsub('-', '_', df_score$patient)
df_score <- merge(df_score, pt_df[,c('patient', 'interval', 'treatment')], by = 'patient') |> distinct(sample, .keep_all = T)
df_score$int_cat <- ifelse(df_score$interval < 21, '< 21d', '>= 21d')
df_score$modality <- ifelse(df_score$treatment %in% c("aPD1","aPDL1","Anti-PD-L1+ Chemo","aPD1(pre-Chemo)"), 'Mono-immune', 'Comb-immune')

uni_models <- lapply(names(gs_list), function(geneset){
  print(geneset)
  df_sub <- df_score[c('patient', 'time_point', 'response', 'dataset', 'int_cat', 'modality', geneset)]
  df_sub$timepoint <- ifelse(df_sub$time_point == 'Pre', 0, 1)
  colnames(df_sub)[7] <- 'score'
  formula <- as.formula("timepoint ~ score + dataset + response + int_cat + modality + (1 | patient)")
  model <- lmer(formula, data = df_sub, REML = FALSE)
  model_summary <- summary(model)
  coef_table <- coef(summary(model))
  confint_table <- confint(model, level = 0.95)
  geneset_results <- coef_table['score', , drop = FALSE]
  geneset_confint <- confint_table['score', , drop = FALSE]
  combined_results <- data.frame(
    Signature = geneset,
    Estimate = geneset_results[1],
    StdError = geneset_results[2],
    tValue = geneset_results[4],
    pValue = coef_table['score', 'Pr(>|t|)'],
    CI_lower = geneset_confint[1],
    CI_upper = geneset_confint[2]
  )
})
res <- do.call(rbind,uni_models) |> arrange(desc(Estimate))
res$p_cat <- '< 0.05'
res$p_cat[res$pValue>0.05] <- '> 0.05'
pdf('/bigdata/zlin/Melanoma_meta/figures/Change/uni_lmer_sig_macro.pdf', height = 4, width = 4.5)
p <- ggplot(res, aes(x= factor(Signature, levels = rev(res$Signature)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(color = "#A44E00",size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(fill = p_cat), shape=21, color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.5), -log10(0.05)), labels = c('<0.5', '<0.05'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-1.1, 3)) + xlab("") +
  scale_fill_manual(values = c("#FF9814", "grey")) +
  coord_flip() + theme_minimal() +
  labs(title = 'Macrophages', size = "P value", fill = "P value") +
  guides(fill = guide_legend(override.aes = list(size = 5)),
         size = "none") +
  theme(plot.margin = unit(c(0.5, 1, 2, 1), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
grid.segments(
  x = unit(0.42, "npc"),
  x1 = unit(0.32, "npc"),
  y = unit(0.22, "npc"),
  y1 = unit(0.22, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.67, "npc"),
  x1 = unit(0.77, "npc"),
  y = unit(0.22, "npc"),
  y1 = unit(0.22, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.37, "npc"), y = unit(0.2, "npc"), gp = gpar(fontsize = 8))
grid.text("Post", x = unit(0.72, "npc"), y = unit(0.2, "npc"), gp = gpar(fontsize = 8))
dev.off()

