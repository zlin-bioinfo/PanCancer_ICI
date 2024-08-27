rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','qs', 'lmerTest','grid','msigdbr','UCell')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

datasets <- c('SKCM_Becker', 'BCC_Yost', 'SCC_Yost', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_IMCISION', 'HNSC_Luoma', 'HNSC_Franken', 'NSCLC_Liu', 'CRC_Li','CRC_Chen','PCa_Hawley')
meta_patient <- read.csv('/bigdata/zlin/PanCancer_ICI/tables/meta_patient.csv')
# CD4T
gs_list <- readxl::read_xlsx('/bigdata/zlin/PanCancer_ICI/tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 7, skip = 1) |> 
  lapply(function(x){return(x[!is.na(x)])})
gs_list$Exhaustion <- NULL
names(gs_list)[names(gs_list  ) == 'OXPHOS'] = 'Oxidative phosphorylation'
names(gs_list)[names(gs_list) == 'Lipid metabolism'] = 'Fatty acid metabolism'
gs_list$`Fatty acid metabolism`[gs_list$`Fatty acid metabolism` == 'PLA2G16'] <- 'PLAAT3'
gs_list$`Fatty acid metabolism`[gs_list$`Fatty acid metabolism` == 'ARNTL'] <- 'BMAL1'
gs_list$`Oxidative phosphorylation` <- gsub('\\.', '-', gs_list$`Oxidative phosphorylation`)

list_cd4 <- lapply(datasets, function(dataset){
  seu <- qread(paste0('/bigdata/zlin/PanCancer_ICI/data/', dataset, '/seu_r2.qs')) |> 
    subset(subset = celltype_main == 'CD4+T') |> 
    subset(subset = patient %in% c("CRC_Li_P28", "HNSC_Franken_P16", "HNSC_Franken_P20", "SKCM_Becker_P10", "SKCM_Becker_P6"), invert = T) |> 
    AggregateExpression(group.by = c('sample','time_point','patient'), return.seurat = T)  
  # seu <- AddModuleScore(seu, features = gs_list, name = "FunctionScore")
  seu <- AddModuleScore_UCell(seu, features = gs_list)
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  return(seu@meta.data)
})
df_score <- do.call(rbind, list_cd4)
rownames(df_score) <- NULL
df_score$response <- meta_patient$response[match(gsub('-','_',df_score$patient), meta_patient$patient)]
df_score$interval <- meta_patient$interval[match(gsub('-','_',df_score$patient), meta_patient$patient)]
df_score$int_cat <- ifelse(df_score$interval < 21, '< 21d', '>= 21d')
df_score$modality <- meta_patient$modality[match(gsub('-','_',df_score$patient), meta_patient$patient)]
df_score$time_point <- factor(df_score$time_point, levels = c('Pre', 'On'))
uni_models <- lapply(names(gs_list), function(geneset){
  print(geneset)
  df_sub <- df_score[c('patient', 'time_point', 'response', 'dataset', 'int_cat', 'modality', geneset)]
  df_sub[,length(df_sub)] <- scales::rescale(df_sub[,length(df_sub)], to=c(0, 1))
  df_sub$timepoint <- ifelse(df_sub$time_point == 'Pre', 0, 1)
  colnames(df_sub)[7] <- 'score'
  formula <- as.formula("timepoint ~ score + dataset + int_cat + modality + (1 | patient)")
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
write.csv(res, file = '/bigdata/zlin/PanCancer_ICI/tables/lmer_cd4t.csv', row.names = F)
res$p_cat <- '< 0.05'
res$p_cat[res$pValue>0.05] <- '> 0.05'
pdf('/bigdata/zlin/PanCancer_ICI/figures/Change/uni_lmer_sig_cd4t.pdf', height = 5, width = 5)
p <- ggplot(res, aes(x= factor(Signature, levels = rev(res$Signature)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(color = "#83000A",size=1, position=position_dodge(width = 0.5)) + # 3E5600
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(fill = p_cat), shape=21, color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.5), -log10(0.05)), labels = c('<0.5', '<0.05'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-1.3, 1.3)) +
  xlab("") +
  scale_fill_manual(values = c("#E22D4B", "grey")) + # 009000
  coord_flip() + theme_minimal() +
  labs(title = 'CD4+T', size = "P value", fill = "P value") +
  guides(fill = guide_legend(override.aes = list(size = 5)),
         size = "none") +
  theme(plot.margin = unit(c(0.5, 2, 0.5, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"),
        legend.justification=c(1,0), legend.position=c(1,0))
grid.draw(p)
grid.segments(
  x = unit(0.52, "npc"),
  x1 = unit(0.42, "npc"),
  y = unit(0.12, "npc"),
  y1 = unit(0.12, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.82, "npc"),
  x1 = unit(0.92, "npc"),
  y = unit(0.12, "npc"),
  y1 = unit(0.12, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.47, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 8))
grid.text("On", x = unit(0.87, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 8))
dev.off()

# CD8T
gs_list <- readxl::read_xlsx('/bigdata/zlin/PanCancer_ICI/tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 5, skip = 1) |> 
  lapply(function(x){return(x[!is.na(x)])})
names(gs_list)[names(gs_list) == 'Activation:Effector function'] <- 'Activation/Effector function'
names(gs_list)[names(gs_list) == 'TCR Signaling'] = 'TCR signaling'
names(gs_list)[names(gs_list) == 'IFN Response'] = 'IFN response'
gs_list$`IFN response`[gs_list$`IFN response` == 'DDX58'] <- 'RIGI'
gs_list$`Oxidative phosphorylation` <- gsub('\\.', '-', gs_list$`Oxidative phosphorylation`)
list_cd8 <- lapply(datasets, function(dataset){
  seu <- qread(paste0('/bigdata/zlin/PanCancer_ICI/data/', dataset, '/seu_r2.qs')) |> 
    subset(subset = celltype_main == 'CD8+T') |> 
    subset(subset =  (celltype_r2 %in% c('gdT', 'MAIT')), invert = T) |> 
    subset(subset = patient %in% c("CRC_Li_P28", "HNSC_Franken_P16", "HNSC_Franken_P20", "SKCM_Becker_P10", "SKCM_Becker_P6"), invert = T) |> 
    AggregateExpression(group.by = c('sample','time_point','patient'), return.seurat = T)  
  seu <- AddModuleScore_UCell(seu, features = gs_list)
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  return(seu@meta.data)
})
df_score <- do.call(rbind, list_cd8)
rownames(df_score) <- NULL
df_score$response <- meta_patient$response[match(gsub('-','_',df_score$patient), meta_patient$patient)]
df_score$interval <- meta_patient$interval[match(gsub('-','_',df_score$patient), meta_patient$patient)]
df_score$int_cat <- ifelse(df_score$interval < 21, '< 21d', '>= 21d')
df_score$modality <- meta_patient$modality[match(gsub('-','_',df_score$patient), meta_patient$patient)]
df_score$time_point <- factor(df_score$time_point, levels = c('Pre', 'On'))
uni_models <- lapply(names(gs_list), function(geneset){
  print(geneset)
  df_sub <- df_score[c('patient', 'time_point', 'response', 'dataset', 'int_cat', 'modality', geneset)]
  df_sub[,length(df_sub)] <- scales::rescale(df_sub[,length(df_sub)], to=c(0, 1))
  df_sub$timepoint <- ifelse(df_sub$time_point == 'Pre', 0, 1)
  colnames(df_sub)[7] <- 'score'
  formula <- as.formula("timepoint ~ score + dataset + int_cat + modality + (1 | patient)")
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
write.csv(res, file = '/bigdata/zlin/PanCancer_ICI/tables/lmer_cd8t.csv', row.names = F)
res$p_cat <- '< 0.05'
res$p_cat[res$pValue>0.05] <- '> 0.05'
pdf('/bigdata/zlin/PanCancer_ICI/figures/Change/uni_lmer_sig_cd8t.pdf', height = 5, width = 5)
p <- ggplot(res, aes(x= factor(Signature, levels = rev(res$Signature)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(color = "#344E7A", size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(fill = p_cat), shape=21, color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.5), -log10(0.05)), labels = c('<0.5', '<0.05'), range = c(2,5)) +
  # scale_y_continuous(name= "Effect Size", limits = c(-2.2, 5.8)) + 
  xlab("") +
  scale_fill_manual(values = c("#1F78B4", "grey")) +
  coord_flip() + theme_minimal() +
  labs(title = 'CD8+T', size = "P value", fill = "P value") +
  guides(fill = guide_legend(override.aes = list(size = 5)),
         size = "none") +
  theme(plot.margin = unit(c(0.5, 2, 0.5, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"),
        legend.justification=c(1,0), legend.position=c(1,0))
grid.draw(p)
grid.segments(
  x = unit(0.52, "npc"),
  x1 = unit(0.42, "npc"),
  y = unit(0.11, "npc"),
  y1 = unit(0.11, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(0.82, "npc"),
  x1 = unit(0.92, "npc"),
  y = unit(0.11, "npc"),
  y1 = unit(0.11, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(0.47, "npc"), y = unit(0.09, "npc"), gp = gpar(fontsize = 8))
grid.text("On", x = unit(0.87, "npc"), y = unit(0.09, "npc"), gp = gpar(fontsize = 8))
dev.off()

# Macrophages
# kegg_antigen_proc_pres <- msigdbr(species = "Homo sapiens", category = 'C2',subcategory = 'CP:KEGG') |> 
#   filter(gs_name  == 'KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION') |> 
#   pull(gene_symbol) |> 
#   unique()
ifna <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_INTERFERON_ALPHA_RESPONSE') |>
  pull(gene_symbol) |>
  unique()
ifng <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_INTERFERON_GAMMA_RESPONSE') |>
  pull(gene_symbol) |>
  unique()
ifng[ifng == 'DDX58'] <- 'RIGI'
mhc1 <- msigdbr(species = "Homo sapiens", category = 'C2') |> 
  filter(gs_name  == 'REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION') |> 
  pull(gene_symbol) |> 
  unique()
mhc2 <- msigdbr(species = "Homo sapiens", category = 'C2') |> 
  filter(gs_name  == 'REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION') |> 
  pull(gene_symbol) |> 
  unique()
# ifn1 <- msigdbr(species = "Homo sapiens", category = 'C5') |> 
#   filter(gs_name  == 'GOBP_RESPONSE_TO_TYPE_I_INTERFERON') |> 
#   pull(gene_symbol) |> 
#   unique()
# ifn2 <- msigdbr(species = "Homo sapiens", category = 'C5') |> 
#   filter(gs_name  == 'GOBP_RESPONSE_TO_INTERFERON_GAMMA') |> 
#   pull(gene_symbol) |> 
#   unique()
complement <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_COMPLEMENT') |> 
  pull(gene_symbol) |> 
  unique()
tnfa <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB') |> 
  pull(gene_symbol) |> 
  unique()
tnfa[tnfa == 'DDX58'] <- 'RIGI'
sig_mac <- readxl::read_xlsx('/bigdata/zlin/PanCancer_ICI/tables/1-s2.0-S0092867421000106-mmc5.xlsx', skip = 1)
sig_mac$M1[sig_mac$M1 == 'IL23'] <- 'IL23A'
sig_mac$M2[sig_mac$M2 == 'FASL'] <- 'FASLG'
gs_list <- list(`M1 Signature` = sig_mac$M1[!is.na(sig_mac$M1)],
                `M2 Signature` = sig_mac$M2[!is.na(sig_mac$M2)],
                Angiogenesis = sig_mac$Angiogenesis[!is.na(sig_mac$Angiogenesis)],
                Phagocytosis = sig_mac$Phagocytosis[!is.na(sig_mac$Phagocytosis)],
                # `Antigen Proc.&Pres.` = kegg_antigen_proc_pres,
                `TNF-alpha Signaling` = tnfa,
                `IFN-alpha Response` = ifna,
                `IFN-gamma Response` = ifng,
                Complement = complement,
                `MHC I Antigen Proc.&Pres.` = mhc1,
                `MHC II Antigen Presentation` = mhc2)
                # `Response to Type I IFN` = ifn1,
                # `Response to Type II IFN` = ifn2)
# Dataset names
datasets <- c('SKCM_Becker', 'BCC_Yost', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_IMCISION', 'HNSC_Luoma', 'HNSC_Franken', 'CRC_Li','CRC_Chen','PCa_Hawley')
list_mac <- lapply(datasets, function(dataset){
  seu <- qread(paste0('/bigdata/zlin/PanCancer_ICI/data/', dataset, '/seu_r2.qs')) |> 
    subset(subset = celltype_main == 'Macro') |> 
    subset(subset = patient %in% c("CRC_Li_P28", "HNSC_Franken_P16", "HNSC_Franken_P20", "SKCM_Becker_P10", "SKCM_Becker_P6"), invert = T) |> 
    AggregateExpression(group.by = c('sample','time_point','patient'), return.seurat = T)  
  seu <- AddModuleScore_UCell(seu, features = gs_list)
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  return(seu@meta.data)
})
df_score <- do.call(rbind, list_mac)
rownames(df_score) <- NULL
# df_score$response <- meta_patient$response[match(gsub('-','_',df_score$patient), meta_patient$patient)]
df_score$interval <- meta_patient$interval[match(gsub('-','_',df_score$patient), meta_patient$patient)]
df_score$int_cat <- ifelse(df_score$interval < 21, '< 21d', '>= 21d')
df_score$modality <- meta_patient$modality[match(gsub('-','_',df_score$patient), meta_patient$patient)]
df_score$time_point <- factor(df_score$time_point, levels = c('Pre', 'On'))

uni_models <- lapply(names(gs_list), function(geneset){
  print(geneset)
  df_sub <- df_score[c('patient', 'time_point', 'dataset', 'int_cat', 'modality', geneset)]
  df_sub[,length(df_sub)] <- scales::rescale(df_sub[,length(df_sub)], to=c(0, 1))
  df_sub$timepoint <- ifelse(df_sub$time_point == 'Pre', 0, 1)
  colnames(df_sub)[6] <- 'score'
  formula <- as.formula("timepoint ~ score + dataset + int_cat + modality + (1 | patient)")
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
write.csv(res, file = '/bigdata/zlin/PanCancer_ICI/tables/lmer_mac.csv', row.names = F)
res$p_cat <- '< 0.05'
res$p_cat[res$pValue>0.05] <- '> 0.05'
pdf('/bigdata/zlin/PanCancer_ICI/figures/Change/uni_lmer_sig_macro.pdf', height = 4, width = 4.5)
p <- ggplot(res, aes(x= factor(Signature, levels = rev(res$Signature)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(color = "#3E5600",size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(aes(fill = p_cat), shape=21, color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.5), -log10(0.05)), labels = c('<0.5', '<0.05'), range = c(2,5)) +
  # scale_y_continuous(name= "Effect Size", limits = c(-1.2, 5.5)) + 
  xlab("") +
  scale_fill_manual(values = c("#33A02C", "grey")) +
  coord_flip() + theme_minimal() +
  labs(title = 'Macrophages', size = "P value", fill = "P value") +
  guides(fill = guide_legend(override.aes = list(size = 5)),
         size = "none") +
  theme(plot.margin = unit(c(0.5, 2, 0.5, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"),
        legend.justification=c(1,0), legend.position=c(1,0))
grid.draw(p)
y_text = 0.13
x1 = 0.45
x2 = 0.85
grid.segments(
  x = unit(x1 + 0.05, "npc"),
  x1 = unit(x1 - 0.05, "npc"),
  y = unit(y_text+0.02, "npc"),
  y1 = unit(y_text+0.02, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.segments(
  x = unit(x2 - 0.05, "npc"),
  x1 = unit(x2 + 0.05, "npc"),
  y = unit(y_text+0.02, "npc"),
  y1 = unit(y_text+0.02, "npc"),
  arrow = arrow(type = "open", length = unit(0.05, "inches"))
)
grid.text("Pre", x = unit(x1, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
grid.text("On", x = unit(x2, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
dev.off()




