rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','escape','qs','ggplot2','rstatix','limma')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

GS.hallmark <- getGeneSets(library = "H")
get_ssgsea <- function(datasets, celltype){
  list_ssgsea <- lapply(datasets, function(dataset){
    print(dataset)
    seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs'))
    seu <- subset(seu, subset = (celltype_main == celltype))
    if (celltype == 'CD8+T'){
      seu <- subset(seu, subset = celltype_r2 %in% c('MAIT','gdT'), invert = T). 
    } 
    seu <- runEscape(seu,
                     method = "ssGSEA",
                     gene.sets = GS.hallmark,
                     min.size = 5,
                     new.assay.name = "escape.ssGSEA")
    seu <- performNormalization(seu,
                                assay = "escape.ssGSEA", 
                                gene.sets = GS.hallmark, 
                                scale.factor = seu$nFeature_RNA,
                                make.positive = T)
    df <- cbind(seu$sample, data.frame(t(seu@assays$escape.ssGSEA_normalized$data), check.names = F)) 
    names(df)[1] <- 'sample'
    df <- df |> 
      group_by(sample) |> 
      summarize_all(mean, .groups = 'drop')
    df_sample <- seu@meta.data |> 
      select(sample, patient, time_point, interval, response, dataset, cancertype) |> 
      distinct(sample, .keep_all = T)
    df <- merge(df_sample,df, by = 'sample')
    return(df)
  })
  df_ssgsea <- do.call(rbind, list_ssgsea)
  write.csv(df_ssgsea, paste0('/bigdata/zlin/Melanoma_meta/tables/ssgsea_', celltype, '.csv'))
  return(df_ssgsea)
}

datasets <- c("SKCM_Becker", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "BCC_Yost", "SCC_Yost", "HNSC_IMCISION", "HNSC_Luoma", "NSCLC_Liu", "CRC_Li", "PCa_Hawley")
df_cd4t <- get_ssgsea(datasets, celltype = 'CD4+T')
df_cd8t <- get_ssgsea(datasets, celltype = 'CD8+T')

datasets <- c("SKCM_Becker", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang", "BCC_Yost", "HNSC_IMCISION", "HNSC_Luoma", "CRC_Li", "PCa_Hawley")
df_nk <- get_ssgsea(datasets, celltype = 'NK')
df_cdc <- get_ssgsea(datasets, celltype = 'cDC')
df_mono <- get_ssgsea(datasets, celltype = 'Mono')
df_macro <- get_ssgsea(datasets, celltype = 'Macro')

datasets <- c("SKCM_Becker", "BRCA_Bassez1", "BRCA_Bassez2", "BCC_Yost", "CRC_Li", "PCa_Hawley")
df_endo <- get_ssgsea(datasets, celltype = 'Endo')
df_caf <- get_ssgsea(datasets, celltype = 'CAF')

celltype <- 'Endo'
df <- read.csv(paste0('/bigdata/zlin/Melanoma_meta/tables/ssgsea_', celltype, '.csv'), check.names = F)[,-1]
df <- df |> 
  mutate(int_cat = case_when(interval < 21 ~ 'short',
                             interval >= 21 ~ 'long'),
         time_point = factor(time_point, levels = c('Pre','Post'))) |> 
  group_by(patient) |> 
  mutate(n = n()) |> 
  filter(n == 2) |> 
  data.frame()
pathways <- colnames(df)[str_detect(colnames(df), 'HALLMARK')]
list_res <- lapply(pathways, function(pathway){
  df_sub <- df[,c('time_point', 'patient', 'int_cat', pathway)]
  df_sub$time_point <- factor(df_sub$time_point, levels = c('Post','Pre'))
  ttest <- t_test(df_sub, as.formula(paste0(pathway, ' ~ time_point')), paired = T)
  cohensd <- cohens_d(df_sub, as.formula(paste0(pathway, ' ~ time_point')), paired = T)
  res <- c(pathway, ttest$p, cohensd$effsize)
  return(res)
})
res <- do.call(rbind, list_res) |> data.frame()
colnames(res) <- c('pathway', 'pvalue', 'cohensd')
res <- res |> 
  mutate(pvalue = as.numeric(pvalue),
         cohensd = as.numeric(cohensd)) |> 
  arrange(desc(cohensd), pvalue)

df |> ggplot(aes(time_point, `HALLMARK.MYC.TARGETS.V1`)) + geom_boxplot() +  geom_point() + geom_line(aes(group = patient)) 

patient <- factor(str_replace(df$patient, '\\/', '_'))
dataset <- factor(df$dataset)
timepoint <- factor(df$time_point)
response <- factor(df$response)
design <- model.matrix(~patient + dataset + response + timepoint)
fit <- lmFit(t(df[,-c(1:7, 58:59)]), design)
fit <- eBayes(fit)
res_limma <- topTable(fit, coef = 'timepointPost')






