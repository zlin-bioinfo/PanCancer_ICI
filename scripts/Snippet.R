rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','scCustomize','ggsci','patchwork','ggplot2','gtools','ComplexHeatmap','dittoSeq','RColorBrewer','ggpubr')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

sce_sub <- subset(sce_filt, idents = c('DC', 'Mono/Macro')) %>%
  ScaleData(features = rownames(.)) %>% 
  RunHarmony(group.by.vars='sample') %>%
  RunTSNE(reduction='harmony', dims = 1:30) %>%
  RunUMAP(reduction='harmony', dims = 1:30) %>%
  FindNeighbors(reduction='harmony') %>%
  FindClusters(resolution = 1)

seu%>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(features = rownames(.))

# excluding samples
skcm_becker <- subset(skcm_becker, subset = patient %in% c('P10', 'P4', 'P6', 'P8'), invert = TRUE)
hnsc_franken <- subset(hnsc_franken, subset = patient %in% c('16', '20'), invert = TRUE)
brca_bassez1 <- subset(brca_bassez1, subset = patient %in% c('BIOKEY_23'), invert = TRUE)
bcc_yost <- subset(bcc_yost, subset = patient %in% c('su007','su009','su012'), invert = TRUE)
tnbc_li <- subset(tnbc_li, subset = patient %in% c('P013'), invert = TRUE)

# Datasets
datasets <- list(skcm_becker, hnsc_franken, brca_bassez1, brca_bassez2, bcc_yost, crc_li, tnbc_li)
convert <- function(datasets){
  for (i in 1:length(datasets)){
    datasets[[i]] <- CreateSeuratObject(counts = datasets[[i]]@assays$RNA@counts, meta.data = datasets[[i]]@meta.data)
  }
  return(datasets)
}

# load testing datasets
skcm_becker <- dior::read_h5(file='/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/non-malignant.h5', target.object = 'seurat')
hnsc_franken <- dior::read_h5(file='/bigdata/zlin/Melanoma_meta/data/HNSC_Franken/non-malignant.h5', target.object = 'seurat')
brca_bassez1 <- dior::read_h5(file='/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/non-malignant.h5', target.object = 'seurat')
brca_bassez2 <- dior::read_h5(file='/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/non-malignant.h5', target.object = 'seurat')
bcc_yost <- dior::read_h5(file='/bigdata/zlin/Melanoma_meta/data/BCC_Yost/non-malignant.h5', target.object = 'seurat') 
crc_li <- dior::read_h5(file='/bigdata/zlin/Melanoma_meta/data/CRC_Li/non-malignant.h5', target.object = 'seurat') 
tnbc_li <- dior::read_h5(file='/bigdata/zlin/Melanoma_meta/data/TNBC_Li/non-malignant.h5', target.object = 'seurat') 
hnsc_franken$patient <- as.character(hnsc_franken$patient)
df_obs <- do.call(rbind,
                  list(select(skcm_becker@meta.data, celltype_major, celltype_r2, sample, time_point, patient, response, batch, dataset),
                       select(hnsc_franken@meta.data, celltype_major, celltype_r2, sample, time_point, patient, response, batch, dataset),
                       select(brca_bassez1@meta.data, celltype_major, celltype_r2, sample, time_point, patient, response, batch, dataset),
                       select(brca_bassez2@meta.data, celltype_major, celltype_r2, sample, time_point, patient, response, batch, dataset),
                       select(bcc_yost@meta.data, celltype_major, celltype_r2, sample, time_point, patient, response, batch, dataset),
                       select(crc_li@meta.data, celltype_major, celltype_r2, sample, time_point, patient, response, batch, dataset),
                       select(tnbc_li@meta.data, celltype_major, celltype_r2, sample, time_point, patient, response, batch, dataset)))
df_clin <- read.csv('/bigdata/zlin/Melanoma_meta/tables/interval.csv')
df_obs$treatment <- df_clin$treatment[match(df_obs$batch, df_clin$batch)]
df_obs$interval <- df_clin$interval[match(df_obs$batch, df_clin$batch)]
df_obs$patient <- paste0(df_obs$dataset, '_', df_obs$patient)
df_obs$time_point <- as.character(df_obs$time_point)
df_obs$time_point[df_obs$time_point %in% c('On','Post')] <- "On/Post"
df_obs$modality <- df_obs$treatment
df_obs$modality <- recode(df_obs$modality, 'aPD1+CTLA4' = 'Dual-immune', 'aPDL1+CTLA4' = 'Dual-immune', 'aPD1+COX2i' = 'Immune+inhibitor', 'aPD1+Hedgehog' = 'Immune+inhibitor',
                          'aPDL1' = 'Mono-immune', 'aPD1' = 'Mono-immune', 'aPD1(pre-Chemo)' = 'Mono-immune', 'aPD1+Chemo' = 'Immune+chemo')
df_obs$response <- as.character(df_obs$response)
df_obs$response[df_obs$response == 'NE' & df_obs$dataset =='SKCM_Becker'| df_obs$response == 'n/a'] <- 'NE(unknown)'
df_obs$response[df_obs$response %in% c('PD','NE')] <- 'Non-response'
df_obs$response[df_obs$response %in% c('CR','PR','E','SD')] <- 'Reseponse'
df_obs$subset <- paste0(df_obs$dataset, '(', df_obs$treatment, ')')

# df_obs <- read.csv('/bigdata/zlin/Melanoma_meta/tables/subtype/obs.csv') 

dataset_names <- unique(df_obs$dataset)  # Get unique dataset names

filtered_celltypes_list <- list()  # Empty list to store filtered cell types for each dataset

ROIE_list <- list()
for (p in unique(df_obs$patient)) {
  patient_df <- subset(df_obs, patient == p)
  cross_table <- table(patient_df$time_point, patient_df$celltype_r2)
  ROIE_list[[p]] <- data.frame(chisq.test(cross_table)$observed/chisq.test(cross_table)$expected)
  names(ROIE_list[[p]])[1] <- 'time'
  ROIE_list[[p]]$patient <- p
  ROIE_list[[p]]$dataset <- paste0(str_split(ROIE_list[[p]]$patient, '_', simplify = T)[,1],'_',str_split(ROIE_list[[p]]$patient, '_', simplify = T)[,2])
  ROIE_list[[p]]$time <- factor(ROIE_list[[p]]$time, levels = c('Pre','On/Post'))
  names(ROIE_list[[p]])[2:3] <- c('celltype','ratio')
}
ROIE_mat <- do.call(rbind, ROIE_list)
ROIE_mat$ratio <- sapply(ROIE_mat$ratio, function(x) log2(x + 0.0001))
ROIE_mat$response <- df_obs$response[match(ROIE_mat$patient, df_obs$patient)]
ROIE_mat$treatment <- df_obs$treatment[match(ROIE_mat$patient, df_obs$patient)]
ROIE_mat$modality <- df_obs$modality[match(ROIE_mat$patient, df_obs$patient)]
ROIE_mat$subset <- paste0(ROIE_mat$dataset, '(', ROIE_mat$treatment, ')')

category <- unique(ROIE_mat$subset)
plots <- list()
df_list <- list()
for (name in category){
  df <- filter(ROIE_mat, subset == name)
  df_check <- df_obs %>% filter(subset == name) %>% 
    group_by(celltype_r2, patient, time_point) %>% 
    summarise(n()) %>% 
    pivot_wider(names_from = time_point, values_from = `n()`)
  df_check[is.na(df_check)] <- 0
  # filter celltype
  celltype_list <- c()
  for (i in unique(df_check$celltype_r2)){
    if(nrow(filter(df_check[df_check$celltype_r2 == i,], Pre>=3, `On/Post`>=3)) >= 3)
      celltype_list <- c(celltype_list, i)
  }
  # make plots for each celltype
  plot_list <- list()
  subtype_sig <- c()
  t_stat <- c()
  t_p <- c()
  n=0
  for (subtype in celltype_list){
    p_filt <- filter(df_check[df_check$celltype_r2 == subtype,], !(Pre<3 & `On/Post`<3))$patient
    d <- filter(df, celltype == subtype, patient %in% p_filt)# filter out patients without sufficient cell counts
    d <- d[order(d$time, decreasing = T),] 
    d_wide <- pivot_wider(d, names_from = time, values_from = ratio)
    t_test <- t.test(d_wide$`On/Post`, d_wide$Pre, paired = T, alternative = 'two.sided')
    t_stat <- c(t_stat, t_test$statistic)
    t_p <- c(t_p, t_test$p.value)
    if (t_test$p.value < 0.05){
      n = n+1
      print(subtype)
      plot_list[[n]] <-
        ggpaired(d, x = "time", y = "ratio", fill = "white", line.color = "gray", line.size = 0.3,
                 palette = "npg", group = 'patient', title = subtype) +
        geom_hline(yintercept = 0, color = 'black', linetype = 'dashed', size = 0.3) +
        # geom_boxplot(fill = 'white') +
        geom_point(aes(color = response)) +
        stat_compare_means(paired = TRUE, method = 't.test', label.y = 2) + 
        theme(legend.direction = "vertical") +
        xlab("Time Point") + ylab("log2(R(O/E)+0.001)")
      celltype_sig <- c(subtype_sig, subtype)
    }
    t_fdr <- p.adjust(t_p, method = 'fdr')
    df$t_stat <- t_stat[match(df$celltype, celltype_list)]
    df$fdr <- t_fdr[match(df$celltype, celltype_list)]
    df_list[[which(category==name)]] <- df %>% distinct(patient, celltype, .keep_all = T)
  }
  plots[[which(category %in% name)]] <- wrap_plots(plot_list) + plot_layout(guides = 'collect') + plot_annotation(title = name)
}

df_list <- lapply(df_list[-7], function(x) {
  print(colnames(x))
  x <- x %>%
    drop_na(t_stat) %>%
    distinct(celltype, .keep_all = TRUE) %>%
    pivot_wider(names_from = subset, values_from = t_stat) %>%
    .[,c(2,10)]
  x$celltype <- as.character(x$celltype)
  return(x)
})
df_heatmap <- Reduce(full_join, df_list) %>% 
  arrange(celltype) %>% 
  column_to_rownames(var = 'celltype') 
Heatmap(df_heatmap,name = "mat", na_col = "white", cluster_rows = F, heatmap_legend_param = list(
  title = "T-score"
))

# foldchange visualization
df <- read.csv('/bigdata/zlin/Melanoma_meta/tables/subtype/myeloids_SingleR.csv')
df_wide <- read.csv('/bigdata/zlin/Melanoma_meta/tables/subtype/myeloids_SingleR_wide.csv', check.names = FALSE)
df_wide$time_point <- factor(df_wide$time_point, levels = c('Pre','On'))
df_fc <- read.csv('/bigdata/zlin/Melanoma_meta/tables/subtype/myeloids_SingleR_fc.csv', check.names = FALSE)
df_fc[,2:(ncol(df_fc)-3)] <- log2(df_fc[,2:(ncol(df_fc)-3)])
annotation <- 'SingleR'
cell_types <- df$Cell_type[df$log10Pvalue > -log10(alpha)]

df <- read.csv('/bigdata/zlin/Melanoma_meta/tables/subtype/myeloids_MetaTiME.csv')
df_wide <- read.csv('/bigdata/zlin/Melanoma_meta/tables/subtype/myeloids_MetaTiME_wide.csv', check.names = FALSE)
df_wide$time_point <- factor(df_wide$time_point, levels = c('Pre','On'))
df_fc <- read.csv('/bigdata/zlin/Melanoma_meta/tables/subtype/myeloids_MetaTiME_fc.csv', check.names = FALSE)
df_fc[,2:(ncol(df_fc)-3)] <- log2(df_fc[,2:(ncol(df_fc)-3)])
annotation <- 'MetaTiME'
cell_types <- df$Cell_type[df$log10Pvalue > -log10(alpha)]

# Set significance threshold (adjust as needed)
alpha <- 0.05

# Create a volcano plot
ggplot(df, aes(x = log2FC, y = log10Pvalue)) +
  geom_point(color = 'gray', alpha = 0.5, size = 3) +
  geom_point(data = subset(df, log10Pvalue > -log10(alpha)), color = 'red', alpha = 0.5, size = 3) +
  geom_hline(yintercept = -log10(alpha), color = 'black', linetype = 'dashed', size = 0.3) +
  geom_vline(xintercept = 0, color = 'black', linetype = 'dashed', size = 0.3) +
  labs(x = 'Comb vs Mono (log2FC)', y = '-log10(p-value)', title = 'Myeloid cells (Mann-Whitney U test)') +
  theme_minimal() +
  theme(axis.line = element_line(color = 'black', size = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  geom_text_repel(data = subset(df, log10Pvalue > -log10(alpha)), aes(label = Cell_type), size = 4, color = 'navy')
ggsave(filename = paste0('/bigdata/zlin/Melanoma_meta/figures/Myeloids/vol_Myeloids_',annotation,'_early_treat.pdf'),width = 5, height = 5)

plots <- list()
# Create a separate plot for each cell type
for (cell_type in cell_types) {
  i <- which(cell_types %in% cell_type)
  # Subset the data for the specific cell type
  subset_df <- df_wide[, c('time_point', 'patient', 'dataset', 'modality', cell_type)]
  colnames(subset_df)[5] <- 'proportion'
  colnames(subset_df)[3] <- 'Source'
  # Create the plot
  plots[[i]] <- ggplot(subset_df, aes(x = time_point, y = proportion, group = patient, color = Source)) +
    geom_point(shape = 16, size = 2) +
    geom_line(color = 'black', size = 0.2) +
    theme_classic() +
    labs(y = "Proportion(%)", x = 'Time point') +
    facet_wrap(vars(modality)) +
    scale_color_nejm() +
    ggtitle(cell_type)
}
wrap_plots(plots) + plot_layout(guides = 'collect')
ggsave(filename = paste0('/bigdata/zlin/Melanoma_meta/figures/Myeloids/',annotation,'_early4_treat.pdf'),width = 10, height = 6)

plots <- list()
# Create a separate plot for each cell type
for (cell_type in cell_types) {
  i <- which(cell_types %in% cell_type)
  # Subset the data for the specific cell type
  subset_df <- df_fc[, c('patient', 'modality', cell_type)]
  colnames(subset_df)[3] <- 'lfc'
  # Create the plot
  plots[[i]] <- ggplot(subset_df, aes(x = modality, y = lfc, fill = modality)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, size=1) +
    geom_hline(yintercept = 0, color = 'black', linetype = 'dashed', size = 0.3) +
    theme_classic() +
    labs(y = "Foldchange (log2)", x = 'Combi vs Mono') +
    scale_fill_nejm() +
    ggtitle(cell_type) + stat_compare_means(method = "wilcox.test")
}
wrap_plots(plots) + plot_layout(guides = 'collect')
ggsave(filename = paste0('/bigdata/zlin/Melanoma_meta/figures/Myeloids/',annotation,'_early4_treat_fc.pdf'),width = 8, height = 6)


filter_patients <- function(adata) {
  # Create the crosstab
  ct <- table(adata$patient, adata$time_point)
  filtered <- rownames(ct)[ct[, 1] < 30 | ct[, 2] < 30]
  return(filtered)
}

# Myeloids
mye_list <- readRDS('/bigdata/zlin/Melanoma_meta/data/Myeloids/sub_list.rds')
datasets_names <- c('SKCM_Becker','HNSC_Franken','BRCA_Bassez1','BRCA_Bassez2','BCC_Yost','CRC_Li','TNBC_Li')
mye_list[[1]]$time_point <- as.character(mye_list[[1]]$time_point)
mye_list[[1]]$time_point[mye_list[[1]]$time_point == 'Before'] <- 'Pre'
mye_list[[2]]$patient <- as.character(mye_list[[2]]$patient)

for (i in 1:length(mye_list)) {
  print(paste0(datasets_names[i], ":"))
  print(paste0(toString(filter_patients(mye_list[[i]]))))
}

# [1] "SKCM_Becker:"
# [1] "P10, P2, P4, P5, P6, P8"
# [1] "HNSC_Franken:"
# [1] "16, 20"
# [1] "BRCA_Bassez1:"
# [1] "BIOKEY_20, BIOKEY_23, BIOKEY_25, BIOKEY_29"
# [1] "BRCA_Bassez2:"
# [1] "BIOKEY_32, BIOKEY_36, BIOKEY_38, BIOKEY_39"
# [1] "BCC_Yost:"
# [1] "su007, su009, su012"
# [1] "CRC_Li:"
# [1] "P28, P30"
# [1] "TNBC_Li:"
# [1] ""

early <- mye_list[1:4]
early_names <- datasets_names[1:4]

select_col <- c('response', 'treatment', 'MetaTiME_overcluster', 'celltype_sc_r2', 'batch', 'dataset', 'time_point', 'patient')
df_early <- rbind(early[[1]]@meta.data[,select_col], early[[2]]@meta.data[,select_col], early[[3]]@meta.data[,select_col], early[[4]]@meta.data[,select_col]) 

# Removing cells identified as other lineages
substrings <- c('T', 'NK', 'Plasma', 'B', 'Fibroblast', 'nan','Endothelial','Erythrocytes','Myofibroblast')
excluded_substring <- 'Macrophage-IL1-NFKB'
df_early <- df_early %>%
  filter(!grepl(paste(substrings, collapse = '|'), MetaTiME_overcluster) & !grepl(excluded_substring, MetaTiME_overcluster))

# Harmonize treatment labels
df_early$modality <- df_early$treatment
df_early$modality <- recode(df_early$modality, 'PD1+CTLA4' = 'Comb', 'PDL1-CTLA4' = 'Comb', 'PDL1' = 'Mono', 'PD1' = 'Mono', 'anti-PD1' = 'Mono', 'Chemo+PD1' = 'Mono')

# Add new columns to df_early
df_early$count_singler <- with(df_early, ave(celltype_sc_r2, batch, celltype_sc_r2, FUN = length))
df_early$count_metatime <- with(df_early, ave(MetaTiME_overcluster, batch, MetaTiME_overcluster, FUN = length))
df_early$sample_total <- with(df_early, ave(celltype_sc_r2, batch, FUN = length))
df_early$celltype_proportion <- as.numeric(df_early$count_singler) / (as.numeric(df_early$sample_total))
pivot_df <- df_early %>%
  select(batch, celltype_sc_r2, celltype_proportion) %>%
  pivot_wider(names_from = celltype_sc_r2, values_from = celltype_proportion)

df_fc$response[df_fc$response == 'NE' & df_fc$dataset =='SKCM_Becker'| df_fc$response == 'n/a'] <- 'NE(unknown)'
df_fc$response[df_fc$response %in% c('PD','SD','NE')] <- 'Non-response'
df_fc$response[df_fc$response %in% c('CR','PR','E')] <- 'Reseponse'



df_fc <- read.csv('/bigdata/zlin/Melanoma_meta/tables/subtype/fc.csv')

library("RColorBrewer")

tr_vec <- dittoColors()[1:length(unique(df_fc$treatment))]
names(tr_vec) <- unique(df_fc$treatment)
t_vec <- dittoColors()[9:12]
names(t_vec) <- unique(df_fc$interval_cat)
r_vec <- pal_npg("nrc")(3)
names(r_vec) <- unique(df_fc$response_cat)
d_vec <- dittoColors()[3:(2+length(unique(df_fc$dataset)))]
names(d_vec) <- unique(df_fc$dataset)
ha_col = HeatmapAnnotation(
  Response = df_fc$response_cat,
  Dataset = df_fc$dataset,
  Interval = df_fc$interval_cat,
  Treatment = df_fc$treatment,
  col = list(Response = r_vec, Dataset = d_vec, Interval = t_vec, Treatment = tr_vec)
)
Heatmap(t(log2(df_fc[, 2:(which(names(df_fc)=='dataset')-1)])), 
        show_column_names = FALSE, cluster_rows = T,
        top_annotation = ha_col, 
        col = rev(brewer.pal(n=10,name = 'RdBu')),
        heatmap_legend_param = list(title = 'log2(FC)'))

df_obs <- read.csv('/bigdata/zlin/Melanoma_meta/tables/subtype/obs.csv') %>% filter(dataset =='SKCM_Becker', patient =='P1')

N.o.byPatient <- table(df_obs$celltype_sc_r2, df_obs$time_point)
R.oe <- apply(N.o.byPatient,1,function(x){
  res.chisq <- chisq.test(x)
  return((res.chisq$observed)/(res.chisq$expected))})
R.oe


tr1_vec <- dittoColors()[1:length(unique(df_obs$treatment))]
names(tr1_vec) <- unique(df_obs$treatment)
tr2_vec <- dittoColors()[8:(7+length(unique(df_obs$treatment_cat)))]
names(tr2_vec) <- unique(df_obs$treatment_cat)
t_vec <- pal_npg("nrc")(length(unique(df_obs$interval_cat)))
names(t_vec) <- unique(df_obs$interval_cat)
type_vec <- pal_d3()(length(unique(df_obs$cancer_type)))
names(type_vec) <- unique(df_obs$cancer_type)
d_vec <- pal_nejm()(length(unique(df_obs$dataset)))
names(d_vec) <- unique(df_obs$dataset)
r_vec <- pal_jco()(length(unique(df_obs$response)))
names(r_vec) <- unique(df_obs$response)
r2_vec <- pal_aaas()(length(unique(df_obs$response_cat)))
names(r2_vec) <- unique(df_obs$response_cat)
ha_col = HeatmapAnnotation(
  Treatment = df_obs$treatment[match(colnames(ROIE_mat), df_obs$patient)],
  Modality = df_obs$treatment_cat[match(colnames(ROIE_mat), df_obs$patient)],
  Dataset = df_obs$dataset[match(colnames(ROIE_mat), df_obs$patient)],
  Interval = df_obs$interval_cat[match(colnames(ROIE_mat), df_obs$patient)],
  `Cancer type` = df_obs$cancer_type[match(colnames(ROIE_mat), df_obs$patient)],
  Response = df_obs$response[match(colnames(ROIE_mat), df_obs$patient)],
  `Response harmonized` = df_obs$response_cat[match(colnames(ROIE_mat), df_obs$patient)],
  col = list(Treatment = tr1_vec, Modality = tr2_vec, Interval = t_vec, Dataset = d_vec, `Cancer type`=type_vec, Response = r_vec,`Response harmonized` = r2_vec)
)
Heatmap(ROIE_mat, 
        show_column_names = FALSE, cluster_rows = T,
        top_annotation = ha_col, 
        col = rev(brewer.pal(n=10,name = 'RdBu')),
        heatmap_legend_param = list(title = 'R o/e'))

# calculate Spearman correlation
mat.cor <- cor(t(ROIE_mat), method = "spearman")
# calculate distance based on correlation matrix
mat.dist <- as.dist(1 - abs(mat.cor))
# hierarchical clustering
mat.hc <- hclust(mat.dist, method = "ward.D2")
# dendrogram
mat.dend <- as.dendrogram(mat.hc)
Heatmap(ROIE_mat, 
        show_column_names = FALSE, cluster_columns = mat.hc, cluster_rows = T,
        top_annotation = ha_col, 
        col = rev(brewer.pal(n=10,name = 'RdBu')),
        heatmap_legend_param = list(title = 'R o/e'))


for (name in category){
  df <- filter(ROIE_mat, eval(parse(text = groupby)) == name)
  df_check <- df_obs %>% filter(eval(parse(text = groupby)) == name) %>% 
    group_by(celltype_r2_comb, patient, time_point) %>% 
    summarise(n()) %>% 
    pivot_wider(names_from = time_point, values_from = `n()`)
  df_check[is.na(df_check)] <- 0
  # filter celltype
  celltype_list <- c()
  for (i in unique(df_check$celltype_r2_comb)){
    if(nrow(filter(df_check[df_check$celltype_r2_comb == i,], Pre>=3, `On/Post`>=3)) >= 3)
      celltype_list <- c(celltype_list, i)
  }
  # make plots for each celltype
  plot_list <- list()
  subtype_sig <- c()
  t_stat <- c()
  t_p <- c()
  n=0
  for (subtype in celltype_list){
    p_filt <- filter(df_check[df_check$celltype_r2_comb == subtype,], !(Pre<3 & `On/Post`<3))$patient
    d <- filter(df, celltype == subtype, patient %in% p_filt)# filter out patients without sufficient cell counts
    d <- d[order(d$time, decreasing = T),] 
    d_wide <- pivot_wider(d, names_from = time, values_from = ratio)
    t_test <- t.test(d_wide$`On/Post`, d_wide$Pre, paired = T, alternative = 'two.sided')
    t_stat <- c(t_stat, t_test$statistic)
    t_p <- c(t_p, t_test$p.value)
    if (t_test$p.value < 0.05){
      n = n+1
      print(subtype)
      plot_list[[n]] <-
        ggpaired(d, x = "time", y = "ratio", fill = "white", line.color = "gray", line.size = 0.3,
                 palette = "npg", group = 'patient', title = subtype) +
        geom_hline(yintercept = 0, color = 'black', linetype = 'dashed', size = 0.3) +
        geom_point(aes(color = cancer_type)) +
        stat_compare_means(paired = TRUE, method = 't.test', label.y = 2) + 
        theme(legend.direction = "vertical", legend.position = 'right') +
        xlab("Time Point") + ylab("log2(Ro/e)")
      celltype_sig <- c(subtype_sig, subtype)
    }
    t_fdr <- p.adjust(t_p, method = 'fdr')
    df$t_stat <- t_stat[match(df$celltype, celltype_list)]
    df$pvalue <- t_p[match(df$celltype, celltype_list)]
    df$fdr <- t_fdr[match(df$celltype, celltype_list)]
    df_list[[which(category==name)]] <- df %>% distinct(patient, celltype, .keep_all = T)
  }
  plots[[which(category %in% name)]] <- wrap_plots(plot_list) + plot_layout(guides = 'collect') + plot_annotation(title = name)
}

# Tell if TCR is pre-existing
# Step 1: Extract the metadata
meta_data <- seu@meta.data

# Step 2: Group, summarize, and determine the status
cdr3s_status <- meta_data %>%
  group_by(patient, cdr3s_nt_unique) %>%
  summarize(has_pre = "pre" %in% time_point,
            has_post = "post" %in% time_point) %>%
  mutate(status = ifelse(has_pre & has_post, "pre-existing", "non-preexisting")) %>%
  select(patient, cdr3s_nt_unique, status)

# Step 3: Merge with the original metadata
new_meta_data <- left_join(meta_data, cdr3s_status, by = c("patient", "cdr3s_nt_unique"))

# Step 4: Update the Seurat object's metadata
seu@meta.data <- new_meta_data

seu$patient <- str_replace(seu$patient, 'BIOKEY', 'BRCA_Bassez2')
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='timepoint')] <- "time_point"
seu$time_point <- ifelse(seu$time_point == 'On', 'Post', 'Pre')
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='expansion')] <- "response"
colnames(seu@meta.data)[which(colnames(seu@meta.data)=='BC_type')] <- "subtype"
seu$sample <- paste0(seu$patient,'_', seu$time_point)
seu$treatment <- 'aPD1'
seu$dataset <- 'BRCA_Bassez2'
seu$cancertype <- 'BRCA'
seu$response <- ifelse(seu$response == 'E', 'RE', 'NR')
seu$res_metric <- 'T-cell expansion'

# Loading dataset
SKCM_Becker <- qread('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/seu_r1.qs')
BRCA_Bassez1 <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r1.qs')
BRCA_Bassez2 <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez2/seu_r1.qs')
TNBC_Zhang <- qread('/bigdata/zlin/Melanoma_meta/data/TNBC_Zhang/seu_r1.qs')
BCC_Yost <- qread('/bigdata/zlin/Melanoma_meta/data/BCC_Yost/seu_r1.qs')
SCC_Yost <- qread('/bigdata/zlin/Melanoma_meta/data/SCC_Yost/seu_r1.qs')
HNSC_IMCISION <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_IMCISION/seu_r1.qs')
HNSC_Luoma <- qread('/bigdata/zlin/Melanoma_meta/data/HNSC_Luoma/seu_r1.qs')
NSCLC_Liu <- qread('/bigdata/zlin/Melanoma_meta/data/NSCLC_Liu/seu_r1.qs')
CRC_Li <- qread('/bigdata/zlin/Melanoma_meta/data/CRC_Li/seu_r1.qs')
PCa_Hawley <- qread('/bigdata/zlin/Melanoma_meta/data/PCa_Hawley/seu_r1.qs')

pkgs <- c('Seurat','SingleCellExperiment','magrittr','qs')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
datasets <- c('SKCM_Becker', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'BCC_Yost', 'SCC_Yost', 'HNSC_IMCISION', 'HNSC_Luoma', 'NSCLC_Liu', 'CRC_Li','PCa_Hawley')
lapply(datasets, function(dataset){
  print(dataset)
  sce <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) %>% as.SingleCellExperiment()
  dior::write_h5(sce, file = paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.h5'), object.type = 'singlecellexperiment')
  })

hvg <- VariableFeatures(seuratObj)
gene.pattern <- c("MALAT1", "^MT-", "^RPL", "^RPS", "^LOC(0-9)", "^TR(A|B|G|D)V", "^MTRNR")
hvg <- hvg[!hvg %in% grep(paste0(gene.pattern, collapse = "|"), hvg, value = T)]

# Mice

stress_genes= c("G0s2", "Jun", "Junb", "Jund", "Fos", "Dusp1", "Cdkn1a", "Fosb", "Btg2", "Klf6", "Klf4")
cc_genes =  unique(c(paste0(toupper(substr(tolower(as.character(unlist(cc.genes.updated.2019))), 1, 1)), substr(tolower(as.character(unlist(cc.genes.updated.2019))), 2, nchar(tolower(as.character(unlist(cc.genes.updated.2019)))))))) #, "1700020L24Rik" ,"5730416F02Rik" ,"Agpat4","Asf1b", "Aspm","Ccdc18","Ccr9","Clspn","Cyp4f17","Dek","Dnmt3b" ,"Dtl","Fancm","Fignl1","Gm14214","Gm14730","Gm15428","Gm15448","Gm21992","Gm23248","Gm26465","Gm5145" ,"Mcm4","Mcm8","Mki67","Oip5","Pcna","Pcna-ps2","Pigv","Rnd1","Snrpa","Ube2c"
IFN_genes = unique(c(grep("Irf", rownames(seurat@assays$RNA@counts), v = T), grep("Ifi", rownames(rownames(seurat@assays$RNA@counts)), v = T), "Cd2","Cd3d", "Cmpk2","Cxcl10", "Isg15","Isg20","Oasl1" ,"Phf11b" ,"Plac8", "Rsad2","Rtp4","Sdf4", "Slfn4","Tnfrsf18" ,"Tnfrsf4","Tnfrsf9","Trex1","Usp18","Xaf1"))
ccl_genes = grep("Ccl", rownames(seurat@assays$RNA@counts), v = T)
MHC_genes =  c(grep("^H2-", rownames(seurat@assays$RNA@counts), v = T), "Cd74", "B2m")
hist_genes = grep("Hist", rownames(seurat@assays$RNA@counts), v = T)
comp_genes =  c(grep("^C1q", rownames(seurat@assays$RNA@counts), v = T))
ig_genes = c(grep("^Igj|^Igh|^Igk|^Igl", rownames(seurat@assays$RNA@counts), v = T))
hb_genes = c(grep("^Hba|^Hbb", rownames(seurat@assays$RNA@counts), v = T))
cc_gene_module = get(load("cc_gene_module.RData"))
bad_features = unique(c(hist_genes, cc_genes, stress_genes, "Malat1", "Neat1", "Xist", "Actb", grep("^mt-|^Mtmr|^Rps|^Rpl|Rik|AL|-rs|-ps|Mir|Atp|Gm|Uqc", rownames(seurat@assays$RNA@counts), v = T)))

seurat = RunPCA(seurat, assay = assay, features = setdiff(seurat@assays[[assay]]@var.features, bad_features))

# Set palette
if (palette == "npg") {
  colors <- colorRampPalette((pal_npg("nrc")(9)))(length(unique(plot_data$seurat_by))) #npg
} else if (palette == "jama") {
  colors <- colorRampPalette((pal_jama("default")(7)))(length(unique(plot_data$seurat_by))) #jama
} else if (palette == "jco") {
  colors <- colorRampPalette((pal_jco("default")(10)))(length(unique(plot_data$seurat_by))) #jco
} else if (palette == "lancet") {
  colors <- colorRampPalette((pal_lancet("default")(8)))(length(unique(plot_data$seurat_by))) #lacent
} else if (palette == "nejm"){
  colors <- colorRampPalette((pal_nejm("default")(8)))(length(unique(plot_data$seurat_by))) #nejm
}

df_check <- meta_combi |>
  select(celltype_r2, patient, time_point, interval, response, treatment, count_r2, res_metric, cancertype, prior) |>
  distinct(celltype_r2, patient, time_point, .keep_all = T) |>
  pivot_wider(names_from = time_point, values_from = count_r2, values_fill = 0) |>
  filter(abs(Pre-Post) >= 3)
subtypes <-
  df_check |>
  group_by(celltype_r2) |>
  dplyr::summarise(count = n()) |>
  filter(count >= n.sample) |>
  pull(celltype_r2)
freq_mat <- distinct(meta_combi, sample, celltype_r2, .keep_all = T)
freq_mat <- meta_combi |>
  distinct(celltype_r2, sample, .keep_all = T) |>
  select(freq_r2_comp, dataset, response, modality, int_cat, patient, time_point, celltype_r2) |>
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  pivot_longer(cols = c('Pre', 'Post'), names_to = 'time_point', values_to = 'freq_r2_comp')

library(msigdbr)
pathwaysDF <- msigdbr("mouse", category="H")
pathways <- split(as.character(pathwaysDF$entrez_gene), pathwaysDF$gs_name)
res <- runGSEA(head(marker_cosg$names$PC_IGHA, n=100) , universe=rownames(seu), category = c('C2'))


##############################################################################
# 0. Load Necessary Packages
##############################################################################
library(Seurat)
library(NMF)        # For NMF analysis
library(dplyr)      # Data manipulation
library(tidyr)      # Data tidying
library(ComplexHeatmap)
library(ggplot2)

##############################################################################
# 1. Load the Pre-built Seurat Object and Extract Malignant Cells
##############################################################################
# Assuming you already have a Seurat object "Seurat_obj" containing all cell information
# The "meta.data" slot contains a column "Harmony_SNN_res.0.6" for cell type annotations
Seurat_obj <- readRDS("seurat.rds")

# Extract the target cell names
Target_cells <- rownames(Seurat_obj@meta.data)[Seurat_obj@meta.data$Harmony_SNN_res.0.6 %in% c("9", "12", "14", "16", "21", "22")]

# Create a new Seurat object containing only the target cells
Target_seurat <- subset(Seurat_obj, cells = Target_cells)

# Ensure RNA expression data is used for subsequent operations
DefaultAssay(Target_seurat) <- "RNA"

# Optionally, check the number of target cells per sample
table(Target_seurat@meta.data$sample)

##############################################################################
# 2. Perform NMF for Each Sample and Extract the Top 30 Genes per Program
##############################################################################

# Set the rank for decomposition (12 programs)
rank_k <- 12

# List to store the top 30 genes for all samples
# Structure: Each list element corresponds to "sample_Program" and contains top 30 genes
all_programs_top30 <- list()

# Get the unique list of samples
sample_list <- unique(Target_seurat@meta.data$sample)

# Define a function to perform NMF on a single sample and extract top 30 genes
run_nmf_extract_top30 <- function(sample_id) {
  
  # 1) Extract target cells for the given sample
  temp_cells <- rownames(Target_seurat@meta.data)[
    Target_seurat@meta.data$sample == sample_id
  ]
  
  # 2) Subset the data for the sample
  temp_seurat <- subset(Target_seurat, cells = temp_cells)
  
  # 3) Standard preprocessing (Normalize, FindVariableFeatures, ScaleData)
  temp_seurat <- NormalizeData(temp_seurat)
  temp_seurat <- FindVariableFeatures(temp_seurat)
  temp_seurat <- ScaleData(temp_seurat, do.center = FALSE)  # Do not center data
  
  # 4) Extract normalized or scaled expression matrix
  exp_matrix <- temp_seurat@assays$RNA@scale.data
  
  # 5) Perform NMF
  #    Use method="snmf/r", set seed = 10, rank = rank_k (12)
  res_nmf <- nmf(exp_matrix, rank_k, seed = 10, method = "snmf/r")
  
  # 6) Extract the "Basis matrix" (gene weights)
  W <- basis(res_nmf)  # Dimensions: genes x rank_k
  
  # 7) Sort genes by weight for each program and extract the top 30
  top30_genes_list <- apply(W, 2, function(x) {
    # Sort weights in descending order
    top_idx <- order(x, decreasing = TRUE)[1:30]
    rownames(W)[top_idx]
  })
  
  # Name columns as sample_Program1, sample_Program2, ..., sample_Program12
  colnames(top30_genes_list) <- paste0(sample_id, "_Program", 1:rank_k)
  
  # Return the top 30 genes for each program
  return(top30_genes_list)
}

# Apply the function to all samples
for(sample_id in sample_list){
  cat("Running NMF for sample:", sample_id, "\n")
  gene_mat <- run_nmf_extract_top30(sample_id)
  
  # Add each program's top 30 genes to the global list
  for(program_name in colnames(gene_mat)){
    all_programs_top30[[program_name]] <- gene_mat[, program_name]
  }
}

##############################################################################
# 3. Calculate Jaccard Similarity Between Programs and Perform Clustering
##############################################################################

# Extract the names of all programs
program_names <- names(all_programs_top30)

# Create an empty matrix for storing Jaccard similarities
n_prog <- length(program_names)
jaccard_mat <- matrix(0, nrow = n_prog, ncol = n_prog,
                      dimnames = list(program_names, program_names))

# Compute pairwise Jaccard similarity
for(i in seq_len(n_prog)){
  for(j in seq_len(n_prog)){
    set1 <- all_programs_top30[[i]]
    set2 <- all_programs_top30[[j]]
    inter_len <- length(intersect(set1, set2))
    union_len <- length(union(set1, set2))
    jaccard_mat[i, j] <- inter_len / union_len
  }
}

# Perform hierarchical clustering on the Jaccard similarity matrix
dist_mat <- as.dist(1 - jaccard_mat)  # Convert similarity to distance
hc <- hclust(dist_mat, method = "ward.D2")

##############################################################################
# 4. Define Meta-Programs Based on Clustering Results and Annotate
##############################################################################
k_meta <- 5
cluster_cut <- cutree(hc, k = k_meta)

# Create a dataframe to store program clusters
program_cluster_df <- data.frame(
  Program = program_names,
  Cluster = cluster_cut
)

# Assign names to meta-programs (based on functional enrichment)
meta_program_names <- c("MetaProg_Inflammation",
                        "MetaProg_CellCycle",
                        "MetaProg_Differentiation",
                        "MetaProg_StressResponse",
                        "MetaProg_Unknown")

program_cluster_df$MetaProgramName <- meta_program_names[program_cluster_df$Cluster]

##############################################################################
# 5. Extract Top 30 Genes for Each Meta-Program
##############################################################################

meta_program_top30 <- list()

for(meta in unique(program_cluster_df$MetaProgramName)){
  progs <- program_cluster_df$Program[program_cluster_df$MetaProgramName == meta]
  genes <- unlist(all_programs_top30[progs])
  gene_freq <- table(genes)
  gene_freq_sorted <- sort(gene_freq, decreasing = TRUE)
  top30 <- names(gene_freq_sorted)[1:30]
  meta_program_top30[[meta]] <- top30
}

##############################################################################
# 6. Optional: Save Results
##############################################################################
# Save results for later analysis or visualization
write.csv(meta_program_top30, file = "MetaPrograms_Top30_Genes.csv", row.names = FALSE)

# check sample quality
df <- data.frame(
  nfeature = sapply(seu_list, function(seu){nfeature <- sum(Matrix::rowSums(seu@assays$RNA$counts>0)>3)
  return(nfeature)}),
  ncell = sapply(seu_list, function(seu){ncell <- ncol(seu)
  return(ncell)})
)
# remove samples with <200 cells or <9000 genes 
sample_rm <- rownames(filter(df, nfeature < 9000 | ncell < 200)); sample_rm





