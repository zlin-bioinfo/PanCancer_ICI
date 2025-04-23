library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(janitor)
library(MultiAssayExperiment)
genelist <- read.csv('tables/cosg_genes.csv', check.names = F) |> as.list()
genelist$`Melanocytes(CNA+)` <- NULL
genelist$`Epithelial(CNA+)` <- NULL
genelist$Cycling <- NULL
# Melanoma 
# GSE91061
expr_mat <- read.csv('data/bulk_datasets/Melanoma/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv')
# convert entrezid to gene symbol
e2s <- bitr(expr_mat$X, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
expr_mat <- merge(e2s, expr_mat, by.x = 'ENTREZID', by.y = 'X')
expr_mat <- expr_mat[!duplicated(expr_mat$SYMBOL), ]
expr_mat <- expr_mat |> tibble::column_to_rownames(var = 'SYMBOL') |> subset(select = -ENTREZID)
expr_mat <- expr_mat[,!str_detect(colnames(expr_mat), 'Pt109_On')]
colnames(expr_mat) <- paste0(str_split(colnames(expr_mat), '_', simplify = T)[,1], '_', str_split(colnames(expr_mat), '_', simplify = T)[,2])
expr_mat_pre <- expr_mat[,str_detect(colnames(expr_mat), '_Pre')] |> log1p()
colnames(expr_mat_pre) <- str_replace(colnames(expr_mat_pre), '_Pre', '')
ssgseaPar <- gsvaParam(as.matrix(expr_mat_pre), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- readxl::read_xlsx('data/bulk_datasets/Melanoma/GSE91061/1-s2.0-S0092867417311224-mmc2.xlsx', skip = 1)
clin_info <- clin_info |> 
  data.frame(check.names = F) |> 
  mutate(time_os = `Time to Death\r\n(weeks)`/4,
         os = ifelse(`Dead/Alive\r\n(Dead = True)` == TRUE, 1, 0),
         time_pfs = NA,
         pfs = NA,
         response = Response,
         patient = Patient,
         sample = patient,
         tx = 'aPD1',
         tx_status = 'Baseline',
         cancertype = 'Melanoma',
         cohort = 'Riaz') |> 
  dplyr::select(patient, sample, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |> 
  filter(patient %in% colnames(score))
comb_pre <- cbind(clin_info, t(score)[clin_info$patient,])
# write.csv(comb_pre, 'data/bulk_datasets/Melanoma/GSE91061/comb_pre.csv', row.names = F)

expr_mat_on <- expr_mat[,str_detect(colnames(expr_mat), '_On')] |> log1p()
colnames(expr_mat_on) <- str_replace(colnames(expr_mat_on), '_On', '')
ssgseaPar <- gsvaParam(as.matrix(expr_mat_on), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- readxl::read_xlsx('data/bulk_datasets/Melanoma/GSE91061/1-s2.0-S0092867417311224-mmc2.xlsx', skip = 1)
clin_info <- clin_info |> 
  data.frame(check.names = F) |> 
  mutate(time_os = `Time to Death\r\n(weeks)`/4,
         os = ifelse(`Dead/Alive\r\n(Dead = True)` == TRUE, 1, 0),
         time_pfs = NA,
         pfs = NA,
         response = Response,
         patient = Patient,
         sample = patient,
         tx = 'aPD1',
         tx_status = 'Treated',
         cancertype = 'Melanoma',
         cohort = 'Riaz') |> 
  dplyr::select(patient, sample, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |> 
  filter(patient %in% colnames(expr_mat_on))
comb_on <- cbind(clin_info, t(score)[clin_info$patient,])
# write.csv(comb_on, 'data/bulk_datasets/Melanoma/GSE91061/comb_on.csv', row.names = F)
comb <- rbind(comb_pre, comb_on)
write.csv(comb, 'data/bulk_datasets/Melanoma/GSE91061/comb.csv', row.names = F)

# Nathanson (27956380)
expr_mat <- readRDS('data/bulk_datasets/TIGER/Melanoma-Nathanson_2017/Melanoma-Nathanson_2017.Response.Rds') |> 
  column_to_rownames(var = 'GENE_SYMBOL') |> 
  drop_na() 
detected_counts <- rowSums(expr_mat > 0)
expr_mat <- expr_mat[detected_counts >= 24, ]
expr_mat <- log1p(expr_mat)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- read.csv('data/bulk_datasets/TIGER/Melanoma-Nathanson_2017/Melanoma-Nathanson_2017.Response.tsv', sep = '\t')
clin_info <- clin_info |> 
  mutate(time_os = overall.survival..days./30,
         os = ifelse(vital.status == 'Dead', 1, 0),
         patient = patient_name,
         sample = as.character(sample_id),
         tx = 'aCTLA4',
         tx_status = ifelse(Treatment == 'PRE', 'Baseline', 'Treated'),
         cancertype = 'Melanoma',
         cohort = 'Nathanson') |> 
  select(patient, sample, response, os, time_os, tx, tx_status, cancertype, cohort) |> 
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/TIGER/Melanoma-Nathanson_2017/comb.csv', row.names = F)

# Van Allen (phs000452.v2.p1)
expr_mat <- readRDS('data/bulk_datasets/TIGER/Melanoma-phs000452/Melanoma-phs000452.Response.Rds') |> column_to_rownames(var = 'GENE_SYMBOL') 
threshold <- 0.8 * ncol(expr_mat)
detected_counts <- rowSums(expr_mat > 0)
expr_mat <- expr_mat[detected_counts >= threshold, ]
expr_mat <- log1p(expr_mat)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- read.csv('data/bulk_datasets/TIGER/Melanoma-phs000452/Melanoma-phs000452.Response.tsv', sep = '\t')
clin_info <- clin_info |> 
  mutate(time_os = overall.survival..days./30,
         os = ifelse(vital.status == 'Dead', 1, 0),
         patient = patient_name,
         sample = sample_id,
         tx = 'aPD1',
         tx_status = 'Baseline',
         cancertype = 'Melanoma',
         cohort = 'Cui') |> 
  dplyr::select(sample, patient, response, os, time_os, tx, tx_status, cancertype, cohort)
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/TIGER/Melanoma-phs000452/comb.csv', row.names = F)

# # Jerby-Arnon(30388455)
# expr_mat <- read.csv('data/bulk_datasets/Melanoma/Jerby-Arnon/ValCo2_tpm_resistance.csv', row.names = 'X')
# clin_info <- read.csv('data/bulk_datasets/Melanoma/Jerby-Arnon/ValCo2_clinicalAnno.csv', check.names = F)
# clin_info <- clin_info |> 
#   drop_na() |> 
#   mutate(time_pfs = `PFS (years)`*12,
#          pfs = `PFS (status, 1 = progression, 0 = no progression)`,
#          patient = Sample,
#          sample = Sample,
#          tx = 'aPD1',
#          cancertype = 'Melanoma',
#          cohort = 'Jerby-Arnon',
#          response = RECIST) |> 
#   select(patient, sample, response, pfs, time_pfs, tx, cancertype, cohort)

# Hugo
expr_mat <- readxl::read_xlsx('data/bulk_datasets/TIGER/Melanoma-GSE78220/GSE78220_PatientFPKM.xlsx') |> column_to_rownames(var = 'Gene')
threshold <- 0.8 * ncol(expr_mat)
detected_counts <- rowSums(expr_mat > 0)
expr_mat <- expr_mat[detected_counts >= threshold, ]
expr_mat <- log1p(expr_mat)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- read.csv('data/bulk_datasets/TIGER/Melanoma-GSE78220/Melanoma-GSE78220.Response.tsv', sep = '\t')
clin_info <- clin_info |> 
  mutate(time_os = overall.survival..days./30,
         os = vital.status,
         patient = patient_name,
         sample = paste0(patient, ifelse(Treatment == 'PRE', '.baseline', '.OnTx')),
         tx = 'aPD1',
         tx_status = 'Baseline',
         cancertype = 'Melanoma',
         cohort = 'Hugo') |> 
  select(patient, response, os, time_os, tx, cancertype, cohort, sample, tx_status) |> 
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/TIGER/Melanoma-GSE78220/comb.csv', row.names = F)

# Auslander (GSE115821)
expr_mat <- readRDS('data/bulk_datasets/TIGER/Melanoma-GSE115821/Melanoma-GSE115821.Response.Rds') |> column_to_rownames(var = 'GENE_SYMBOL')
threshold <- 0.9 * ncol(expr_mat)
detected_counts <- rowSums(expr_mat > 0)
expr_mat <- expr_mat[detected_counts >= threshold, ]
expr_mat <- log1p(expr_mat)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
colnames(score) <- str_replace(colnames(score), '.x', '')
score <- score[, !str_detect(colnames(score),'.y')]
clin_info <- read.csv('data/bulk_datasets/TIGER/Melanoma-GSE115821/Melanoma-GSE115821.Response.tsv', sep = '\t')
clin_info <- clin_info |> 
  mutate(time_os = NA,
         os = NA,
         patient = patient_name,
         sample = sample_id,
         tx = ifelse(Therapy == "anti-PD-1", 'aPD1', ifelse(Therapy == "anti-CTLA-4", 'aCTLA4', 'aPD1+CTLA4')),
         tx_status = ifelse(Treatment == 'PRE', 'Baseline','Treated'),
         cancertype = 'Melanoma',
         cohort = 'Auslander') |> 
  select(sample, patient, response, os, time_os, tx, cancertype, cohort, tx_status) |> 
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/TIGER/Melanoma-GSE115821/comb.csv', row.names = F)

# phs000452.v3.p1(31792460 Liu)
expr_mat <- read.table('data/bulk_datasets/Melanoma/Liu/41591_2019_654_MOESM3_ESM.txt') |> t() |> log1p()
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- readxl::read_xlsx('data/bulk_datasets/Melanoma/Liu/41591_2019_654_MOESM4_ESM.xlsx', skip = 2, sheet = 1, n_max = 144)
clin_info <- clin_info |> 
  mutate(time_os = OS/30,
         os = dead,
         time_pfs = PFS/30,
         pfs = progressed,
         sample = patient,
         tx = 'aPD1',
         tx_status = 'Baseline',
         cancertype = 'Melanoma',
         cohort = 'Liu',
         response = BR) |> 
  data.frame() |> 
  select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |> 
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/Melanoma/Liu/comb.csv', row.names = F)

# Ribas
expr_mat <- read.csv('data/bulk_datasets/Melanoma/Ribas/GSE158403_FPKM_matrix.txt.gz', sep = '\t', check.names = F) |> column_to_rownames(var = 'gene_id')
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- readxl::read_xlsx('data/bulk_datasets/Melanoma/Ribas/41467_2020_19810_MOESM3_ESM.xlsx')
clin_info <- clin_info |> pivot_wider(names_from = 'PARAMCD', values_from = c('AVAL','CNSR'))
clin_info <- clin_info |> 
  filter(Visit != 'NA') |> 
  mutate(time_os = AVAL_OS/30,
         os = CNSR_OS,
         time_pfs = AVAL_PFS/30,
         pfs = CNSR_PFS,
         sample = Sample.ID,
         patient = SUBJECT,
         tx = 'aPDL1',
         tx_status = ifelse(Visit == 'SCREENING', 'Baseline', 'Treated'),
         cancertype = 'Melanoma',
         cohort = 'Ribas',
         response = BESTRESP) |> 
  data.frame() |> 
  select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |> 
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/Melanoma/Ribas/comb.csv', row.names = F)

# Gide
expr_mat <- readRDS('data/bulk_datasets/TIGER/Melanoma-PRJEB23709/Melanoma-PRJEB23709.Response.Rds') |> 
  column_to_rownames(var = 'GENE_SYMBOL') |> 
  drop_na() 
detected_counts <- rowSums(expr_mat > 0)
expr_mat <- expr_mat[detected_counts >= 90, ]
expr_mat <- log1p(expr_mat)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- read.csv('data/bulk_datasets/TIGER/Melanoma-PRJEB23709/Melanoma-PRJEB23709.Response.tsv', sep = '\t')
clin_info <- clin_info |> 
  mutate(time_os = overall.survival..days./30,
         os = vital.status,
         time_pfs = NA,
         pfs = NA,
         sample = sample_id,
         patient = patient_name,
         tx = ifelse(Therapy == 'anti-PD-1', 'aPD1','aPD1+CTLA4'),
         tx_status = ifelse(Treatment == 'PRE', 'Baseline', 'Treated'),
         cancertype = 'Melanoma',
         cohort = 'Gide') |> 
  data.frame() |> 
  select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |> 
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/TIGER/Melanoma-PRJEB23709/comb.csv', row.names = F)

# LUD2015-005(EAC)
expr_mat <- readRDS('data/carroll_etal_2023/processed_data/bulk_RNAseq_counts/LUD2015-005_RNAseq_featureCounts.Rds')
e2s <- bitr(expr_mat$Geneid, fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
e2s <- e2s[!duplicated(e2s$SYMBOL),]
commongene <- intersect(expr_mat$Geneid, e2s$ENSEMBL)
expr_mat <- expr_mat[expr_mat$Geneid %in% commongene,-c(2:5)]
rownames(expr_mat) <- NULL
expr_mat <- expr_mat |> column_to_rownames(var = 'Geneid')
rownames(expr_mat) <- e2s$SYMBOL[match(rownames(expr_mat), e2s$ENSEMBL)]
genelength <- expr_mat$Length
expr_mat$Length <- NULL
expr_mat <- expr_mat/genelength
sum_column <- colSums(expr_mat)
expr_mat <- sweep(x = expr_mat*10e6, MARGIN = 2, STATS = sum_column, FUN = '/')
expr_mat <- log1p(expr_mat)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- readxl::read_xlsx('data/carroll_etal_2023/supplementary_files/Table_S8_papermetadata.xlsx')
# baseline
score_pre <- score[,str_detect(colnames(score), 'PreTx_Tumor')]
colnames(score_pre) <- str_replace(colnames(score_pre), '_Tumor', '')
clin_info_baseline <- clin_info |> 
  mutate(time_os = OS/30,
         os = Status,
         time_pfs = PFS/30,
         pfs = Progressed,
         patient = Patient,
         sample = paste0(patient, '_PreTx'),
         tx = 'aPDL1',
         tx_status = 'Baseline',
         cancertype = 'Esophageal',
         cohort = 'Carrol',
         response = Response_binary) |> 
  data.frame() |> 
  dplyr::select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |> 
  filter(sample %in% colnames(score_pre))
comb_baseline <- cbind(clin_info_baseline, t(score_pre)[clin_info_baseline$sample,])
# 4w
score_4w <- score[,str_detect(colnames(score), '_ICI-4W_Tumor')]
colnames(score_4w) <- str_replace(colnames(score_4w), '_Tumor', '')
clin_info_4w <- clin_info |> 
  mutate(time_os = OS/30,
         os = Status,
         time_pfs = PFS/30,
         pfs = Progressed,
         patient = Patient,
         sample = paste0(patient, '_ICI-4W'),
         tx = 'aPDL1',
         tx_status = '4w',
         cancertype = 'Esophageal',
         cohort = 'Carrol',
         response = Response_binary) |> 
  data.frame() |> 
  dplyr::select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |> 
  filter(sample %in% colnames(score_4w))
comb_4w <- cbind(clin_info_4w, t(score_4w)[clin_info_4w$sample,])
# 4w
score_post <- score[,str_detect(colnames(score), '_PostTx_Tumor')]
colnames(score_post) <- str_replace(colnames(score_post), '_Tumor', '')
clin_info_post <- clin_info |> 
  mutate(time_os = OS/30,
         os = Status,
         time_pfs = PFS/30,
         pfs = Progressed,
         patient = Patient,
         sample = paste0(patient, '_PostTx'),
         tx = 'aPDL1+chemotherapy',
         tx_status = 'Post',
         cancertype = 'Esophageal',
         cohort = 'Carrol',
         response = Response_binary) |> 
  data.frame() |> 
  dplyr::select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |> 
  filter(sample %in% colnames(score_post))
comb_post <- cbind(clin_info_post, t(score_post)[clin_info_post$sample,])
comb <- rbind(comb_baseline, comb_4w, comb_post)
write.csv(comb, 'data/carroll_etal_2023/processed_data/comb.csv', row.names = F)

# IMvigor210 (IMvigor210CoreBiologies)
expr_mat <- read.csv('data/bulk_datasets/IMvigor210/counts.csv')
gene_anno <- read.csv('data/bulk_datasets/IMvigor210/annotation.csv')
clin_info <- read.csv('data/bulk_datasets/IMvigor210/coldata.csv')
expr_mat$gene_symbol <- gene_anno$Symbol
expr_mat <- expr_mat[!(duplicated(expr_mat$gene_symbol)|is.na(expr_mat$gene_symbol)),]
gene_anno <- gene_anno[gene_anno$X %in% expr_mat$X,]
expr_mat$X <- NULL
rownames(expr_mat) <- NULL
expr_mat <- expr_mat |> column_to_rownames(var = 'gene_symbol')
expr_mat <- expr_mat/gene_anno$length
sum_column <- colSums(expr_mat)
expr_mat <- sweep(x = expr_mat*10e6, MARGIN = 2, STATS = sum_column, FUN = '/')
expr_mat <- log1p(expr_mat)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- clin_info |>
  mutate(time_os = os,
         os = censOS,
         time_pfs = NA,
         pfs = NA,
         patient = X,
         sample = X,
         response = Best.Confirmed.Overall.Response,
         tx = 'aPDL1',
         tx_status = 'Baseline',
         cancertype = 'mUC',
         cohort = 'IMvigor210') |>
  data.frame() |>
  dplyr::select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |>
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/IMvigor210/comb.csv', row.names = F)

# IMmotion150 (IMvigor210CoreBiologies)
expr_mat <- read.csv('data/bulk_datasets/IMmotion150/counts.csv')
gene_anno <- read.csv('data/bulk_datasets/IMmotion150/annotation.csv')
clin_info <- read.csv('data/bulk_datasets/IMmotion150/coldata.csv')
expr_mat$gene_symbol <- gene_anno$Symbol
expr_mat <- expr_mat[!(duplicated(expr_mat$gene_symbol)|is.na(expr_mat$gene_symbol)),]
gene_anno <- gene_anno[gene_anno$X %in% expr_mat$X,]
expr_mat$X <- NULL
rownames(expr_mat) <- NULL
expr_mat <- expr_mat |> column_to_rownames(var = 'gene_symbol')
expr_mat <- expr_mat/gene_anno$width
sum_column <- colSums(expr_mat)
expr_mat <- sweep(x = expr_mat*10e6, MARGIN = 2, STATS = sum_column, FUN = '/')
expr_mat <- log1p(expr_mat)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
colnames(score) <- str_replace(colnames(score), 'X','')
clin_info <- clin_info |>
  mutate(time_os = NA,
         os = NA,
         time_pfs = NA,
         pfs = NA,
         patient = X,
         sample = X,
         response = BestResponse,
         tx = 'aPDL1',
         tx_status = 'Baseline',
         cancertype = 'RCC',
         cohort = 'IMmotion150') |>
  data.frame() |>
  dplyr::select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |>
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/IMmotion150/comb.csv', row.names = F)

# RCC Braun
expr_mat <- readxl::read_xlsx('data/bulk_datasets/RCC_Braun/41591_2020_839_MOESM2_ESM.xlsx', sheet = 5, skip = 1) 
expr_mat <- expr_mat[!duplicated(expr_mat$gene_name),]
expr_mat <- expr_mat |> column_to_rownames(var = 'gene_name')
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- readxl::read_xlsx('data/bulk_datasets/RCC_Braun/41591_2020_839_MOESM2_ESM.xlsx', sheet = 1, skip = 1)
clin_info <- clin_info[clin_info$RNA_ID %in% colnames(expr_mat),] 
clin_info <- clin_info |> 
  mutate(time_os = OS,
         os = OS_CNSR,
         time_pfs = PFS,
         pfs = PFS_CNSR,
         patient = RNA_ID,
         sample = RNA_ID,
         response = ORR,
         tx = 'aPD1',
         tx_status = 'Baseline',
         cancertype = 'RCC',
         cohort = Cohort) |>
  data.frame() |>
  dplyr::select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |>
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/RCC_Braun/comb.csv', row.names = F)

# NSCLC_Kang(GSE218989)
expr_mat <- read.table('data/bulk_datasets/NSCLC_Kang/GSE218989_SMC_KAIST_TPM_matrix_pc_rmdupli.txt.gz', header = T) 
expr_mat <- expr_mat[!duplicated(expr_mat$Gene_Name),]
expr_mat <- expr_mat |> column_to_rownames(var = 'Gene_Name')
expr_mat <- log1p(expr_mat)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- readxl::read_xlsx('data/bulk_datasets/NSCLC_Kang/41467_2025_58068_MOESM3_ESM.xlsx')
clin_info <- clin_info |> 
  data.frame(check.names = F) |>
  mutate(time_os = `Overall survival (days)`/30,
         os = Death,
         time_pfs = `Progression-free survival (days)`/30,
         patient = PatientID,
         sample = PatientID,
         response = Responder,
         tx = 'aPD1',
         tx_status = 'Unspecified',
         cancertype = 'NSCLC',
         cohort = 'NSCLC_Kang') |>
  mutate(pfs = ifelse(time_os>time_pfs,1,0)) |> 
  dplyr::select(sample, patient, response, os, time_os, pfs, time_pfs, tx, tx_status, cancertype, cohort) |>
  filter(sample %in% colnames(score))
comb <- cbind(clin_info, t(score)[clin_info$sample,])
write.csv(comb, 'data/bulk_datasets/NSCLC_Kang/comb.csv', row.names = F)










