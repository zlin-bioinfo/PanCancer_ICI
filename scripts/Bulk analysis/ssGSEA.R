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
# Melanoma 
# GSE91061
expr_mat <- read.csv('data/bulk_datasets/TIGER/Melanoma-GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz')
# convert entrezid to gene symbol
e2s <- bitr(expr_mat$X, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
expr_mat <- merge(e2s, expr_mat, by.x = 'ENTREZID', by.y = 'X')
expr_mat <- expr_mat[!duplicated(expr_mat$SYMBOL), ]
expr_mat <- expr_mat |> tibble::column_to_rownames(var = 'SYMBOL') |> subset(select = -ENTREZID)
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
score <- score[,!str_detect(colnames(score), 'Pt109_On')]
colnames(score) <- paste0(str_split(colnames(score), '_', simplify = T)[,1], '_', str_split(colnames(score), '_', simplify = T)[,2])

clin_info <- readxl::read_xlsx('data/bulk_datasets/Melanoma/GSE91061/1-s2.0-S0092867417311224-mmc2.xlsx', skip = 1)
clin_info <- clin_info |> 
  data.frame(check.names = F) |> 
  mutate(time_os = `Time to Death\r\n(weeks)`/4,
         os = ifelse(`Dead/Alive\r\n(Dead = True)` == TRUE, 1, 0),
         response = Response,
         patient = Patient,
         tx = 'aPD1',
         cancertype = 'Melanoma',
         cohort = 'Riaz') |> 
  select(patient, response, os, time_os, tx, cancertype, cohort)

# Nathanson(27956380)
expr_mat <- readRDS('data/bulk_datasets/TIGER/Melanoma-Nathanson_2017/Melanoma-Nathanson_2017.Response.Rds') |> column_to_rownames(var = 'GENE_SYMBOL')
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- read.csv('data/bulk_datasets/TIGER/Melanoma-Nathanson_2017/Melanoma-Nathanson_2017.Response.tsv', sep = '\t')
# No paired samples
clin_info <- clin_info |> 
  mutate(time_os = overall.survival..days./30,
         os = ifelse(vital.status == 'Dead', 1, 0),
         patient = patient_name,
         sample = sample_id,
         tx = 'aCTLA4',
         tx_status = ifelse(Treatment == 'PRE', 'Baseline', 'Treated'),
         cancertype = 'Melanoma',
         cohort = 'Nathanson') |> 
  select(patient, sample, response, os, time_os, tx, tx_status, cancertype, cohort)

# Van Allen (phs000452.v2.p1)
expr_mat <- readRDS('data/bulk_datasets/TIGER/Melanoma-phs000452/Melanoma-phs000452.Response.Rds') |> column_to_rownames(var = 'GENE_SYMBOL')
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- read.csv('data/bulk_datasets/TIGER/Melanoma-phs000452/Melanoma-phs000452.Response.tsv', sep = '\t')
clin_info <- clin_info |> 
  mutate(time_os = overall.survival..days./30,
         os = ifelse(vital.status == 'Dead', 1, 0),
         patient = patient_name,
         tx = 'aPD1',
         cancertype = 'Melanoma',
         cohort = 'Cui') |> 
  select(patient, response, os, time_os, tx, cancertype, cohort)

# Jerby-Arnon(30388455)
expr_mat <- read.csv('data/bulk_datasets/Melanoma/Jerby-Arnon/ValCo2_tpm_resistance.csv', row.names = 'X')
ssgseaPar <- gsvaParam(as.matrix(expr_mat), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
clin_info <- read.csv('data/bulk_datasets/Melanoma/Jerby-Arnon/ValCo2_clinicalAnno.csv', check.names = F)
clin_info <- clin_info |> 
  drop_na() |> 
  mutate(time_pfs = `PFS (years)`*12,
         pfs = `PFS (status, 1 = progression, 0 = no progression)`,
         patient = Sample,
         sample = Sample,
         tx = 'aPD1',
         cancertype = 'Melanoma',
         cohort = 'Jerby-Arnon',
         response = RECIST) |> 
  select(patient, sample, response, pfs, time_pfs, tx, cancertype, cohort)
# Hugo
expr_mat <- readxl::read_xlsx('data/bulk_datasets/TIGER/Melanoma-GSE78220/GSE78220_PatientFPKM.xlsx') |> column_to_rownames(var = 'Gene')
clin_info <- read.csv('data/bulk_datasets/TIGER/Melanoma-GSE78220/Melanoma-GSE78220.Response.tsv', sep = '\t')
clin_info <- clin_info |> 
  mutate(time_os = overall.survival..days./30,
         os = vital.status,
         patient = patient_name,
         sample = paste0(patient, ifelse(Treatment == 'PRE', '.baseline', '.OnTx')),
         tx = 'aPD1',
         cancertype = 'Melanoma',
         cohort = 'Hugo') |> 
  select(patient, response, os, time_os, tx, cancertype, cohort)

# Gide(30753825)
expr_mat <- readRDS('data/bulk_datasets/TIGER/Melanoma-GSE115821/Melanoma-GSE115821.Response.Rds') |> column_to_rownames(var = 'GENE_SYMBOL')
clin_info <- read.csv('data/bulk_datasets/TIGER/Melanoma-GSE115821/Melanoma-GSE115821.Response.tsv', sep = '\t')
clin_info <- clin_info |> 
  mutate(time_os = NA,
         os = NA,
         patient = patient_name,
         sample = sample_id,
         tx = ifelse(Therapy == "anti-PD-1", 'aPD1', ifelse(Therapy == "anti-CTLA-4", 'aCTLA4', 'aPD1+CTLA4')),
         cancertype = 'Melanoma',
         cohort = 'Gide') |> 
  select(sample, patient, response, os, time_os, tx, cancertype, cohort)

# phs000452.v3.p1(31792460 Liu)
# No paired
expr_mat <- read.table('data/bulk_datasets/Melanoma/Liu/41591_2019_654_MOESM3_ESM.txt') |> t() |> log1p()
clin_info <- readxl::read_xlsx('data/bulk_datasets/Melanoma/Liu/41591_2019_654_MOESM4_ESM.xlsx', skip = 2, sheet = 1, n_max = 144)
clin_info <- clin_info |> 
  mutate(time_os = OS,
         os = dead,
         time_pfs = PFS,
         pfs = progressed,
         sample = patient,
         tx = ifelse(priorCTLA4 == 1, 'aCTLA4', 'aPD1'),
         tx_status = ifelse(priorCTLA4 == 1, 'Treated', 'Baseline'),
         cancertype = 'Melanoma',
         cohort = 'Liu',
         response = BR) |> 
  data.frame() |> 
  select(sample, patient, response, os, time_os, pfs, time_pfs, tx, cancertype, cohort)




df_expr <- read.csv('/home/zlin/workspace/PanCancer_ICI/data/bulk_datasets/Melanoma/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv', row.names = 'X')
clin <- readxl::read_xlsx('/home/zlin/workspace/PanCancer_ICI/data/bulk_datasets/Melanoma/GSE91061/1-s2.0-S0092867417311224-mmc2.xlsx', skip = 1)
# convert entrezid to gene symbol
e2s <- bitr(rownames(df_expr), fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
e2s <- e2s[!duplicated(e2s),]
df_expr <- df_expr[e2s$ENTREZID,]
rownames(df_expr) <- e2s$SYMBOL
df_expr <- log1p(df_expr)
genelist <- c(read.csv('/home/zlin/workspace/PanCancer_ICI/tables/marker_cd4t.csv', check.names = F) |> as.list(),
              read.csv('/home/zlin/workspace/PanCancer_ICI/tables/marker_cd8tnk.csv', check.names = F) |> as.list(),
              read.csv('/home/zlin/workspace/PanCancer_ICI/tables/marker_bplasma.csv', check.names = F) |> as.list(),
              read.csv('/home/zlin/workspace/PanCancer_ICI/tables/marker_myeloids.csv', check.names = F) |> as.list(),
              read.csv('/home/zlin/workspace/PanCancer_ICI/tables/marker_nonimmune.csv', check.names = F) |> as.list())
genelist$Cycling <- NULL
ssgseaPar <- gsvaParam(as.matrix(df_expr), genelist, minSize=5, maxSize=500)
score <- gsva(ssgseaPar)
score_pre <- score[,str_detect(colnames(score), pattern = 'Pre')]
colnames(score_pre) <- str_split(colnames(score_pre), '_', simplify = T)[,1]
common_pt <- intersect(colnames(score_pre), clin$Patient)

t(score_pre)[match(clin$Patient, common_pt), ]

df <- cbind(clin[match(common_pt, clin$Patient), ], t(score_pre))
df$purity <- TumorPurity

# Tumor purity
df_expr <- df_expr[,str_detect(colnames(df_expr), pattern = 'Pre')]
estimate <- function(dat,pro){
  input.f=paste0('/home/zlin/workspace/PanCancer_ICI/data/bulk_datasets/Melanoma/', pro,'_estimate_input.txt')
  output.f=paste0('/home/zlin/workspace/PanCancer_ICI/data/bulk_datasets/Melanoma/', pro,'_estimate_gene.gct')
  output.ds=paste0('/home/zlin/workspace/PanCancer_ICI/data/bulk_datasets/Melanoma/', pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## platform
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
scores=estimate(df_expr, 'GSE91061')
TumorPurity = cos(0.6049872018 + 0.0001467884 * scores[,3])
df$purity <- TumorPurity
library(survminer)
library(survival)
table(df$`Dead/Alive\r\n(Dead = True)`, useNA='ifany')
names(df)[4] <- 'status'
names(df)[5] <- 'time'
df$status <- ifelse(df$status == TRUE, 1, 0)
df$time <- as.numeric(df$time)
df$surv_obj <- Surv(df$time, df$status == 1)

celltypes <- rownames(score)
list_res <- lapply(celltypes, function(celltype){
  formula <- as.formula(paste0('surv_obj ~ ', celltype))
  cox_model <- coxph(formula, data = df)
  model_summary <- summary(cox_model)
  # Extracting values
  hr <- exp(model_summary$coefficients[, "coef"])  # Exponentiate coef
  p_value <- model_summary$coefficients[, "Pr(>|z|)"]
  conf_int <- model_summary$conf.int[, c("lower .95", "upper .95")]
  # Combine into a named vector
  c(HR = hr, lower95 = conf_int[1], upper95 = conf_int[2], p_value = p_value)
})
names(list_res) <- celltypes
res <- do.call(rbind,list_res)
# View the summary of the Cox model












