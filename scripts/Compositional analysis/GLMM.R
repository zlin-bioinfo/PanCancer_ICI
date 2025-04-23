rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','patchwork','ggplot2','tibble','qs','MetBrewer','forcats','grid','effsize','lme4','rstatix')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

meta_int <- read.csv('tables/meta_all.csv') 
unwanted_celltypes <- c('Melanocytes(CNA-)', 'Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Epithelial(CNA-)', 'Malignant(CNA+)','Cycling T','Cycling NK', "GCB-cycling", "PC-cycling", 'Cycling myeloids','Cycling non-immune')

# Only include paired samples
pt <- meta_int |>
  distinct(sample, .keep_all = TRUE) |>
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |>
  pull(patient)

freq_mat <- meta_int |> 
  filter(!celltype_r2 %in% unwanted_celltypes, patient %in% pt) |>
  distinct(celltype_r2, sample, .keep_all = T) |> 
  select(freq_r2_comp, cohort, patient, tx_status, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = 'tx_status', values_fill = 0) |> 
  pivot_longer(cols = c('Baseline', 'Treated'), names_to = 'tx_status', values_to = 'freq_r2_comp') |> 
  group_by(celltype_r2) |> 
  mutate(freq_r2_comp_scalelog = scale(log(freq_r2_comp + 1e-5))) |>
  ungroup()
freq_mat$tx_status <- factor(freq_mat$tx_status, c('Baseline','Treated'))
# Generalized linear mixed-effects models (GLMM)
uni_glmer <- function(freq_mat) {
  subtypes <- unique(freq_mat$celltype_r2)
  
  # Apply glmer model for each cell type
  uni_models <- lapply(subtypes, function(subtype) {
    print(subtype)
    
    # Filter for the current cell subtype
    subtype_data <- freq_mat |> filter(celltype_r2 == subtype)
  
    # Fit the model
    model <- tryCatch(
      glmer(tx_status ~ freq_r2_comp_scalelog + (1 | patient), 
            family = binomial(link = "logit"), 
            data = subtype_data),
      error = function(e) {
        warning(paste("Model failed for cell type:", subtype))
        return(NULL)
      }
    )
    
    # Skip subtypes without a fitted model
    if (is.null(model)) return(NULL)
    
    # Extract coefficients and confidence intervals
    coef_table <- coef(summary(model))
    confint_table <- confint(model, parm = "freq_r2_comp_scalelog", method = "Wald")
    
    # Extract values for freq_r2_comp_scalelog
    if (!"freq_r2_comp_scalelog" %in% rownames(coef_table)) return(NULL)
    
    subtype_results <- coef_table["freq_r2_comp_scalelog", , drop = FALSE]
    subtype_confint <- confint_table["freq_r2_comp_scalelog", , drop = FALSE]
    
    combined_results <- data.frame(
      Celltypes = subtype,
      Estimate = exp(subtype_results[1]),  # Odds ratio
      pValue = subtype_results[4],
      CI_lower = exp(subtype_confint[1]),
      CI_upper = exp(subtype_confint[2])
    )
    
    return(combined_results)
  })
  uni_models <- Filter(Negate(is.null), uni_models)
  
  # Combine results if non-empty
  if (length(uni_models) == 0) {
    message("No successful models fit.")
    return(NULL)
  }
  
  # Combine results
  df <- do.call(rbind, uni_models) |> data.frame()
  
  # Adjust for FDR
  df$fdr <- p.adjust(df$pValue, method = "fdr", n = nrow(df))
  
  # Categorize FDR results
  df <- df |> arrange(desc(Estimate))
  df$fdr_cat <- ifelse(df$fdr <= 0.05, "<0.05", ">0.05")
  
  return(df)
}

res_glmer <- uni_glmer(freq_mat)
write.csv(res_glmer, file = 'tables/GLMM_all.csv', row.names = F)
res_glmer1 <- read.csv('tables/GLMM_all.csv')
res_glmer <- filter(res_glmer, pValue < 0.05)

pdf('figures/Change/uni_glmer_r2.pdf', height = 5, width = 5)
p <- ggplot(res_glmer, aes(x= factor(Celltypes, levels = rev(res_glmer$Celltypes)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(aes(color = fdr_cat), size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.001), -log10(0.005), -log10(0.01)), labels = c('<0.001', '<0.005', '<0.01'), range = c(2,5)) +
  scale_y_continuous(name= "Odd Ratio", limits = c(0.5, 2.5)) +
  xlab("") +
  scale_color_manual(values = c("black", "grey")) +
  coord_flip() + theme_minimal() +
  labs(size = "P value", color = "FDR") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 10, colour = "black"))
grid.draw(p)
y_text = 0.14
x1 = 0.3
x2 = 0.75
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
grid.text("Baseline", x = unit(x1, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
grid.text("Treated", x = unit(x2, "npc"), y = unit(y_text, "npc"), gp = gpar(fontsize = 8))
dev.off()

#
subtypes <- res_glmer |> 
  filter(fdr<0.05) |> 
  select(Celltypes) |> 
  pull()
list_res <- lapply(subtypes, function(subtype){
  print(subtype)
  res_glmer |> filter(Celltypes == subtype)
  freq_mat_sub <- meta_int |> 
    filter(!celltype_r2 %in% unwanted_celltypes, patient %in% pt) |> 
    distinct(celltype_r2, sample, .keep_all = T) |> 
    select(freq_r2_comp, cohort, patient, time_point, celltype_r2) |> 
    pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
    pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
    filter(celltype_r2 == subtype) |> 
    group_by(cohort) |> 
    mutate(n_sample = n()/2) |> 
    filter(n_sample > 2) |>
    ungroup()
  uni_ds <- lapply(unique(freq_mat_sub$cohort), function(ds){
    mat <- freq_mat_sub[freq_mat_sub$cohort == ds,]
    print(ds)
    mat <- mat |> mutate(freq_r2_comp_scalelog = scale(log(freq_r2_comp + 1e-5)))
    mat$time_point <- ifelse(mat$time_point == 'Pre', 0, 1)
    model <- glmer(time_point ~ freq_r2_comp_scalelog + (1 | patient), 
                   family = binomial(link = "logit"), 
                   data = mat)
    model_summary <- summary(model)
    
    # Extract coefficients and confidence intervals
    coef_table <- coef(summary(model))
    confint_table <- confint(model, parm = "beta_", method = "Wald")
    
    # Extract values for freq_r2_comp_scale
    if (!"freq_r2_comp_scalelog" %in% rownames(coef_table)) return(NULL)
    
    subtype_results <- coef_table["freq_r2_comp_scalelog", , drop = FALSE]
    subtype_confint <- confint_table["freq_r2_comp_scalelog", , drop = FALSE]
    
    combined_results <- data.frame(
      Celltypes = subtype,
      Estimate = exp(subtype_results[1]),  # Odds ratio
      pValue = subtype_results[4],
      CI_lower = exp(subtype_confint[1]),
      CI_upper = exp(subtype_confint[2]),
      cohort = ds
    )
    return(combined_results)
  })
  res_df <- do.call(rbind, uni_ds)
  res_subtype_all <- res_glmer |> filter(Celltypes == subtype) |> select(!fdr_cat) 
  res_subtype_all$fdr <- 'Overall'
  names(res_subtype_all)[[6]] <- 'cohort'
  print(res_df)
  res_df <- rbind(res_subtype_all, res_df) 
})
names(list_res) <- subtypes

i <- 4
subtype <- subtypes[[i]]; subtype
res <- list_res[[i]]
min(res$CI_lower)
max(res$CI_upper)
res$Estimate[res$Estimate > 5] <- 5
res$CI_upper[res$CI_upper > 5] <- 5
order <- c("Overall", "SKCM_this study", "SKCM_Plozniak", "BCC&SCC_Yost",
           "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Shiao", "TNBC_Zhang",
           "HNSC_Franken", "HNSC_vanderLeun", "HNSC_Luoma",
           "CRC_Li","CRC_Chen", "NSCLC_Yan", "NSCLC_Liu","HCC_Guo")
pdf(paste0('figures/Change/', subtype, '.pdf'), height = 4, width = 5)
p <- ggplot(res, aes(x= factor(cohort, levels = rev(order)), y= Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, color = '#A4DBFF', position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=1, lty=2) +
  geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  scale_size(breaks = c(-log10(0.05), -log10(0.1), -log10(0.5)), labels = c('<0.05', '<0.1', '<0.5'), range = c(2,5)) +
  scale_y_continuous(name= "Odd Ratio", limits = c(0, 5),
                     breaks = c(0, 1, 2, 3, 4, 5), 
                     labels = c("0", "1", "2", "3", "4", ">5")) + 
  coord_flip() + theme_minimal() + ggtitle(unique(res$Celltypes)) +
  labs(size = "P value") + xlab("") +
  guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
y_text = 0.18
x1 = 0.3
x2 = 0.75
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

# boxplot
# ISG15 specific
subtypes <- as.factor(c('CD8_Tpex', 'Macro_TREM2', 'CD4_Tfh', 'Macro_SPP1'))
df_add <- meta_int |> 
  filter(patient %in% pt) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% subtypes) |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(direction = ifelse(On > Pre*1.1, 'Up', ifelse(On < Pre*0.9, 'Down', 'not'))) |>
  group_by(celltype_r2) |> 
  mutate(pt_count = n()) |> 
  group_by(direction, .add = T) |> 
  mutate(dir_count = n()) |> 
  distinct(celltype_r2, direction, pt_count, dir_count) |> 
  filter(!direction == 'not') |> 
  pivot_wider(values_from = 'dir_count', names_from = 'direction')
stat.test <- meta_int |> 
  filter(patient %in% pt) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% subtypes) |>  
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  group_by(celltype_r2) |> 
  wilcox_test(freq_r2_comp~time_point, paired = T) |> 
  adjust_pvalue(p.col = 'p', output.col = 'adj.p', method = "BH") |> 
  add_significance("adj.p")
stat.test <- merge(stat.test, df_add[df_add$celltype_r2 %in% subtypes ,], by = 'celltype_r2') 
stat.test <- stat.test |> 
  mutate(strip_label = paste0('Up: ', Up, '/', pt_count, ', Down: ', Down, '/', pt_count, '\n Adjusted p=', p)) 

meta_int |> 
  filter(patient %in% pt) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% subtypes) |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  select(freq_r2_comp, time_point, patient, celltype_r2) |> 
  merge(stat.test, by = 'celltype_r2') |> 
  ggplot(aes(x = factor(time_point, levels = c('Pre', 'On')), y = freq_r2_comp)) +
  geom_violin(aes(fill = time_point), alpha = 0.6) +
  geom_boxplot(width = 0.3, alpha = 0.2) +
  geom_line(aes(group = patient), linetype = "dashed", size = 0.2, alpha = 0.3) +
  geom_point(aes(color = time_point), size = 0.1) +
  scale_color_manual(values = c('#154999', '#CF0034')) +
  scale_fill_manual(values = c('#154999', '#CF0034')) +
  facet_wrap(.~factor(celltype_r2, levels = subtypes) + strip_label, scales = 'free', ncol = 6) +
  xlab("") + ylab("Relative Frequency") + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 12, colour = 'black'),
        strip.text = element_text(size = 10, margin = margin(2,1,2,1)),
        panel.spacing.x = unit(2, "lines")) +
  guides(fill="none", color = "none",
         color = guide_legend(override.aes = list(size=6))) 
ggsave('figures/Change/boxplot_top4.pdf', width =12, height = 3.5)

# Scissor validation
library(Scissor)
setwd("/home/zlin/workspace/PanCancer_ICI")
Scissor_m <- function(bulk_dataset, sc_dataset, phenotype, tag = NULL, alpha = NULL, 
                    cutoff = 0.2, family = c("gaussian", "binomial", "cox"), 
                    Save_file = "Scissor_inputs.RData", Load_file = NULL) 
{
  library(Seurat)
  library(Matrix)
  library(preprocessCore)
  if (is.null(Load_file)) {
    common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
    if (length(common) == 0) {
      stop("There is no common genes between the given single-cell and bulk samples.")
    }
    if (class(sc_dataset) == "Seurat") {
      sc_exprs <- as.matrix(sc_dataset@assays$RNA$data)
      network <- as.matrix(sc_dataset@graphs$RNA_snn)
    }
    else {
      sc_exprs <- as.matrix(sc_dataset)
      Seurat_tmp <- CreateSeuratObject(sc_dataset)
      Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", 
                                         verbose = F)
      Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
      Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), 
                           verbose = F)
      Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, 
                                  verbose = F)
      network <- as.matrix(Seurat_tmp@graphs$RNA_snn)
    }
    diag(network) <- 0
    network[which(network != 0)] <- 1
    dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, 
    ])
    print(dim(dataset0))
    dataset0 <- as.matrix(dataset0)
    dataset1 <- normalize.quantiles(dataset0)
    rownames(dataset1) <- rownames(dataset0)
    colnames(dataset1) <- colnames(dataset0)
    Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
    Expression_cell <- dataset1[, (ncol(bulk_dataset) + 1):ncol(dataset1)]
    X <- cor(Expression_bulk, Expression_cell)
    quality_check <- quantile(X)
    print("|**************************************************|")
    print("Performing quality-check for the correlations")
    print("The five-number summary of correlations:")
    print(quality_check)
    print("|**************************************************|")
    if (quality_check[3] < 0.01) {
      warning("The median correlation between the single-cell and bulk samples is relatively low.")
    }
    if (family == "binomial") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        print(sprintf("Current phenotype contains %d %s and %d %s samples.", 
                      z[1], tag[1], z[2], tag[2]))
        print("Perform logistic regression on the given phenotypes:")
      }
    }
    if (family == "gaussian") {
      Y <- as.numeric(phenotype)
      z <- table(Y)
      if (length(z) != length(tag)) {
        stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
      }
      else {
        tmp <- paste(z, tag)
        print(paste0("Current phenotype contains ", paste(tmp[1:(length(z) - 
                                                                   1)], collapse = ", "), ", and ", tmp[length(z)], 
                     " samples."))
        print("Perform linear regression on the given phenotypes:")
      }
    }
    if (family == "cox") {
      Y <- as.matrix(phenotype)
      if (ncol(Y) != 2) {
        stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
      }
      else {
        print("Perform cox regression on the given clinical outcomes:")
      }
    }
    save(X, Y, network, Expression_bulk, Expression_cell, 
         file = Save_file)
  }
  else {
    load(Load_file)
  }
  if (is.null(alpha)) {
    alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
               0.6, 0.7, 0.8, 0.9)
  }
  for (i in 1:length(alpha)) {
    set.seed(123)
    fit0 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, nlambda = 100, 
                  nfolds = min(10, nrow(X)))
    fit1 <- APML1(X, Y, family = family, penalty = "Net", 
                  alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
    if (family == "binomial") {
      Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
    }
    else {
      Coefs <- as.numeric(fit1$Beta)
    }
    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    percentage <- (length(Cell1) + length(Cell2))/ncol(X)
    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", 
                  length(Cell1), length(Cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", 
                  formatC(percentage * 100, format = "f", digits = 3)))
    if (percentage < cutoff) {
      break
    }
    cat("\n")
  }
  print("|**************************************************|")
  return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, 
                          family = family), Coefs = Coefs, Scissor_pos = Cell1, 
              Scissor_neg = Cell2))
}
pkgs <- c('GEOquery','clusterProfiler','org.Hs.eg.db','Seurat')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
# query <- 'GSE91061'
# getGEOSuppFiles(query, makeDirectory = TRUE, baseDir = 'data/',
#                 fetch_files = TRUE, filter_regex = 'fpkm')
# R.utils::gunzip('data/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz')
# df_expr <- read.csv('data/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv')
# 
# # convert entrezid to gene symbol
# e2s <- bitr(df_expr$X, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
# df_expr <- filter(df_expr, X %in% e2s$ENTREZID) %>% tibble::column_to_rownames(var = 'X')
# if (all.equal(rownames(df_expr),e2s$ENTREZID)){
#   rownames(df_expr) <- e2s$SYMBOL
# }
library(stringr)
count_mat <- data.table::fread('data/bulk_datasets/Melanoma/GSE115821/GSE115821/GSE115821_MGH_counts.csv.gz') 
count_mat <- count_mat[,-c(2:5)]
count_mat <- count_mat[!duplicated(count_mat$Geneid),]
count_mat <- count_mat |> tibble::column_to_rownames(var = 'Geneid')
gene_length <- count_mat$Length
count_mat <- count_mat[,-1]
expr_tpm <- apply(count_mat, 2, function(x) x/gene_length) |> 
  apply(2, function(x) x / sum(as.numeric(x)) * 10e6)
colnames(expr_tpm) <- colnames(expr_tpm) |> str_replace('MGH','') |> str_replace('.bam','') |> str_replace_all('-','') |> str_replace_all('_','')

phenotype <- readxl::read_xlsx('data/bulk_datasets/Melanoma/GSE115821/GSE115821/GSE115821.xlsx')
phenotype$Sample <- phenotype$Sample |> str_replace('MGH','') |> str_replace('.bam','') |> str_replace_all('-','') |> str_replace_all('_','')
phenotype <- phenotype[!duplicated(phenotype$Sample),]
phenotype <- tibble::column_to_rownames(phenotype, var = 'Sample')
phenotype$`Treatnent state` <- ifelse(str_detect(phenotype$`Treatnent state`, 'PRE'), 0, 1)
common <- intersect(colnames(expr_tpm), rownames(phenotype))
expr_tpm <- expr_tpm[,common]
phenotype <- phenotype[common,]
# phenotype <- data.frame(time_point = str_split(colnames(df_expr),'_',simplify = T)[,2])
# phenotype <- ifelse(phenotype$time_point == 'Pre', 0, 1)
seu <- qs_read('data/SKCM_Becker/seu_r2.qs2')
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample)
seu <- seu |> 
  NormalizeData() |>
  FindVariableFeatures()  |>
  ScaleData() |>
  RunPCA(verbose=FALSE) |>
  IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca",
                  new.reduction = 'harmony') |> 
  FindNeighbors(reduction = "harmony", dims = 1:20) |>
  FindClusters() |> 
  RunUMAP(dims = 1:20, reduction = 'harmony') |> 
  JoinLayers()
seu <- subset(seu, subset = celltype_main %in% c('CD4+T','CD8+T','Mono/macro') )
qs_save(seu, 'data/seu_scissor.qs2')
seu <- qs_read('data/seu_scissor.qs2')
infos1 <- Scissor_m(expr_tpm, seu, phenotype, alpha = 0.05, tag = c('Pre','On'),
                  family = "binomial", Save_file = 'data/Scissor_skcm_becker.RData')
infos1

seu$scissor <- 'NA'
seu$scissor[colnames(seu) %in% infos1$Scissor_pos] <- 'scissor+'
seu$scissor[colnames(seu) %in% infos1$Scissor_neg] <- 'scissor-'

seu_sub1 <- seu |> subset(subset = celltype_main == 'CD8+T' & time_point=='On' & scissor != 'NA')
df <- seu_sub@meta.data |> tabyl(celltype_r2,scissor)
df$`scissor+`>3
celltype <- df$celltype_r2[df$`scissor+`>3]
seu_sub1 <- seu_sub1 |> subset(subset = celltype_r2 %in% celltype)
DotPlot(seu_sub, group.by = 'scissor', features = c('PDCD1','LAG3','HAVCR2','CTLA4','SELL','CCR7','TCF7')) + RotatedAxis()

seu_sub2 <- seu |> subset(subset = celltype_main == 'CD4+T' & time_point=='On' & scissor != 'NA')
df <- seu_sub@meta.data |> tabyl(celltype_r2,scissor)
df$`scissor+`>3
celltype <- df$celltype_r2[df$`scissor+`>3]
seu_sub2 <- seu_sub2 |> subset(subset = celltype_r2 %in% celltype)
DotPlot(seu_sub, group.by = 'scissor', features = c('CXCL13','GNG4','CD200','CXCR5','BCL6','ICOS')) + RotatedAxis()

seu_sub3 <- seu |> subset(subset = celltype_main == 'Mono/macro' & time_point=='On' & scissor != 'NA')
df <- seu_sub@meta.data |> tabyl(celltype_r2,scissor)
df$`scissor+`>3
celltype <- df$celltype_r2[df$`scissor+`>3]
seu_sub3 <- seu_sub3 |> subset(subset = celltype_r2 %in% celltype)
DotPlot(seu_sub, group.by = 'scissor', features = c('C1QC','APOE','TREM2','GPNMB','MMP9','SPP1')) + RotatedAxis()
seu_merge <- merge(seu_sub2, y = c(seu_sub1, seu_sub3))
seu_merge$celltype_main <- factor(seu_merge$celltype_main, levels = c("CD4+T", "CD8+T","Mono/macro"))
grouped_features <- list(
  "Tfh markers" = c('CXCR5', 'BCL6', 'ICOS', 'SH2D1A'),  
  "Stem-like" = c('TCF7', 'SELL', 'CCR7', 'IL7R', 'LEF1'),  
  "Exhaustion" = c('PDCD1', 'LAG3', 'HAVCR2', 'CTLA4', 'ENTPD1', 'TOX'),  
  "SPP1-related" = c('SPP1', 'MMP9', 'GPNMB', 'FABP5', 'CHI3L1'),  
  "TREM2-related" = c('TREM2', 'APOE', 'C1QC', 'GPR34')  
)

DotPlot2(seu_merge, group.by = 'celltype_main', split.by = 'scissor', split.by.method = 'color', show_grid = F, 
         features = grouped_features,flip = TRUE,color_scheme='iwh_intense', dot.scale = 4) + 
  RotatedAxis() + ggtitle('Melanoma') + 
  scale_alpha_continuous(name = "Expression\n(scaled)") +
  guides(color = guide_legend(ncol=1, override.aes = list(size = 4)),
         alpha = guide_legend(ncol=4, override.aes = list(size = 4))) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        axis.text.x = element_text(color = "black", size = 10),  # X axis title color and size
        axis.text.y = element_text(color = "black", size = 12))
ggsave('figures/Scissor.pdf', height = 4, width = 8)


