rm(list=ls())
pkgs <- c('qs2','tidyr','dplyr','plyr','stringr','tibble','janitor','Seurat', 'clusterProfiler','enrichplot','ggplot2','RColorBrewer','ComplexHeatmap','MetBrewer','GSVA','glmnet')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
source('scripts/Celltype_classification.R')
sample_info <- read.csv('tables/sample_info_time.csv')
# sample_info$sample <- str_replace_all(sample_info$sample, '_', '-')
datasets <- sample_info |> filter(subset == 'TME') |>  select(cohort) |> pull() |> unique()
datasets[which(datasets == "BCC&SCC_Yost")] <- 'BCC_Yost'
list_pseudobulk <- lapply(datasets, function(dataset) {
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = sample %in% sample_info$sample & celltype_r2 %in% immune) |> 
    # subset(subset = celltype_r2 %in% c('Epithelial(CNA-)','Epithelial(CNA+)','Melanocytes(CNA-)','Melanocytes(CNA+)','Malignant(CNA+)'), invert = T) |> 
    # subset(subset = celltype_main %in% c('Cycling','Cycling T/NK'), invert = T) |>
    # subset(subset = celltype_main %in% immune) |>
    NormalizeData() |> 
    AverageExpression(group.by = 'sample', layer = 'data', return.seurat = T) 
  seu$group <- sample_info$group[match(seu$sample, str_replace_all(sample_info$sample, '_', '-'))]
  return(seu)
})
seu <- merge(x = list_pseudobulk[[1]], y=list_pseudobulk[2:length(list_pseudobulk)])
seu <- JoinLayers(seu)
pseudobulk <- GetAssayData(seu, layer = 'data')

gs_list <- readxl::read_xlsx('tables/Signature.xlsx') |>
  lapply(function(x){return(x[!is.na(x)])})
# gs_list <- gs_list[-which(names(gs_list) %in% c("Th1","Th2","Effector cell traffic", "Cytotoxic cell inactivation", "Treg and Th2 traffic", "Metabolic suppression of CTL", "Breg cells", 
#                                                 "Myeloid cell traffic", "Macrophage and DC traffic", "pDC", "cDC1", "cDC2", "Neutrophil", "Granulocyte traffic",
#                                                 "Phagocytosis inhibition", "Myeloid suppression", "Immune suppression by myeloids", "Adipocytes", "Angiogenesis", "Stromal suppression", "Stromal immiune exclusion","Proliferation rate",
#                                                 "Hypoxia factors", "Glycolysis", "Acidosis", "Autophagy", "Senescence", "EMT Signature", "Metastasis Signature"))]
gs_list <- gs_list[which(names(gs_list) %in% c("T cells", "Effector cells", "Treg cells", "CD8+T cells", "NK cells", "B cells","Pan-macrophage","M1-like", "TAM",
                                               "Co-activation molecules", "Checkpoint molecules", "Interferon response"))] # ,
ssgseaPar <- ssgseaParam(as.matrix(pseudobulk), gs_list)
score <- gsva(ssgseaPar)
# gsvaPar <- gsvaParam(as.matrix(pseudobulk), gs_list)
# score <- gsva(gsvaPar)
df <- t(score) |> scale() |> data.frame(check.names = F)
df$group <- sample_info$group[match(rownames(df), str_replace_all(sample_info$sample, '_', '-'))]
df_mean <- aggregate(. ~ group, data = df, FUN = median) |> column_to_rownames(var = 'group') |> apply(1, as.numeric)
rownames(df_mean) <- names(gs_list)
Heatmap(df_mean)

# Prepare data
X <- as.matrix(df[, -which(names(df) %in% c("group"))])  # predictor matrix
y <- as.factor(df$group)  # target variable

# Train lasso multinomial logistic regression with cross-validation
cv_fit <- cv.glmnet(
  x = X,
  y = y,
  family = "multinomial",
  alpha = 1,            # alpha = 1 for lasso
  type.multinomial = "ungrouped",  # more flexible modeling
  nfolds = 10            # or set a higher number if you have more data
)

# Plot the cross-validation curve
plot(cv_fit)

# Best lambda value
best_lambda <- cv_fit$lambda.min
cat("Best lambda:", best_lambda, "\n")

# Coefficients at best lambda
coef(cv_fit, s = "lambda.min")



expr_mat <- read.csv('data/bulk_datasets/Melanoma/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv')
e2s <- bitr(expr_mat$X, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
expr_mat <- merge(e2s, expr_mat, by.x = 'ENTREZID', by.y = 'X')
expr_mat <- expr_mat[!duplicated(expr_mat$SYMBOL), ]
expr_mat <- expr_mat |> tibble::column_to_rownames(var = 'SYMBOL') |> subset(select = -ENTREZID)
expr_mat <- expr_mat[,!str_detect(colnames(expr_mat), 'Pt109_On')]
expr_mat <- log1p(expr_mat)
ssgseaPar <- ssgseaParam(as.matrix(expr_mat), gs_list)
score <- gsva(ssgseaPar)
score <- score |> t() |> scale() |> data.frame(check.names = F)



pred_probs <- predict(cv_fit, newx = as.matrix(score), s = "lambda.min", type = "response")
pred_labels <- predict(cv_fit, newx = as.matrix(score), s = "lambda.min", type = "class")
score <- score[order(pred_labels),]

col_ha = HeatmapAnnotation(
  group = sort(pred_labels),
  col = list(`TIME Subtype` = structure(names = as.character(unique(pred_labels)), met.brewer("Juarez", length(unique(pred_labels))))
  )
)
# genelist <- readxl::read_xlsx('tables/1-s2.0-S1535610821002221-mmc2.xlsx', sheet = 1, skip = 1)
# genelist_clean <- genelist[ , -1]
# gs_list <- apply(genelist_clean, 1, function(row) {
#   genes <- as.character(row)
#   genes <- genes[-1]
#   genes <- genes[!is.na(genes)]
#   return(genes)
# })
# names(gs_list) <- genelist_clean[[1]]

Heatmap(t(score), name = "Abundance\n(z-score)", column_title = '',
        cluster_rows = T, cluster_columns = F,
        show_row_names = T, show_column_names = F, 
        show_row_dend = F, show_column_dend = F,
        # column_order = order(sample_assignment$Subtype),
        # col = circlize::colorRamp2(c(-0.3, 0, 0.3), c("#154999", "white", "#CF0034")),
        na_col = 'lightgray',
        # heatmap_legend_param = list(
        #   title_position = "topcenter",
        #   at = c(-4, 0, 4), labels = c("Min", "", "Max"),
        #   legend_direction = "horizontal"),
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 8))

Heatmap(a |> t() |> scale() |> t(), name = "Abundance\n(z-score)", column_title = '',
        cluster_rows = T, cluster_columns = F,
        show_row_names = T, show_column_names = F, 
        show_row_dend = F, show_column_dend = F,
        # column_order = order(sample_assignment$Subtype),
        # col = circlize::colorRamp2(c(-0.3, 0, 0.3), c("#154999", "white", "#CF0034")),
        na_col = 'lightgray',
        # heatmap_legend_param = list(
        #   title_position = "topcenter",
        #   at = c(-4, 0, 4), labels = c("Min", "", "Max"),
        #   legend_direction = "horizontal"),
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 8))









# genelist <- readxl::read_xlsx('tables/1-s2.0-S1535610821002221-mmc2.xlsx', sheet = 1, skip = 1)
# genelist_clean <- genelist[ , -1]
# gs_list <- apply(genelist_clean, 1, function(row){
#   genes <- as.character(row)
#   genes <- genes[-1]
#   genes <- genes[!is.na(genes)]
#   return(genes)
# })
# gs_list
# names(gs_list) <- genelist_clean[[1]]
# gs_list$`Treg and Th2 traffic` <- NULL
# gs_list$`Neutrophil signature` <- NULL
# gs_list$`Tumor proliferation rate` <- NULL
# gs_list$`EMT signature` <- NULL
# top200_immune <- read.csv('tables/top200_deg_TIME.csv')
# top200_immune <- top200_immune |> filter(cluster != 'TIME-IE') |> select(gene) |> pull()
# top200_stromal <- read.csv('tables/top200_deg_TIME_nonimmune.csv')
# ie <- top200_stromal |> filter(cluster == 'TIME-IE') |> select(gene) |> pull()
# top200 <- c(top200_immune, ie) |> unique()


mtx_combined <- cbind(df[, -which(names(df) %in% c("group"))] |> t(), t(score))
M <- cor(mtx_combined)
mtx <- M[rownames(df), rownames(score)] |> data.frame(check.names = F)
mtx$group <- sample_info$group[match(rownames(mtx), str_replace_all(sample_info$sample, '_', '-'))]

means_group <- list()
for (i in 1:length(unique(mtx$group))){
  means <- mtx |> filter(group == unique(mtx$group)[[i]]) |> select(!group) |> miscTools::colMedians()
  means_group[[i]] <- means
}
names(means_group) <- unique(mtx$group)
summary_means <- do.call(rbind, means_group)

transposed_data <- as.data.frame(t(summary_means))

sample_assignment <- transposed_data |> 
  scale() |>
  data.frame(check.names = F) |> 
  rownames_to_column(var = "sample") |> 
  pivot_longer(
    cols = -sample,
    names_to = "Subtype",
    values_to = "Score"
  ) |> 
  group_by(sample) |> 
  slice_max(Score, n = 1) |>   # Keeps only the row with the highest score per patient
  ungroup()
sample_assignment <- sample_assignment |> 
  arrange(Subtype) |> 
  filter(Score > 0.25)

col_ha = HeatmapAnnotation(
  group = sample_assignment$Subtype,
  col = list(`TIME Subtype` = structure(names = as.character(unique(sample_assignment$Subtype)), met.brewer("Juarez", length(unique(sample_assignment$Subtype))))
  )
)
# genelist <- readxl::read_xlsx('tables/1-s2.0-S1535610821002221-mmc2.xlsx', sheet = 1, skip = 1)
# genelist_clean <- genelist[ , -1]
# gs_list <- apply(genelist_clean, 1, function(row) {
#   genes <- as.character(row)
#   genes <- genes[-1]
#   genes <- genes[!is.na(genes)]
#   return(genes)
# })
# names(gs_list) <- genelist_clean[[1]]
mat <-expr_mat[intersect(rownames(expr_mat),c(gs_list$`B cells`, gs_list$Endothelium, gs_list$`Tumor-associated Macrophages`, gs_list$`Effector cell traffic`, gs_list$Treg)), sample_assignment$sample] |> t() |> scale() |> t()

Heatmap(a[,sample_assignment$sample] |> t() |> scale() |> t(), name = "Abundance\n(z-score)", column_title = '',
        cluster_rows = T, cluster_columns = F,
        show_row_names = T, show_column_names = F, 
        show_row_dend = F, show_column_dend = F,
        # column_order = order(sample_assignment$Subtype),
        # col = circlize::colorRamp2(c(-0.3, 0, 0.3), c("#154999", "white", "#CF0034")),
        na_col = 'lightgray',
        # heatmap_legend_param = list(
        #   title_position = "topcenter",
        #   at = c(-4, 0, 4), labels = c("Min", "", "Max"),
        #   legend_direction = "horizontal"),
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 8))




library(xCell)
a <- xCellAnalysis(expr_mat)

gs_list <- readxl::read_xlsx('tables/Signature.xlsx') |>
  lapply(function(x){return(x[!is.na(x)])})

