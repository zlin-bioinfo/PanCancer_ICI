rm(list=ls())
pkgs <- c('dplyr','SingleCellExperiment','muscat','edgeR','qs','Seurat','scales','stringr','ggplot2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

sce <- qread('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/seu_r2.qs') |>
  subset(subset = celltype_main == 'Macro') |>
  as.SingleCellExperiment() 
# scale_exp <- assay(pb) |> apply(1, log1p) |> apply(1, rescale, to=c(0, 1)) 
# scale_pre <- scale_exp[, str_detect(colnames(scale_exp), 'Pre')] |> rowMeans()
# scale_post <- scale_exp[, str_detect(colnames(scale_exp), 'Post')] |> rowMeans()
# scale_mean <- cbind(scale_pre, scale_post) |> data.frame()

seu <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r2.qs') |>
  subset(subset = celltype_main == 'Macro') |> 
  NormalizeData() %>%
  FindVariableFeatures()%>%
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters()

datasets_mye <- c("SKCM_Becker", "BRCA_Bassez1", "BRCA_Bassez2", "TNBC_Zhang", "BCC_Yost", "HNSC_IMCISION", "HNSC_Luoma", "CRC_Li", "PCa_Hawley")
list_seu <- lapply(datasets_mye, function(dataset){
  seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) %>% 
    subset(subset = celltype_main %in%c('Macro'))
  return(seu)
})
list_seu[[5]]$dataset <- 'BCC_Yost'
seu <- merge(x=list_seu[[1]], y=list_seu[2:length(list_seu)]) 
seu <- JoinLayers(seu)
seu$int_cat <- ifelse(seu$interval > 21, 'Post', 'Early on')
sce <- as.SingleCellExperiment(seu)
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("sample"))

pbMDS(pb, )




pt_paired <- seu@meta.data |>
  group_by(sample) |> 
  mutate(count = n()) |> 
  ungroup() |> 
  filter(! count < 5) |> 
  distinct(sample, .keep_all = T) |> 
  group_by(patient) |> 
  mutate(sample_count = n()) |> 
  ungroup() |> 
  filter(sample_count == 2) |> 
  pull(patient) %>% unique()
pb <- seu |> subset(subset = patient %in% pt_paired) |> 
  Seurat2PB(sample = 'sample', cluster = 'celltype_main')
head(pb$samples, n=10L)
summary(pb$samples$lib.size)
keep.samples <- pb$samples$lib.size > 5e4
table(keep.samples)
pb <- pb[, keep.samples]
keep.genes <- filterByExpr(pb, group=pb$samples$cluster, 
                           min.count=5, min.total.count=10)
table(keep.genes)
pb <- pb[keep.genes, , keep=FALSE]
pb <- normLibSizes(pb)
head(pb$samples, n=10L)
df <- seu@meta.data |> distinct(sample, .keep_all = T) |> select(sample, response)
pb$sample$response <- df$response[match(pb$samples$sample, df$sample)]
pb$sample <- pb$sample |> mutate(response = case_when(response == 'RE'~ '1',
                                                               response == 'NR'~ '2'))

summary(pb$samples$norm.factors)
group <- as.factor(pb$samples$cluster)
plotMDS(pb, pch=16,  main="MDS")
legend("bottomright", legend=paste0("group", levels(group)),
       pch=16, col=2:8, cex=0.8)

pb <- estimateDisp(pb, design, robust=TRUE)


sce <- qread('/bigdata/zlin/Melanoma_meta/data/BRCA_Bassez1/seu_r2.qs') |>
  as.SingleCellExperiment() 
sce <- sce[rowSums(counts(sce)) > 5,]
sce <- pseudobulk(sce, group_by = vars(patient, time_point, celltype_main))
colData(sce)

colData(sce) %>%
  as_tibble() %>%
  group_by(patient, condition = time_point, celltype_main) %>%
  summarize(n_cells = n(), .groups = "drop") 
fit <- glm_gp(sce, design = ~ time_point + patient, size_factor = "ratio", verbose = TRUE)
res <- test_de(fit, contrast = cond(celltype_main = "CD4+T", time_point = "Post") - cond(celltype_main = "CD4+T", time_point = "Pre"))

ggplot(res, aes(x = lfc, y = - log10(pval))) +
  geom_point(aes(color = adj_pval < 0.01), size = 0.5)



