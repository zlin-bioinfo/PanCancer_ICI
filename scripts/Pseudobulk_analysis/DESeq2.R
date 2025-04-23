pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer','DESeq2','apeglm')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(max.print = 10000)
options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)

# 
cd4_t <- c('CD4_Tcm','CD4_Treg','CD4_T-ISG','CD4_Tfh','CD4_Tstr','CD4_Tctl','CD4_Th17','CD4_T-naive')
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 'SCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'NSCLC_Liu',
              'PCa_Hawley', 'HCC_Guo','RCC_Bi')
seu_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) |> 
    subset(subset = celltype_r2 %in% cd4_t) 
  for .sample in (unique(seu$sample)){
    if(seu@meta.data |> filter(sample == .sample) |> nrow() > 1){
      
    }
  }
  return(seu)
})

subs <- s@meta.data[s@meta.data$scGate_multi == c,]


seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) |> 
  subset(subset = celltype_main %in% c("NK","CD8+T","CD4+T")) 
metadata <- seu@meta.data
available_sample <- metadata |> 
  group_by(sample, celltype_main) |> 
  mutate(cell_count = n()) |> 
  filter(cell_count >= 30) |> 
  pull(sample) |> 
  unique()
seu <- subset(seu, subset = sample %in% available_sample)
Idents(seu) <- seu$response
seu <- NormalizeData(seu)
res.de <- seu |> FindMarkers( ident.1 = 'RE', ident.2 = 'NR', test.use = 'wilcox')

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
sce <- seu |> as.SingleCellExperiment()
reduced_sce <- pseudobulk(sce, group_by = vars(sample, patient, response, condition = time_point), n_cells = n())
reduced_sce <- reduced_sce[,reduced_sce$n_cells>=50]
fit <- glm_gp(reduced_sce, design = ~ condition + response, size_factor = "ratio", verbose = TRUE)
summary(fit)
colnames(fit$Beta)
res <- test_de(fit, contrast = cond(condition = "On", response = 'RE') - cond(condition = "On", response='NR'))
ggplot(res, aes(x = lfc, y = - log10(pval))) +
  geom_point(aes(color = adj_pval < 0.05), size = 0.5)
up <- res |> filter(adj_pval<0.05 & lfc>1) |> dplyr::select(name) |> pull()
gg <- bitr(up,'SYMBOL','ENTREZID','org.Hs.eg.db')
ggo <- enrichKEGG(gg$ENTREZID,organism = 'hsa')
dotplot(ggo)
library(msigdbr)
h_gene_sets = msigdbr(species = "human", category = "H")
library(fgsea)

msigdbr_list = split(x = h_gene_sets$gene_symbol, f = h_gene_sets$gs_name)
res.ora <- fora(genes=up, universe=res$name, pathways=msigdbr_list)
res.ora |> arrange(desc(foldEnrichment)) |> head()

kk <- compareCluster(ENTREZID~cluster, data = top20, fun=enrichKEGG)

avg_expr <- AverageExpression(seu, layer = 'counts', group.by = 'sample')
coldata <- metadata |> 
  distinct(sample, .keep_all = T) |> 
  select(cohort, sample, patient, response, time_point) |> 
  filter(sample %in% available_sample)
coldata$sample <- str_replace_all(coldata$sample, '_', '-')
avg_expr <- avg_expr$RNA |> data.frame(check.names = F) |> mutate_if(is.numeric, ~ifelse(is.na(.), 0, .)) |> data.frame(check.names = F)
avg_expr <- avg_expr[, coldata$sample]
rownames(coldata) <- coldata$sample
dds <- DESeq2::DESeqDataSetFromMatrix(countData = avg_expr,
                                      colData = coldata,
                                      design = ~ sample + response + time_point) 
#Remove genes with very low counts:
dds <- dds[rowSums(counts(dds)) > 6, ]

# do Differential expression
de <- DESeq2::DESeq(dds)

# shrink using apelgm
res <- DESeq2::lfcShrink(de,
                         coef = DESeq2::resultsNames(de)[3],
                         type = "apeglm")

res.list[[c]] <- as.data.frame(res) 

