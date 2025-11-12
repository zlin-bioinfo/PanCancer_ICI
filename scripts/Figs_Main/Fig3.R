pkgs <- c('qs2','tidyr','dplyr','plyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','tibble','pheatmap','MetBrewer','viridis','ComplexHeatmap','colorRamp2','corrr','ggnewscale','NMF','ggpubr','corrplot','rstatix','janitor','Seurat','survminer','survival')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

nmf.res <- qs_read('data/res_nmf_rank4_seed_1234.qs2')
w <- basis(nmf.res)
colnames(w) = c('TIME-Mye','TIME-Q','TIME-I', 'TIME-B')
w_t = w |> 
  t() |>  
  scale()

pdf('figures/Fig3/ht_rank4.pdf', height = 2.5, width = 13)
Heatmap(w_t, width = 11, height = 4,  cluster_rows = T, border = F,
        row_names_side = "right", 
        column_names_rot = 45, clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2',
        col = circlize::colorRamp2(c(-1.1, 0, 1.1), c("#154999", "white", "#CF0034")), na_col = 'lightgray',
        heatmap_legend_param = list(title = "Loading \n(scaled)", 
                                    at = c(-1.1, 0, 1.1), 
                                    labels = c("Min", "", "Max"),
                                    border = "black", 
                                    legend_direction = "horizontal",
                                    title_position = "topcenter"),
        rect_gp = gpar(col = "white", lwd = 0.5),
        column_names_gp = gpar(fontsize = 8.5),
        row_names_gp = gpar(fontsize = 11)) 
dev.off()

sample_group <- predict(nmf.res) |> as.data.frame()
names(sample_group) <- 'TIME'
sample_info <- metadata |> distinct(sample, .keep_all = T) |> filter(sample %in% rownames(sample_group))
sample_info$TIME <- sample_group$TIME[match(sample_info$sample, rownames(sample_group))]
sample_info <- sample_info |> 
  mutate(TIME = factor(case_when(TIME == 1 ~ 'TIME-Mye', TIME == 3 ~ 'TIME-I', TIME == 2 ~ 'TIME-Q', TIME == 4 ~ 'TIME-B'), 
                       levels = c('TIME-Q', 'TIME-I', 'TIME-B', 'TIME-Mye')))
sample_info$TIME <- factor(sample_info$TIME, levels = c('TIME-Q', 'TIME-I', 'TIME-B', 'TIME-Mye'))
sample_info$group <- sample_info$TIME
pt <- sample_info |>
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |> 
  pull(patient) 
sample_info$paired <- ifelse(sample_info$patient %in% pt, 'Yes', 'No')
write.csv(sample_info, 'tables/sample_info_time_updated.csv', row.names = F)

sample_info <- read.csv('tables/sample_info_time_updated.csv')
unmatched_pt <- c("SKCM_this study_Patient2", "SKCM_this study_Patient3")
sample_info <- sample_info |> 
  mutate(response = factor(response, levels = c('R','NR','NE')),
         tx_status = factor(tx_status, levels =  c('Baseline','Treated')),
         subtype = factor(subtype, levels = c("SKCM", "BCC", "BRCA(ER/HER+)", "TNBC", "HNSC",  "NSCLC", "HCC", "CRC", "RCC", "PCa"))) |> 
  arrange(TIME, tx_status, response, subtype)
col_ha = HeatmapAnnotation(
  `Cancer Type` = sample_info$subtype,
  Response = sample_info$response,
  `Treatment Status` = sample_info$tx_status,
  `TIME Subtype` = sample_info$TIME,
  col = list(`Cancer Type` = structure(met.brewer('Austria', length(unique(sample_info$subtype))), names = levels(sample_info$subtype)),
             Response = c('R' = '#CC0C00FF','NR' = '#5C88DAFF','NE' = '#84BD00FF'),
             `Treatment Status` = c('Baseline' = "#04a3bd", 'Treated' = "#f0be3d"),
             `TIME Subtype` = structure(met.brewer('Juarez', length(unique(sample_info$TIME))), names = levels(sample_info$TIME))
  ),
  simple_anno_size = unit(0.4, "cm"), annotation_name_gp= gpar(fontsize = 10)
)
pal <- colorRampPalette(brewer.pal(10, "RdYlBu"))
celltype <- c(
  "CD4_T-naive","B-naive","CD8_T-naive","B-memory", "CD4_Tctl","CD8_Trm", "CD8_Temra", "NK_CD56loCD16hi", 
  "CD4_T-ISG", "CD8_T-ISG","B-ISG","mregDC", "CD4_Tfh","CD4_Treg", 'CD8_Tpex',"CD8_Tex_GZMK","CD8_Tex_CXCL13",
  "CD4_Tstr", "CD8_Tstr", "B-HSP", "ACB_NR4A2", "ACB_EGR1",'ACB_CCR7', "B-AtM", "PC_IGHA", "PC_IGHG",
  "Cycling myeloids", 'cDC1',"cDC2", "Mono_CD14", "Mono_CD16", "Macro_LYVE1", "Macro_SPP1","Macro_TREM2", "Macro_C1QC"
)
celltype.group <- factor(c(
  rep("TIME-Q",8),
  rep("TIME-I",9),
  rep("TIME-B",9),
  rep("TIME-Mye",9)), levels = c('TIME-Q', 'TIME-I', 'TIME-B', 'TIME-Mye'))
zscale_mtx <- apply(mtx_immune[sample_info$sample, celltype], MARGIN = 2, scale)
zscale_mtx <- t(zscale_mtx)

df.row<- data.frame(group = celltype.group); rownames(df.row) <- celltype
row_annotation <- HeatmapAnnotation(df = df.row, 
                                    col = list(group = structure(met.brewer('Juarez', length(levels(sample_info$TIME))), names = levels(sample_info$TIME))),
                                    which = 'row',
                                    show_legend = F, show_annotation_name = F,
                                    simple_anno_size = grid::unit(3, "cm"))

# heatmap
pdf('figures/Fig3/ht_sample_rank4.pdf', height = 5.5, width = 12)
Heatmap(zscale_mtx, name = "Abundance\n(z-score)", column_title = 'Pan-cancer TIME Classification (n=418)', 
        column_split = sample_info$TIME, 
        row_split = celltype.group,
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F, 
        show_row_dend = F, show_column_dend = F,
        column_order = order(sample_info$subtype, sample_info$response, sample_info$tx),
        col = circlize::colorRamp2(c(-3, 0, 3), c("#154999", "white", "#CF0034")),
        na_col = 'lightgray',
        heatmap_legend_param = list(
          border = T,
          title_position = "topcenter",
          at = c(-3, 0, 3), labels = c("Min", "", "Max"),
          legend_direction = "horizontal"),
        top_annotation = col_ha,
        right_annotation = row_annotation,
        row_names_gp = gpar(fontsize = 8))
dev.off()

# TCGA
expr_mat <- data.table::fread('data/TCGA/tcga_RSEM_gene_tpm.gz', header = T)
expr_mat <- column_to_rownames(expr_mat, var = 'sample')
annotation <- data.table::fread('data/TCGA/probeMap_gencode.v23.annotation.gene.probemap', header = T)
annotation <- annotation[!duplicated(annotation$gene),]
expr_mat <- expr_mat[intersect(rownames(expr_mat), annotation$id),]
rownames(expr_mat) <- annotation$gene[match(rownames(expr_mat), annotation$id)]
expr_mat <- expr_mat[,str_detect(colnames(expr_mat), '-01')]
colnames(expr_mat) <- str_replace(colnames(expr_mat), '-01', '')
expr_mat <- expr_mat[, !colSums(is.na(expr_mat)) > 0]
expr_mat <- expr_mat[,-which(colnames(expr_mat) %in% 'TCGA-25-1870')]

# Thorsson et al 2018(https://doi.org/10.1016/j.immuni.2018.03.023)
immune_subtype <- data.table::fread('data/TCGA/Subtype_Immune_Model_Based.txt.gz', header = T)
immune_subtype <- immune_subtype[str_detect(immune_subtype$sample, '-01'),]
immune_subtype$sample <- str_replace(immune_subtype$sample, '-01', '')

# Bagaev et al 2021(https://doi.org/10.1016/j.ccell.2021.04.014)
df_Bagaev <- readxl::read_xlsx('data/TCGA/1-s2.0-S1535610821002221-mmc6.xlsx', sheet = 2, skip = 1)
# df_Bagaev <- df_Bagaev |> filter(Purity < 0.9)
sample_common <- intersect(colnames(expr_mat), df_Bagaev$Sample)
df_Bagaev <- df_Bagaev[df_Bagaev$Sample %in% sample_common,]

# Ecotype
df_ecotype <- readxl::read_xlsx('data/TCGA/1-s2.0-S0092867421010618-mmc6.xlsx', sheet = 2, skip = 4)
colnames(df_ecotype)[1] <- 'sample'
colnames(df_ecotype)[12] <- 'assignment'
df_ecotype <- df_ecotype |> filter(str_detect(sample, '01A'))
df_ecotype$sample <- paste(str_split(df_ecotype$sample, '\\.', simplify = T)[, 1],
                           str_split(df_ecotype$sample, '\\.', simplify = T)[, 2],
                           str_split(df_ecotype$sample, '\\.', simplify = T)[, 3], sep = '-')

gs_list <- readxl::read_xlsx('tables/Signature.xlsx') |>
  lapply(function(x){return(x[!is.na(x)])})
sample_info <- read.csv('tables/sample_info_time_updated.csv')
seu <- qs_read('data/pseudobulk_filtered.qs2')

# GSVA
gs_sub <- c("T cells","Treg cells", "Effector T cells","NK cells", "Th1","Th2","Th17", 
            "B cells", "Plasma cells", "MHC-I", "MHC-II", "TLS formation", 
            "DC", "Pan-macrophage", "M1-like", "M2-like",
            "Co-activation molecules", "Checkpoint molecules", "Interferon response",
            "Endothelium", "Angiogenesis", "CAF", "Matrix remodeling","Stromal immiune exclusion","Stromal suppression")
ssgseaPar <- ssgseaParam(GetAssayData(seu, layer = 'data'), gs_list)
score_pb <- gsva(ssgseaPar)
score_pb <- score_pb |> t() |> data.frame(check.names = F)
score_pb$group <- sample_info$group[match(rownames(score_pb), str_replace_all(sample_info$sample, '_', '-'))]
score_mean <- aggregate(. ~ group, data = score_pb, FUN = mean) |> column_to_rownames(var = 'group') |> apply(1, as.numeric)
rownames(score_mean) <- names(gs_list)
Heatmap(score_mean[gs_sub,] |> t() |> scale() |> t(), col = rev(brewer.pal(10, 'RdBu')), name = 'Enrichment')

expr_mat <- expr_mat[, sample_common]
ssgseaPar <- ssgseaParam(as.matrix(expr_mat), gs_list)
score_tcga <- gsva(ssgseaPar)

reference <- score_mean[gs_sub,] |> t() |> scale() |> t()
samples_by_cancer_type <- split(df_Bagaev$Sample, df_Bagaev$TCGA_project)
scaled_score_list <- lapply(samples_by_cancer_type, function(samples) {
  t(scale(t(score_tcga[gs_sub, samples])))
})
query <- do.call(cbind, scaled_score_list)

# Assignment
mtx_combined <- cbind(reference, query) 
M <- cor(mtx_combined, method = 'spearman')
mtx_raw_cor <- M[colnames(reference), colnames(query)] |>
  data.frame(check.names = F)
mtx_normalized_score <- t(mtx_raw_cor) |>
  scale() |>
  data.frame(check.names = F)
df_raw_cor_long <- mtx_raw_cor |>
  rownames_to_column(var = "Subtype") |>
  pivot_longer(
    cols = -Subtype,
    names_to = "sample",
    values_to = "Actual_Correlation"
  )
df_norm_score_long <- mtx_normalized_score |>
  rownames_to_column(var = "sample") |>
  pivot_longer(
    cols = -sample,
    names_to = "Subtype",
    values_to = "Normalized_Score"
  )
df_combined <- df_raw_cor_long |>
  inner_join(df_norm_score_long, by = c("sample", "Subtype"))
sample_assignment <- df_combined |>
  group_by(sample) |>
  slice_max(Actual_Correlation, n = 1, with_ties = FALSE) |>
  ungroup()  |>
  filter(Actual_Correlation > 0.25) 
sample_assignment |> 
  tabyl(Subtype)

sample_assignment <- merge(sample_assignment, df_Bagaev, by.x = 'sample', by.y = 'Sample')
immune_subtype <- immune_subtype[immune_subtype$sample %in% sample_assignment$sample]
sample_assignment <- merge(sample_assignment, immune_subtype, by = 'sample')
sample_assignment$Subtype_Immune_Model_Based <- factor(sample_assignment$Subtype_Immune_Model_Based, 
                                                       levels = c('Wound Healing (Immune C1)', 'IFN-gamma Dominant (Immune C2)', 'Inflammatory (Immune C3)', 
                                                                  'Lymphocyte Depleted (Immune C4)', 'Immunologically Quiet (Immune C5)', 'TGF-beta Dominant (Immune C6)'))
sample_assignment <- sample_assignment |> 
  mutate(MFP = case_when(MFP == 'D' ~ 'Depleted',
                         MFP == 'F' ~ 'Fibrotic',
                         MFP == 'IE' ~ 'Immune-enriched',
                         MFP == 'IE/F' ~ 'Immune-enriched, Fibrotic')) |> 
  mutate(MFP = factor(MFP, levels = c('Immune-enriched, Fibrotic', 'Immune-enriched', 'Fibrotic','Depleted'))) |> 
  mutate(Subtype = case_when(Subtype == 'Immune Quiescent' ~ 'TIME-Q',
                             Subtype == 'Immune Inflamed' ~ 'TIME-I',
                             Subtype == 'B Cell-enriched' ~ 'TIME-B',
                             Subtype == 'Myeloid-enriched' ~ 'TIME-Mye')) |> 
  mutate(Subtype = factor(Subtype, levels = c('TIME-Q','TIME-I','TIME-B','TIME-Mye')))
sample_assignment <- arrange(sample_assignment, Subtype, TCGA_project, MFP)
sample_assignment$ecotype <- df_ecotype$assignment[match(sample_assignment$sample, df_ecotype$sample)]

sample_assignment |>
  tabyl(TCGA_project, Subtype) |>
  pivot_longer(col = -TCGA_project, values_to = 'number', names_to = 'Subtype') |>
  ggplot(aes(x = TCGA_project, y = number, fill = Subtype)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = structure(names = as.character(unique(sample_assignment$Subtype)), met.brewer("Juarez", length(unique(sample_assignment$Subtype))))) +
  labs(y = "Relative Proportion") +
  theme_classic() +
  theme()

col_ha = HeatmapAnnotation(
  `Cancer Type` = sample_assignment$TCGA_project,
  # `Thorsson et al 2018` = sample_assignment$Subtype_Immune_Model_Based,
  # `Bagaev et al 2021` = sample_assignment$MFP,
  `TIME Subtype` = sample_assignment$Subtype,
  col = list(
    `TIME Subtype` = structure(names = as.character(unique(sample_assignment$Subtype)), met.brewer("Juarez", length(unique(sample_assignment$Subtype)))),
    `Cancer Type` = structure(names = as.character(unique(sample_assignment$TCGA_project)), met.brewer("Juarez", length(unique(sample_assignment$TCGA_project))))
  )
)# row split
celltype <- c(
  "T cells","Treg cells", "Effector T cells","NK cells", "Th1","Th2","Th17", 
  "B cells", "Plasma cells", "MHC-I", "MHC-II", "TLS formation", 
  "DC", "Pan-macrophage", "M1-like", "M2-like",
  "Co-activation molecules", "Checkpoint molecules", "Interferon response",
  "Endothelium", "Angiogenesis", "CAF", "Matrix remodeling","Stromal immiune exclusion","Stromal suppression"
)
celltype.group <- factor(c(
  rep("T/NK",7),
  rep("B/plasma",5),
  rep('myeloid',4),
  rep("Signature",3),
  rep("Stromal",6)), levels = c("T/NK", "B/plasma", 'myeloid', "Signature", "Stromal"))
df.row<- data.frame(group = celltype.group); rownames(df.row) <- celltype
pdf('figures/Fig3//ht_TCGA.pdf', height = 6, width = 12)
samples_by_cancer_type <- split(df_Bagaev$Sample, df_Bagaev$TCGA_project)
scaled_score_list <- lapply(samples_by_cancer_type, function(samples) {
  t(scale(t(score_tcga[gs_sub, samples])))
})
scores <- do.call(cbind, scaled_score_list) 

Heatmap(scores[,sample_assignment$sample], name = "Enrichment Score\n(z-score)", 
        column_title = paste0('TCGA pan-cancer (n=', nrow(sample_assignment),')'),
        column_title_gp = gpar(fontsize = 14),
        clustering_method_columns = 'ward.D2', clustering_method_rows = 'ward.D2',
        column_split = sample_assignment$Subtype,
        row_split = df.row,
        cluster_rows = F, cluster_columns = F,
        show_row_names = T, show_column_names = F, 
        show_row_dend = T, show_column_dend = T,
        col = circlize::colorRamp2(c(-2, 0, 2), c("#154999", "white", "#CF0034")),
        na_col = 'lightgray',
        heatmap_legend_param = list(
          title_position = "topcenter",
          at = c(-2, 0, 2), labels = c("Min", "", "Max"),
          legend_direction = "horizontal", border = T),
        top_annotation = col_ha,
        row_names_gp = gpar(fontsize = 12), use_raster=T)
dev.off()

df_os <- sample_assignment |> 
  dplyr::select(sample, Subtype,MFP, OS, OS_FLAG) |> 
  drop_na() |> 
  mutate(OS = OS/30)
df_os$Subtype <- factor(df_os$Subtype, levels = c('TIME-Q','TIME-I','TIME-B','TIME-Mye'))
fit <- survfit(Surv(OS, OS_FLAG) ~ Subtype, data = df_os)
p <- ggsurvplot(fit, data = df_os, conf.int = F, pval = T, 
                legend = "top",
                legend.title = "",
                legend.labs = c('TIME-Q','TIME-I','TIME-B','TIME-Mye'),
                surv.median.line = "hv",
                risk.table = T,
                risk.table.col = "strata",
                risk.table.height = 0.3,
                ggtheme = theme_classic(),
                palette =  met.brewer("Juarez",4),
                xlab = 'Overall Survival(month)', 
                xlim = c(0,144),
                break.time.by = 24) ; p

pdf('figures/Fig3/surv_os.pdf', height = 5, width = 4.5)
print(p, newpage = FALSE)
dev.off()
df_pfs <- sample_assignment |> 
  dplyr::select(sample, Subtype, MFP, PFS, PFS_Flag) |> 
  drop_na() |> 
  mutate(PFS = PFS/30)
df_pfs$Subtype <- factor(df_pfs$Subtype, levels = c('TIME-Q','TIME-I','TIME-B','TIME-Mye'))
fit <- survfit(Surv(PFS, PFS_Flag) ~ Subtype, data = df_pfs)
p <- ggsurvplot(fit, data = df_pfs, conf.int = F, pval = T, 
                legend = "top",
                legend.title = "",
                legend.labs = c('TIME-Q','TIME-I','TIME-B','TIME-Mye'),
                surv.median.line = "hv",
                risk.table = TRUE,
                risk.table.col = "strata",
                risk.table.height = 0.3,
                ggtheme = theme_classic(),
                palette =  met.brewer("Juarez",4),
                xlab = 'Prognosis-free Survival(month)'
);p
pdf('figures/Fig3/surv_pfs.pdf', height = 5, width = 4.5)
print(p, newpage = FALSE)
dev.off()

pdf('figures/Fig3/ht_compare_2021.pdf', height = 3, width = 5.5)
sample_assignment$MFP <- sample_assignment$MFP |> recode(Fibrotic = 'F', Depleted = 'D', `Immune-enriched, Fibrotic` = 'IE/F', `Immune-enriched` = 'IE') 
table(sample_assignment$MFP, sample_assignment$Subtype) |>
  t() |>
  scale() |>
  t() |>
  Heatmap(column_names_rot = 30, column_title = 'TIME Subtype', row_title = 'Bagaev et al 2021',
          col = brewer.pal(8, 'Blues'),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.0f", pindex(table(sample_assignment$MFP, sample_assignment$Subtype), i, j)), x, y, gp = gpar(fontsize = 10))
          },
          show_heatmap_legend = FALSE, border = T,
          width = 5*unit(8, "mm"), 
          height = 5*unit(8, "mm")
  )
dev.off()
pdf('figures/Fig3/ht_compare_2018.pdf', height = 4, width = 6.5)
table(sample_assignment$Subtype_Immune_Model_Based, sample_assignment$Subtype) |>
  t() |>
  scale() |>
  t() |>
  Heatmap(column_names_rot = 30, column_title = 'TIME Subtype', row_title = 'Thorsson et al 2018',
          col = brewer.pal(8, 'Blues'),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.0f", pindex(table(sample_assignment$Subtype_Immune_Model_Based, sample_assignment$Subtype), i, j)), x, y, gp = gpar(fontsize = 10))
          },
          show_heatmap_legend = FALSE, border = T,
          width = 4*unit(12, "mm"), 
          height = 6*unit(8, "mm")
  )
dev.off()

pdf('figures/Fig3/ht_compare_ecotype.pdf', height = 4.5, width = 4.5)
table(sample_assignment$ecotype, sample_assignment$Subtype) |>
  t() |>
  scale() |>
  t() |>
  Heatmap(column_names_rot = 30, column_title = 'TIME Subtype', row_title = 'Ecotype',
          col = brewer.pal(8, 'Blues'),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.0f", pindex(table(sample_assignment$ecotype, sample_assignment$Subtype), i, j)), x, y, gp = gpar(fontsize = 10))
          },
          show_heatmap_legend = FALSE, border = T,
          width = 4*unit(12, "mm"), 
          height = 12*unit(6, "mm"))
dev.off()














