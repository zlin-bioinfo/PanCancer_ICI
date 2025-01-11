rm(list=ls())
pkgs <- c('GEOquery','clusterProfiler','org.Hs.eg.db','dplyr','tidyr','stringr','xCell','readxl','patchwork','ComplexHeatmap','RColorBrewer','ggpubr','ggplot2', 'ggsci','estimate')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

query <- 'GSE91061'
getGEOSuppFiles(query, makeDirectory = TRUE, baseDir = '/bigdata/zlin/Melanoma_meta/data/',
                fetch_files = TRUE, filter_regex = 'fpkm')
R.utils::gunzip('/bigdata/zlin/Melanoma_meta/data/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz')
df_expr <- read.csv('/bigdata/zlin/Melanoma_meta/data/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv')

# convert entrezid to gene symbol
e2s <- bitr(df_expr$X, fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Hs.eg.db')
df_expr <- filter(df_expr, X %in% e2s$ENTREZID) %>% tibble::column_to_rownames(var = 'X')
if (all.equal(rownames(df_expr),e2s$ENTREZID)){
  rownames(df_expr) <- e2s$SYMBOL
}

# Tumor purity
estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
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
TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])

# Gene expression
metadata <- data.frame(t(df_expr['CLEC9A',]))
metadata$purity <- TumorPurity
metadata$pt <- str_split(rownames(metadata), '_', simplify = T)[,1]
metadata$tp <- str_split(rownames(metadata), '_', simplify = T)[,2]
patient <- metadata |> 
  group_by(pt) |> 
  summarize(n = n()) |>
  filter(n == 2) |> 
  pull(pt)
metadata <- metadata |> filter(pt %in% patient)
df_clin <- read_xlsx('/bigdata/zlin/Melanoma_meta/data/GSE91061/1-s2.0-S0092867417311224-mmc2.xlsx', skip = 2) %>% tibble::column_to_rownames(var = 'Patient')
pt_nr <- intersect(rownames(df_clin %>% filter(Response %in% c('PD','SD'))), patient)
pt_r <- intersect(rownames(df_clin %>% filter(Response %in% c('PR','CR'))), patient)
metadata$response <- 'NE'
metadata$response[metadata$pt %in% pt_nr] <- 'NR'
metadata$response[metadata$pt %in% pt_r] <- 'RE'
metadata$tp <- factor(metadata$tp, levels = c('Pre', 'On'))
metadata$response <- factor(metadata$response, levels = c('RE', 'NR'))
metadata <- filter(metadata, pt != 'Pt109', response != 'NE')
metadata$prior <- df_clin$Cohort[match(metadata$pt, rownames(df_clin))]
metadata$prior <- df_clin$Cohort[match(metadata$pt, rownames(df_clin))]

metadata <- metadata |> group_by(pt) |> 
  mutate(diff = CLEC9A[tp == "On"] - CLEC9A[tp == "Pre"],
        purity_0.7  = all(purity > 0.7)) 
metadata |> group_by(pt) |> 
  filter(response != 'NE', diff > 0.2 | diff < -0.2, purity_0.7 == FALSE) |> 
  ggplot(aes(x = tp, y = CLEC9A)) +
  geom_boxplot(aes(fill = tp), outlier.shape = NA) +
  geom_line(aes(group = pt), color = "gray", linetype = 'dashed') +
  geom_point(aes(color = prior), size = 2, alpha = 0.8) +
  scale_color_jama() +
  scale_fill_aaas() +
  facet_wrap(~response, scales = "free") +
  stat_compare_means(paired = T, method = 't.test') +
  theme_classic2() + xlab("") + ylab("CLEC9A Expression") + ggtitle("GSE91061 (Melanoma)") + labs(color = "Prior ICB", fill = "Time") 

ggsave('/bigdata/zlin/Melanoma_meta/figures/CLEC9A.pdf')

df_comp <- xCellAnalysis(df_expr) %>% t() %>%  data.frame(check.names = F)
df_comp <- df_comp %>% 
  mutate(patient = str_split(rownames(.), '_', simplify = T)[,1]) %>% 
  mutate(time_point = str_split(rownames(.), '_', simplify = T)[,2])

pt_paired <- c()
for (pt in unique(df_comp$patient)){
  tp <- df_comp %>% filter(patient == pt) %>% dplyr::select(time_point)
  if ('Pre' %in% tp$time_point & 'On' %in% tp$time_point){
    pt_paired <- c(pt_paired, pt)
  }
}
df_clin <- read_xlsx('/bigdata/zlin/Melanoma_meta/data/GSE91061/1-s2.0-S0092867417311224-mmc2.xlsx', skip = 2) %>% tibble::column_to_rownames(var = 'Patient')
pt_nr <- intersect(rownames(df_clin %>% filter(Response %in% c('PD','SD'))), pt_paired)
pt_r <- intersect(rownames(df_clin %>% filter(Response %in% c('PR','CR'))), pt_paired)
pt_paired <- setdiff(pt_paired, intersect(rownames(df_clin %>% filter(Response %in% c('NE'))), pt_paired))

change_val <- function(df_comp, pt){
  df_pre <- df_comp %>% 
    filter(patient %in% pt, time_point == 'Pre') %>% 
    dplyr::select(!time_point) %>% 
    pivot_longer(cols = -patient, names_to = 'celltype', values_to = 'pre')
  df_on <- df_comp %>% 
    filter(patient %in% pt, time_point == 'On') %>% 
    dplyr::select(!time_point) %>% 
    pivot_longer(cols = -patient, names_to = 'celltype', values_to = 'on')
  val <- list()
  celltype <- colnames(df_comp)[-c(65:69)]
  for (i in 1:length(celltype)){
    df <- merge(df_pre[df_pre$celltype == celltype[i],], df_on[df_on$celltype == celltype[i],], by='patient', all=TRUE) 
    df <- replace(df, is.na(df), 0)
    t_test <- t.test(df$on, df$pre, paired = T, alternative = 'two.sided')
    val[[i]] <- c(t_test$p.value, t_test$statistic, p.adjust(t_test$p.value, method = 'fdr', n=length(celltype)))
  }
  mat <- do.call(rbind, val) %>% data.frame()
  colnames(mat) <- c('p_value','t_score','FDR')
  rownames(mat) <- colnames(df_comp)[-c(65:69)]
  return(mat)
}

mat_r <- change_val(df_comp, pt_r)
mat_nr <- change_val(df_comp, pt_nr)

celltype_sig <- union(rownames(filter(mat_r, `p_value` < 0.05)), rownames(filter(mat_nr, `p_value` < 0.05)))
# Heatmap
mat_val <- cbind('Non-responder' = mat_nr[celltype_sig, 't_score'], 
                 'Responder' = mat_r[celltype_sig, 't_score'])
rownames(mat_val) <- celltype_sig
mat_val[mat_val > 3] <- 3
mat_val[mat_val < -3] <- -3
mat_p <- cbind('Non-responder' = mat_nr[celltype_sig, 'p_value'], 
               'Responder' = mat_r[celltype_sig, 'p_value'])
rownames(mat_p) <- celltype_sig

# fdr_mat <- Reduce(full_join, p_list) %>% 
#   arrange(celltype) %>% 
#   column_to_rownames(var = 'celltype')
Heatmap(mat_val, name = "mat", na_col = "white", col = rev(brewer.pal(n=8, name = 'RdBu')),
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, row_names_side = "left", 
        heatmap_legend_param = list(title = "t-value"),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        width = ncol(mat_val)*unit(6, "mm"), 
        height = nrow(mat_val)*unit(5, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(!is.na(mat_p[i, j]) & mat_p[i, j] < 0.001) {
            grid.text('***', x, y,gp=gpar(fontsize=14))
          }
          else if(!is.na(mat_p[i, j]) & mat_p[i, j] < 0.01) {
            grid.text('**', x, y,gp=gpar(fontsize=14))
          }
          else if(!is.na(mat_p[i, j]) & mat_p[i, j] < 0.05) {
            grid.text('*', x, y,gp=gpar(fontsize=14))
          }
        })

# paired boxplot
check <- 'Tregs'
df_pre <- df_comp %>% 
  filter(patient %in% pt_nr, time_point == 'Pre') %>% 
  dplyr::select(!time_point) %>% 
  pivot_longer(cols = -patient, names_to = 'celltype', values_to = 'pre') %>% 
  filter(celltype == check)
df_on <- df_comp %>% 
  filter(patient %in% pt_nr, time_point == 'On') %>% 
  dplyr::select(!time_point) %>% 
  pivot_longer(cols = -patient, names_to = 'celltype', values_to = 'on') %>% 
  filter(celltype == check)
df <- merge(df_pre, df_on, by='patient', all=TRUE) 
df$response <- df_clin[df$patient, 'Response']
p_nr <- ggpaired(df, cond1 = "pre", cond2 = "on",fill = "white", line.color = "gray", line.size = 0.3, group = 'patient', title = 'Non-responder') +
  geom_point(aes(color = response)) + scale_color_aaas(name='Response') +
  stat_compare_means(paired = TRUE, method = 't.test', label.y = 0.075) + 
  theme(legend.direction = "vertical", legend.position = 'right') +
  xlab("Time Point") + ylab('Enrichment') 

df_pre <- df_comp %>% 
  filter(patient %in% pt_r, time_point == 'Pre') %>% 
  dplyr::select(!time_point) %>% 
  pivot_longer(cols = -patient, names_to = 'celltype', values_to = 'pre') %>% 
  filter(celltype == check)
df_on <- df_comp %>% 
  filter(patient %in% pt_r, time_point == 'On') %>% 
  dplyr::select(!time_point) %>% 
  pivot_longer(cols = -patient, names_to = 'celltype', values_to = 'on') %>% 
  filter(celltype == check)
df <- merge(df_pre, df_on, by='patient', all=TRUE) 
df$response <- df_clin[df$patient, 'Response']
p_r <- ggpaired(df, cond1 = "pre", cond2 = "on",fill = "white", line.color = "gray", line.size = 0.3, group = 'patient', title = 'Responder') +
  geom_point(aes(color = response)) + scale_color_aaas(name='Response') +
  stat_compare_means(paired = TRUE, method = 't.test', label.y = 0.075) + 
  theme(legend.direction = "vertical", legend.position = 'right') +
  xlab("Time Point") + ylab('Enrichment') 
p_nr + p_r + plot_annotation(title = check) + plot_layout(guides = 'collect')


df_pre <- df_comp %>% 
  filter(patient %in% pt_paired, time_point == 'Pre') %>% 
  dplyr::select(!time_point) %>% 
  pivot_longer(cols = -patient, names_to = 'celltype', values_to = 'pre') %>% 
  filter(celltype == check)
df_on <- df_comp %>% 
  filter(patient %in% pt_paired, time_point == 'On') %>% 
  dplyr::select(!time_point) %>% 
  pivot_longer(cols = -patient, names_to = 'celltype', values_to = 'on') %>% 
  filter(celltype == check)
df <- merge(df_pre, df_on, by='patient', all=TRUE) 
df$response <- df_clin[df$patient, 'Response']
df$res_cat[df$response %in% c('PD','SD')] <- 'Non-responder'
df$res_cat[df$response %in% c('PR','CR')] <- 'Responder'
ggpaired(df, cond1 = "pre", cond2 = "on",fill = "white", line.color = "gray", line.size = 0.3, group = 'patient', title = check) +
  geom_point(aes(color = response)) + scale_color_d3(name='Response') +
  stat_compare_means(paired = TRUE, method = 't.test', label.y = 0.07) + 
  facet_grid(.~res_cat) + theme_minimal() +
  theme(legend.direction = "vertical", legend.position = 'right', strip.text.x = element_text(size = 12)) +
  xlab("Time Point") + ylab('Enrichment') 








