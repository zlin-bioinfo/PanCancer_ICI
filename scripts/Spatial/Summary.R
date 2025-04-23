library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
clin_info <- readxl::read_xlsx('data/NSCLC_Yan/14728962/41588_2024_1998_MOESM3_ESM.xlsx',skip = 1)
clin_info <- clin_info |> drop_na(`ST \nbarcode`)

samples <- list.files('data/NSCLC_Yan/14728962/ST-visium/')
list_file <- list.files('data/NSCLC_Yan/MISTy/colocalization/')[list.files('data/NSCLC_Yan/MISTy/colocalization/') |> str_detect('interactions_RF')]
interaction_list <- lapply(samples, function(sample){
  interaction_sample <- read.csv(paste0('data/NSCLC_Yan/MISTy/colocalization/',sample,'_interactions_LM.csv'))
  interaction_sample$sample <- sample
  interaction_sample$response <- clin_info$`Pathologic\nresponse`[match(interaction_sample$sample, clin_info$`ST \nbarcode`)]
  interaction_sample$response <- ifelse(interaction_sample$response %in% c('MPR','pCR'), 'R', 'NR')
  return(interaction_sample)
}) 
interaction_aggregated <- do.call(rbind, interaction_list)
interaction_aggregated

interaction_aggregated |> 
  filter(view == 'intra') |> 
  group_by(target, predictor) |> 
  summarize(mean_importance = mean(importances, na.rm = TRUE), .groups = "drop") |> 
  pivot_wider(values_from = 'mean_importance', names_from = 'predictor') |> column_to_rownames(var = 'target') |> Heatmap(col = circlize::colorRamp2(c(-10, 0, 10), c("#154999", "white", "#CF0034")),
                                                                                                                          column_names_gp = gpar(fontsize = 7),
                                                                                                                          row_names_gp = gpar(fontsize = 8))
mat_nr <- interaction_aggregated |> 
  filter(view == 'intra', response == 'NR') |> 
  group_by(target, predictor) |> 
  summarize(mean_importance = mean(importances, na.rm = TRUE), .groups = "drop") |> 
  pivot_wider(values_from = 'mean_importance', names_from = 'predictor') |> column_to_rownames(var = 'target')
mat_r <- interaction_aggregated |> 
  filter(view == 'intra', response == 'R') |> 
  group_by(target, predictor) |> 
  summarize(mean_importance = mean(importances, na.rm = TRUE), .groups = "drop") |> 
  pivot_wider(values_from = 'mean_importance', names_from = 'predictor') |> column_to_rownames(var = 'target')
(mat_r-mat_nr) |> Heatmap(col = circlize::colorRamp2(c(-10, 0, 10), c("#154999", "white", "#CF0034")),
                          column_names_gp = gpar(fontsize = 7),
                          row_names_gp = gpar(fontsize = 10))


inter_mat_list <- lapply(samples, function(sample){
  interaction_sample <- read.csv(paste0('data/NSCLC_Yan/MISTy/colocalization/',sample,'_interactions_RF.csv'))
  interaction_sample$sample <- sample
  return(interaction_sample)
}) 

a <- read.csv(paste0('data/NSCLC_Yan/MISTy/colocalization/',sample,'_interactions_RF.csv'))
a <- a |> filter(view == 'para')
a <- a[,-which(colnames(a)=='view')]
a <- pivot_wider(a, values_from = 'importances', names_from = 'predictor') |> column_to_rownames(var = 'target')
pheatmap::pheatmap(a)


samples <- c("HCC1R",  "HCC2R", "HCC3R",  "HCC4R", "HCC5NR",  "HCC6NR", "HCC7NR")
list_file <- list.files('data/GSE238264_RAW/MISTy/colocalization/')[list.files('data/GSE238264_RAW/MISTy/colocalization/') |> str_detect('interactions_RF')]
interaction_list <- lapply(samples, function(sample){
  interaction_sample <- read.csv(paste0('data/GSE238264_RAW/MISTy/colocalization/',sample,'_interactions_LM.csv'))
  interaction_sample$sample <- sample
  return(interaction_sample)
}) 
interaction_aggregated <- do.call(rbind, interaction_list)
interaction_aggregated$response <- 'R'
interaction_aggregated$response[str_detect(interaction_aggregated$sample, 'NR')] <- 'NR'

mat_nr <- interaction_aggregated |> 
  filter(view == 'intra', response == 'NR') |> 
  group_by(target, predictor) |> 
  summarize(mean_importance = mean(importances, na.rm = TRUE), .groups = "drop") |> 
  pivot_wider(values_from = 'mean_importance', names_from = 'predictor') |> column_to_rownames(var = 'target')
mat_r <- interaction_aggregated |> 
  filter(view == 'intra', response == 'R') |> 
  group_by(target, predictor) |> 
  summarize(mean_importance = mean(importances, na.rm = TRUE), .groups = "drop") |> 
  pivot_wider(values_from = 'mean_importance', names_from = 'predictor') |> column_to_rownames(var = 'target')
(mat_r-mat_nr) |> pheatmap::pheatmap()

interaction_aggregated |> 
  filter(view == 'intra') |> 
  group_by(target, predictor) |> 
  summarize(mean_importance = mean(importances, na.rm = TRUE), .groups = "drop") |> 
  pivot_wider(values_from = 'mean_importance', names_from = 'predictor') |> column_to_rownames(var = 'target') |> 
  pheatmap::pheatmap(scale = 'row')



df <- read.csv('tables/cosg_genes.csv', check.names = F)

