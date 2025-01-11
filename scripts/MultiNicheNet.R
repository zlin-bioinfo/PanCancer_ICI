#!/usr/bin/env Rscript
rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','SingleCellExperiment','multinichenetr','dior','qs','scRNAtoolVis','clusterProfiler','org.Hs.eg.db','SignatuR','RColorBrewer')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

lr_network = readRDS("/bigdata/zlin/Melanoma_meta/data/lr_network_human/lr_network_human_21122021.rds")
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))
ligand_target_matrix = readRDS("/bigdata/zlin/Melanoma_meta/data/lr_network_human/ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()

datasets <- c('SKCM_Becker', 'BCC_Yost', 'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Shiao', 'TNBC_Zhang', 'HNSC_Franken', 'HNSC_IMCISION', 'HNSC_Luoma', 'CRC_Li', 'PCa_Hawley')
list_seu <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qread(paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu_r2.qs')) 
  seu$response_timepoint <- paste0(seu$response, '_', seu$time_point)
  return(seu)})
names(list_seu) <- datasets
list_seu[['BCC_Yost']] <- subset(list_seu[['BCC_Yost']], subset = patient %in% c("BCC_Yost_su009", "BCC_Yost_su012"), invert = T)
seu <- merge(x=list_seu[[1]], y=list_seu[2:length(list_seu)], add.cell.ids = names(list_seu))
seu <- subset(seu, subset = response == 'NE', invert = T)
seu$modality <- ifelse(seu$treatment == 'aPD1+CTLA4', 'Dual', 'Mono')
seu$interval <- as.numeric(seu$interval)
seu$int_cat <- ifelse(seu$interval >14, 'early', 'late')
genes_to_keep <- rownames(seu)[Matrix::rowSums(seu@assays$RNA$counts > 0) > 30]
hb <- setdiff(rownames(seu)[grepl("^HB[^P]", rownames(seu))], c("HBEGF", "HBS1L"))
ig <- rownames(seu)[grepl("^(IGH|IGK|IGL)", rownames(seu))]
tcr <- rownames(seu)[grepl("^(TRAV|TRBV|TRDV|TRGV)", rownames(seu))]
rm_genes <- c(hb, ig, tcr)
seu <- seu[setdiff(genes_to_keep, rm_genes), ]
seu <- JoinLayers(seu)
sce <- as.SingleCellExperiment(seu) %>% alias_to_symbol_SCE("human") %>% makenames_SCE()
sce$celltype_main <- sce$celltype_main %>% make.names()
rm(list_seu, seu)

sample_id = "sample"
group_id = "timepoint"
celltype_id = "celltype_main"
covariates = c("int_cat", "modality", "patient", "response")
batches = NA

senders_oi = unique(sce$celltype_main)
receivers_oi = unique(sce$celltype_main)

min_cells = 3

# Step1: Abundance and expression
sce$sample <- sce$sample %>% make.names()
abundance_expression_info = get_abundance_expression_info(sce = sce, sample_id = sample_id, group_id = group_id, 
                                                          celltype_id = celltype_id, min_cells = min_cells, senders_oi = senders_oi,
                                                          receivers_oi = receivers_oi, lr_network = lr_network, batches = batches)
Sys.time()
abundance_expression_info$abund_plot_sample

# Step2: Differential expression analysis
contrasts_oi = c("'(RE_Post-RE_Pre)-(NR_Post-NR_Pre)','(NR_Post-NR_Pre)-(RE_Post-RE_Pre)'")
contrast_tbl = tibble(contrast = c("(RE_Post-RE_Pre)-(NR_Post-NR_Pre)", "(NR_Post-NR_Pre)-(RE_Post-RE_Pre)"),
                      group = c("RE_Post","NR_Post")) 
DE_info = get_DE_info(sce = sce, sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, batches = batches, covariates = covariates, contrasts_oi = contrasts_oi, min_cells = min_cells)
# Check DE results
DE_info$celltype_de$de_output_tidy %>% arrange(p_val) %>% head()
DE_info$hist_pvals
# Empirical Null procedure
empirical_pval = TRUE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
}
DE_info_emp$hist_pvals_emp
# Compare p-values
comparison_plots = compare_normal_emp_pvals(DE_info, DE_info_emp, adj_pval = FALSE)
comparison_plots
# Normal/empirical p-values
if(empirical_pval == FALSE){
  celltype_de = DE_info$celltype_de$de_output_tidy
} else {
  celltype_de = DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>% dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
}
# Combine DE&LR
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
# Step3: Prediction
logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05
# p_val_adj = TRUE 
p_val_adj = FALSE 
empirical_pval = FALSE
top_n_target = 250
verbose = TRUE
cores_system = 20
n.cores = min(cores_system, union(senders_oi, receivers_oi) %>% length()) # use one core per receiver cell type
# Run NicheNet
ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(get_ligand_activities_targets_DEgenes(
  receiver_de = celltype_de,
  receivers_oi = receivers_oi,
  ligand_target_matrix = ligand_target_matrix,
  logFC_threshold = logFC_threshold,
  p_val_threshold = p_val_threshold,
  p_val_adj = p_val_adj,
  top_n_target = top_n_target,
  verbose = verbose, 
  n.cores = n.cores
)))

# Step4: Prioritizing all LR pairs
prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)

prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de = sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  contrast_tbl = contrast_tbl,
  sender_receiver_tbl = sender_receiver_tbl,
  grouping_tbl = grouping_tbl,
  prioritizing_weights = prioritizing_weights,
  fraction_cutoff = fraction_cutoff, 
  abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
  abundance_data_sender = abundance_expression_info$abundance_data_sender
))
# Step 5: Add info of prior knowledge and expression correlation between LR and target expression
lr_target_prior_cor = lr_target_prior_cor_inference(prioritization_tables$group_prioritization_tbl$receiver %>% unique(), abundance_expression_info, celltype_de, grouping_tbl, prioritization_tables, ligand_target_matrix, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold, p_val_adj = p_val_adj)

# Save ouput

path = "/bigdata/zlin/Melanoma_meta/data/"

multinichenet_output = list(
  celltype_info = abundance_expression_info$celltype_info,
  celltype_de = celltype_de,
  sender_receiver_info = abundance_expression_info$sender_receiver_info,
  sender_receiver_de =  sender_receiver_de,
  ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
  prioritization_tables = prioritization_tables,
  grouping_tbl = grouping_tbl,
  lr_target_prior_cor = lr_target_prior_cor
) 
multinichenet_output = make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output_", Sys.Date(), ".rds"))
}

multinichenet_output <- readRDS('/bigdata/zlin/Melanoma_meta/data/multinichenet_output.rds')

# Top 50 predictions
prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, rank_per_group = FALSE)
prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
# Top 30 predictions
prioritized_tbl_oi_RE_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = "RE_Post")
prioritized_tbl_oi_NR_30 = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 30, groups_oi = "NR_Post")
circos_E = make_circos_one_group(prioritized_tbl_oi_RE_30, colors_sender, colors_receiver)
circos_NE = make_circos_one_group(prioritized_tbl_oi_NR_30, colors_sender, colors_receiver)

plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_RE_30)
plot_oi

plot_oi = make_sample_lr_prod_activity_plots(multinichenet_output$prioritization_tables, prioritized_tbl_oi_NR_30)
plot_oi

# Volcano plot
res_de <- multinichenet_output$celltype_de 
colnames(res_de)[which(colnames(res_de) == 'logFC')] <- 'avg_log2FC'
colnames(res_de)[which(colnames(res_de) == 'p_adj')] <- 'p_val_adj'
colnames(res_de)[which(colnames(res_de) == 'cluster_id')] <- 'cluster'
res_de <- filter(res_de, contrast == "(RE_Post-RE_Pre)-(NR_Post-NR_Pre)")
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
jjVolcano(diffData = res_de,
          log2FC.cutoff = 0.5,
          col.type = "adjustP", adjustP.cutoff = 0.05, 
          topGeneN = 5, tile.col = getPalette(length(unique(res_de$cluster))))

ids=bitr(unique(res_de$gene),'SYMBOL','ENTREZID','org.Hs.eg.db') 
res_de=merge(res_de,ids,by.x='gene',by.y='SYMBOL')



