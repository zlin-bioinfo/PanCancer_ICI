pkgs <- c('qs2','tidyr','dplyr','plyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','tibble','pheatmap','MetBrewer','viridis','ComplexHeatmap','colorRamp2','corrr','ggnewscale','NMF','ggpubr','corrplot','rstatix','janitor','Seurat')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

metadata <- read.csv('tables/meta_all.csv') 
metadata$freq_r2_comp[metadata$celltype_r2 == 'Malignant(CNA+)'] <- metadata$freq_r2[metadata$celltype_r2 == 'Malignant(CNA+)']
source('scripts/Celltype_classification.R')
mtx_immune <- metadata |> 
  filter(celltype_r2 %in% immune, 
         count_immune >= 50, non_malignant_count >= 100,
         subset %in% c('TME', 'CD45+sorted'),
         cohort != 'SKCM_this study') |> 
  select(sample,  freq_r2_comp, celltype_r2) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = celltype_r2, values_fill = 0) |> 
  column_to_rownames(var = 'sample') 

# Sample classification
mtx <- mtx_immune
scale_mtx <- apply(mtx, MARGIN = 2, function(x) (x - min(x)) / (max(x) - min(x)))
scale_mtx <- t(scale_mtx)
res.nmf <- nmf(scale_mtx, rank = 3:10, nrun = 100, seed = 1234, method = 'lee')
pdf('figures/NMF/rank_survey_2707.pdf', height = 6, width =8)
plot(res.nmf)
dev.off()

nrank=4
nmf.res = nmf(scale_mtx, rank=nrank, nrun=100, seed=1234, method = 'lee')
qs_save(nmf.res,'data/res_nmf_rank4_seed_1234.qs2')
nmf.res <- qs_read('data/res_nmf_rank4_seed_1234.qs2')

w <- basis(nmf.res)
colnames(w) = c('TIME-Mye','TIME-Q','TIME-I', 'TIME-B')
w_t = w |> 
  t() |>  
  scale()
# 
# # Initialize a list to store cell subtypes assigned to each NMF program,
# score_threshold <- 0
# assigned_and_scored_subtypes <- setNames(lapply(rownames(w_t), function(x) numeric(0)), rownames(w_t))
# for (j in 1:ncol(w_t)) {
#   cell_subtype_name <- colnames(w_t)[j]
#   coefficients_for_subtype <- w_t[, cell_subtype_name]
#   assigned_nmf_index <- which.max(coefficients_for_subtype)
#   assigned_nmf_name <- rownames(w_t)[assigned_nmf_index]
#   value_in_assigned_nmf <- coefficients_for_subtype[assigned_nmf_index]
#   other_nmf_values_for_subtype <- coefficients_for_subtype[names(coefficients_for_subtype) != assigned_nmf_name]
#   avg_other_nmfs <- mean(other_nmf_values_for_subtype)
#   score_difference <- value_in_assigned_nmf - avg_other_nmfs
#   if (score_difference >= score_threshold) {
#     assigned_and_scored_subtypes[[assigned_nmf_name]] <- c(
#       assigned_and_scored_subtypes[[assigned_nmf_name]],
#       setNames(score_difference, cell_subtype_name)
#     )
#   }
# }
# 
# final_prioritized_results <- lapply(assigned_and_scored_subtypes, function(x) {
#   if (length(x) > 0) {
#     sort(x, decreasing = TRUE)
#   } else {
#     x 
#   }
# })
# print(final_prioritized_results)
# qs_save(final_prioritized_results, 'data/prioritized_subtype.qs2')

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
pdf('figures/NMF/ht_sample_rank4.pdf', height = 5.5, width = 12)
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

# DE analysis on pseudobulk
gs_list <- readxl::read_xlsx('tables/Signature.xlsx') |>
  lapply(function(x){return(x[!is.na(x)])})
sample_info <- read.csv('tables/sample_info_time_updated.csv')
datasets <- sample_info |> dplyr::filter(subset == 'TME') |>  dplyr::select(cohort) |> pull() |> unique()
datasets[which(datasets == "BCC&SCC_Yost")] <- 'BCC_Yost'
sample_included <- sample_info |> pull(sample)
list_pseudobulk <- lapply(datasets, function(dataset) {
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |>
    subset(sample %in% sample_included & celltype_r2 %in% c(immune, nonimmune)) |> 
    NormalizeData() |>
    AverageExpression(group.by = 'sample', layer = 'data', return.seurat = T)
  seu$group <- sample_info$group[match(seu$sample, str_replace_all(sample_info$sample, '_', '-'))]
  return(seu)
})
seu <- merge(x = list_pseudobulk[[1]], y=list_pseudobulk[2:length(list_pseudobulk)])
seu <- JoinLayers(seu)
pseudobulk <- GetAssayData(seu, layer = 'data')
sample_info$sample <- str_replace_all(sample_info$sample, '_', '-')
pseudobulk <- pseudobulk[,intersect(sample_info$sample, colnames(pseudobulk))]
qs_save(pseudobulk, 'data/pseudobulk_TME.qs2')

# Recovery on in-house dataset
# This function creates and returns another function that nmf will call.
.custom_nmf_seeding <- function(fixed_W_matrix, data_matrix) {
  n_components <- ncol(fixed_W_matrix)
  # Initialize H randomly with a fixed seed for reproducibility of H's starting point
  dummy_H <- rnmf(n_components, data_matrix)@H
  # The actual seeding function to be passed to nmf
  function(model, target){
    basis(model) <- fixed_W_matrix # Set the basis matrix (W) to the pre-computed one
    coef(model) <- dummy_H        # Initialize the coefficient matrix (H) randomly
    return(model)
  }
}
# Helper function for custom NMF update logic.
# This function creates and returns another function that NMFStrategy will call.
# It only updates H, keeping W fixed, using the Euclidean update rule.
.custom_nmf_h_update <- function() {
  # The actual update function to be passed to NMFStrategy
  function (i, v, x, copy = FALSE, eps = .Machine$double.eps, ...) {
    w <- .basis(x) # Access the current (fixed) W matrix from the NMF model object
    h <- .coef(x)  # Access the current H matrix from the NMF model object
    
    # Use the internal Euclidean update rule for H
    h <- NMF:::std.euclidean.update.h(v, w, h,
                                      nbterms = NMF:::nbterms(x), # Number of basis terms
                                      ncterms = NMF:::ncterms(x), # Number of coefficient terms
                                      copy = copy)
    
    # Apply numerical stability constraints if necessary
    if (i %% 10 == 0) { # Perform this check periodically, e.g., every 10 iterations
      h <- pmax.inplace(h, eps, NMF:::icterms(x))
    }
    
    # If 'copy' is true, update the H slot of the NMF model object
    if (copy) {
      .coef(x) <- h
    }
    return(x) # Return the updated NMF model object
  }
}

#' Recover cellular modules in scRNA_seq data
#' @description Recover cellular module abundances in scRNA-seq data using a pre-fitted NMF model
#'              that was trained with the 'lee' method (Euclidean objective).
#' @param nmf_result_obj The fitted NMF result object (output from NMF::nmf with method = "lee").
#' @param mat_fq Normalized sample-cell frequency matrix for the new data.
#'
#' @return A matrix containing the abundance of the cellular modules in each sample.
#' @export
#'
#' @examples
#' # First, load the necessary library
#' # library(NMF)
#' #
#' # # Create dummy data for example: 50 features (celltypes), 20 samples
#' # V_train <- abs(matrix(rnorm(50*20), 50, 20))
#' # rownames(V_train) <- paste0("subCluster", 1:50)
#' # colnames(V_train) <- paste0("Sample_Train", 1:20)
#' #
#' # # Fit an NMF model with the 'lee' method
#' # my_nmf_result <- nmf(V_train, rank = 5, method = "lee", .opt='v')
#' #
#' # # Create new frequency data (e.g., from 15 new samples)
#' # new_mat_fq <- abs(matrix(rnorm(50*15), 50, 15))
#' # rownames(new_mat_fq) <- paste0("subCluster", 1:50)
#' # colnames(new_mat_fq) <- paste0("Sample_New", 1:15)
#' #
#' # # Recover cellular module abundances
#' # cm_abundances <- sc_recover(my_nmf_result, new_mat_fq)
#' # print(head(cm_abundances))
#' # print(dim(cm_abundances)) # Should be K x Number_of_New_Samples
sc_recover<-function(nmf_result_obj, mat_fq, cm_names){
  # --- 1. Extract and prepare the fixed W matrix (Basis Matrix) ---
  W_fixed <- basis(nmf_result_obj) |> as.data.frame()
  colnames(W_fixed) <- cm_names
  W_fixed_matrix <- as.matrix(W_fixed) # Convert to matrix
  # --- 2. Prepare the new frequency data matrix (Target Matrix) ---
  to_predict_matrix <- as.matrix(mat_fq)
  # --- 3. Filter common celltypes and remove zero-variance rows for robust calculation ---
  common_celltypes <- intersect(rownames(W_fixed_matrix), rownames(to_predict_matrix))
  if(length(common_celltypes) == 0){
    stop("Error: No common celltypes found between the NMF basis matrix (W) and the new frequency data (mat_fq). Please ensure row names (cell subset identifiers) match.")
  }
  
  W_fixed_filtered <- W_fixed_matrix[common_celltypes, , drop = FALSE]
  to_predict_filtered <- to_predict_matrix[common_celltypes, , drop = FALSE]
  # Check for and handle celltypes with zero variance in the new data
  row_variances <- apply(to_predict_filtered, 1, var) 
  row_variances[is.na(row_variances)] <- 0
  valid_rows <- row_variances > 0
  if(sum(valid_rows) == 0){
    warning("Warning: All matched celltypes in the new data have zero variance. Cannot perform NMF prediction. Returning a matrix of zeros.")
    # Return a zero matrix with correct dimensions and naming for the output
    res_h <- matrix(0, nrow = ncol(W_fixed_matrix), ncol = ncol(to_predict_matrix))
    colnames(res_h) <- colnames(to_predict_matrix)
    rownames(res_h) <- colnames(W_fixed_matrix) # Assign CM names
    return(res_h)
  }
  W_fixed_filtered <- W_fixed_filtered[valid_rows, , drop = FALSE]
  to_predict_filtered <- to_predict_filtered[valid_rows, , drop = FALSE]
  # --- 4. Define NMF prediction strategy (fixed W, Euclidean objective) ---
  n_components <- ncol(W_fixed_filtered)
  # Create the custom NMF strategy
  nmf_strategy <- NMFStrategy(
    'fixed_W_prediction',     # A unique name for this custom strategy
    method = "lee",           # We specify "lee" as the underlying method
    Update = .custom_nmf_h_update(), # Use our custom H-only update function
    objective = 'euclidean',  # Specifies the objective function to minimize
    Stop = 'connectivity'     # A stopping criterion for the optimization
  )
  # --- 5. Perform NMF prediction using the fixed W and custom strategy ---
  new_nmf_model <- nmf(
    x = to_predict_filtered,       # The data to be factorized (new frequencies)
    rank = n_components,           # The number of components (K) to predict
    nrun = 1,                      # Only one run is needed as W is fixed and H is initialized
    method = nmf_strategy,         # Our custom NMF strategy
    seed = .custom_nmf_seeding(W_fixed_filtered, to_predict_filtered), # Our custom seeding function
    .opt = 'P1'                    # Option to suppress progress messages
  )
  # --- 6. Extract and format the predicted H matrix (Coefficient Matrix) ---
  h <- coef(new_nmf_model) # Extract the H matrix from the resulting NMF model
  # Assign CM names (from W_fixed) to the rows of H for clarity
  rownames(h) <- colnames(W_fixed_matrix)
  # Ensure the order of CMs in the output H matrix is consistent with W_fixed_matrix
  h <- h[cm_names, , drop = FALSE]
  return(h)
}

source('scripts/Celltype_classification.R')
metadata <- read.csv('tables/meta_all.csv') 
metadata$freq_r2_comp[metadata$celltype_r2 == 'Malignant(CNA+)'] <- metadata$freq_r2[metadata$celltype_r2 == 'Malignant(CNA+)']
sample_info <- read.csv('tables/sample_info_time_updated.csv')
nmf.res <- qs_read('data/res_nmf_rank4_seed_1234.qs2')
celltype_included <- basis(nmf.res) |> rownames()
celltype_included <- celltype_included[-which(celltype_included %in% c('CD4_Tcm','CD8_Tem-early'))]
sample_included <- metadata |> 
  filter(cohort == 'SKCM_this study', celltype_r2 %in% celltype_included, 
         count_immune >= 50, non_malignant_count >= 100) |> 
  pull(sample) |> unique()
mat_fq <- metadata |> 
  filter(celltype_r2 %in% celltype_included, 
         count_immune >= 50, non_malignant_count >= 100) |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  select(sample, celltype_r2, freq_r2_comp) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = sample, values_fill = 0) |> 
  tibble::column_to_rownames(var = 'celltype_r2')
mat_fq <- apply(mat_fq, MARGIN = 1, function(x) (x - min(x)) / (max(x) - min(x))) |> t()
mat_cm <- sc_recover(nmf.res, mat_fq[, sample_included], 
                     cm_names = c('TIME-Mye','TIME-Q','TIME-I', 'TIME-B' ))
sample_assignment <- apply(mat_cm, 2, function(x) rownames(mat_cm)[which.max(x)]) |> data.frame(check.names = F)
colnames(sample_assignment) <- 'Subtype'
sample_assignment <- rownames_to_column(sample_assignment, var = 'sample')
sample_assignment$response <- metadata$response[match(sample_assignment$sample, metadata$sample)]
sample_assignment$tx_status <- metadata$tx_status[match(sample_assignment$sample, metadata$sample)]
sample_assignment |> tabyl(Subtype, response, tx_status)
sample_assignment$Subtype <- factor(sample_assignment$Subtype, levels = c('TIME-Q', 'TIME-I', 'TIME-B', 'TIME-Mye'))
sample_assignment$group <- sample_assignment$Subtype
sample_assignment$cohort <- 'SKCM_this study'
sample_assignment$patient <- metadata$patient[match(sample_assignment$sample, metadata$sample)]
pt <- sample_assignment |>
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |> 
  pull(patient) 
unmatched_pt <- c("SKCM_this study_Patient2", "SKCM_this study_Patient3")
sample_assignment$paired <- ifelse(sample_assignment$patient %in% setdiff(pt, unmatched_pt), 'Yes', 'No')
write.csv(sample_assignment, 'tables/sample_info_time_validation.csv', row.names = F)




