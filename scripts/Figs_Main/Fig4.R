pkgs <- c('qs2','tidyr','dplyr','plyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','tibble','MetBrewer','viridis','rstatix','janitor','ggalluvial')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
metadata <- read.csv('tables/meta_all.csv') 
metadata$freq_r2_comp[metadata$celltype_r2 == 'Malignant(CNA+)'] <- metadata$freq_r2[metadata$celltype_r2 == 'Malignant(CNA+)']
source('scripts/Celltype_classification.R')
sample_info <- read.csv('tables/sample_info_time_updated.csv')
sample_info |>
  mutate(group = case_when(group == 'Immune Quiescent' ~ 'TIME-Q',
                           group == 'Immune Inflamed' ~ 'TIME-I',
                           group == 'B Cell-enriched' ~ 'TIME-B',
                           group == 'Myeloid-enriched' ~ 'TIME-Mye')) |> 
  mutate(group = factor(group, levels = c('TIME-Q','TIME-I','TIME-B','TIME-Mye'))) |> 
  filter(paired == 'Yes') |>
  select(patient, tx_status, group, response, subtype, subtype) |> 
  pivot_wider(values_from = 'group', names_from = 'tx_status') |>
  arrange(response, subtype) |>
  mutate(row_value = 1) |>
  ggplot(aes(y = row_value, axis1 = Baseline, axis2 = Treated)) + # , axis3 = response
  geom_alluvium(aes(fill = Baseline), curve_type = "cubic", alpha=0.8) +
  geom_stratum(alpha = 0.95, width = 0.4) +
  geom_text(aes(label = after_stat(stratum)), stat = 'stratum', size = 5.5) +
  scale_x_discrete(limits = c("Baseline", "Treated"), expand = c(0.2, 0.2)) +
  scale_fill_manual(values = structure(names = levels(sample_info$group), met.brewer('Juarez',5))) +
  theme_void() +
  ggtitle("") + ylab("") + xlab("") + labs(fill = 'Cancer Type') +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black")) + theme(legend.position = "none")
ggsave('figures/Fig4/alluviu_bl_tx.pdf', height = 4.5, width = 6)

# Shift
pt_paired <- sample_info |>
  filter(paired == 'Yes') |>
  select(patient, tx_status, group, response, subtype, interval) |> 
  pivot_wider(values_from = 'group', names_from = 'tx_status') |> 
  mutate(dynamics = ifelse(Baseline == Treated, 'Stable', 'Shifted'),
         int_cat = ifelse(interval <= 14, '<=2w',
                          ifelse(interval <= 28, '2-4w', 
                                 ifelse(interval <= 98, '4-14w', '>8w'))))
test_result <- pt_paired |> filter(!is.na(int_cat), int_cat %in% c('<=2w','2-4w', '4-14w')) |> tabyl(dynamics, int_cat) |> chisq.test()
pt_paired |> filter(!is.na(int_cat), int_cat!= '>8w') |> 
  group_by(int_cat, dynamics) |>
  tally(name = "n") |> 
  group_by(int_cat) |>
  mutate(prop = n / sum(n)) |>
  ungroup() |> 
  ggplot(aes(x = int_cat, y = prop, fill = dynamics)) +
  geom_col(position = "fill", color='black', width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(
    aes(label = n), 
    position = position_stack(vjust = 0.5), 
    size = 4, 
    color = "black"
  ) +
  labs(
    x = "",
    y = "",
    fill = "Dynamics",
    title = ""
  ) +
  theme_classic2() + scale_fill_manual(values = c('Shifted' = '#1B9E77', 'Stable' = '#D95F02')) +
  theme(axis.text.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        plot.title = element_text(hjust = 0.5)) 
ggsave('figures/Fig4/barplot_dyn_interval.pdf', height = 4, width = 5)

# Response (baseline)
sample_info |> 
  filter(cohort != 'HCC_Guo') |> 
  mutate(group = case_when(group == 'Immune Quiescent' ~ 'TIME-Q',
                           group == 'Immune Inflamed' ~ 'TIME-I',
                           group == 'B Cell-enriched' ~ 'TIME-B',
                           group == 'Myeloid-enriched' ~ 'TIME-Mye')) |> 
  mutate(group = factor(group, levels = c('TIME-Q','TIME-I','TIME-B','TIME-Mye'))) |> 
  # mutate(group = factor(group, levels = c('Immune Quiescent', 'Immune Inflamed', 'B Cell-enriched', 'Myeloid-enriched'))) |> 
  select(patient, tx_status, group, response) |> 
  mutate(response = factor(response, levels = rev(c('R','NE','NR')))) |>
  # filter(response %in% c('R','NR')) |> 
  group_by(group, tx_status, response) |>
  tally(name = "n") |> 
  group_by(group, tx_status) |>
  mutate(prop = n / sum(n)) |>
  ggplot(aes(x = tx_status, y = prop, fill = response)) +
  geom_col(position = "fill", color='black', width = 0.8) 
  facet_wrap(~ group, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(
    aes(label = n), 
    position = position_stack(vjust = 0.5), 
    size = 4, 
    color = "black"
  ) +
  labs(
    x = "",
    y = "",
    fill = "Response",
    title = ""
  ) +
  theme_classic2() + scale_fill_manual(values = c('R' = '#CC0C00FF','NR' = '#5C88DAFF','NE' = '#84BD00FF')) +
  theme(axis.text.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black")) 
ggsave('figures/Fig4/barplot_tx_res.pdf', height = 4.5, width = 7)

# Signature score
pseudobulk <- qs_read('data/pseudobulk_TME.qs2')
datasets <- sample_info |> dplyr::filter(subset == 'TME') |>  dplyr::select(cohort) |> pull() |> unique()
datasets[which(datasets == "BCC&SCC_Yost")] <- 'BCC_Yost'
sample_included <- sample_info |> pull(sample)
list_pseudobulk <- lapply(datasets, function(dataset) {
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |>
    subset(sample %in% sample_included) |> 
    NormalizeData() |>
    AverageExpression(group.by = 'sample', layer = 'data', return.seurat = T)
  seu$sample <- colnames(seu)
  seu$group <- sample_info$group[match(seu$sample, str_replace_all(sample_info$sample, '_', '-'))]
  return(seu)
})
seu <- merge(x = list_pseudobulk[[1]], y=list_pseudobulk[2:length(list_pseudobulk)])
seu <- JoinLayers(seu)

# Others to Myeloid-enriched
gs_time <- qs_read('data/sig_time.qs2')
gs_list <- readxl::read_xlsx('tables/Signature.xlsx') |>
  lapply(function(x){return(x[!is.na(x)])})
gs_list <- c(gs_time, gs_list)
ssgseaPar <- ssgseaParam(pseudobulk, gs_list, normalize = T)
score <- gsva(ssgseaPar)|> t() |> data.frame(check.names = F)
pt <- sample_info |>
  filter(paired == "Yes", subset == 'TME') |> 
  dplyr::select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Treated == 'Myeloid-enriched') |> pull(patient)
pt_mye_bt <- sample_info |>
  filter(patient %in% pt) |> 
  dplyr::select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Baseline  == 'Myeloid-enriched', Treated == 'Myeloid-enriched') |> pull(patient)
sample_included <- sample_info |> filter(patient %in% setdiff(pt, pt_mye_bt)) |> pull(sample) |> str_replace_all('_','-')
score <- score[sample_included,]
score$group <- sample_info$group[match(rownames(score), str_replace_all(sample_info$sample, '_','-'))]
score$tx_status <- sample_info$tx_status[match(rownames(score), str_replace_all(sample_info$sample, '_','-'))]
score$patient <- sample_info$patient[match(rownames(score), str_replace_all(sample_info$sample, '_','-'))]
p1 <- score |> select(`Myeloid-enriched`, tx_status, patient) |>
  arrange(patient) |> 
  ggplot(aes(y=`Myeloid-enriched`, x = tx_status)) +
  geom_boxplot(aes(fill = tx_status), width=0.5, outliers = F, show.legend = F) + 
  geom_line(aes(group = patient), color = "lightgray",linetype = "dashed", size= 0.4) +
  geom_point(size = 0.3) +
  theme_classic() + ylab('') + xlab('') + ggtitle('Myeloid-enriched') +
  theme(axis.text.x = element_text(color = 'black', size =10),
        strip.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Baseline" = "Baseline\n(Not Mye)",
                              "Treated" = "Treated\n(Mye)")) +
  scale_fill_manual(values = c('Baseline' = '#637b31', 'Treated' = '#cf5e4e'), name = '') +
  stat_compare_means(method = 'wilcox', paired = T, label='p.signif', label.x = 1.4, size = 7, vjust=1); p1

p2 <- score |> select(`M2-like`, tx_status, patient) |> 
  arrange(patient) |> 
  ggplot(aes(y=`M2-like`, x = tx_status)) +
  geom_boxplot(aes(fill = tx_status), width=0.5, outliers = F, show.legend = F) + 
  geom_line(aes(group = patient), color = "lightgray",linetype = "dashed", size= 0.4) +
  geom_point(size = 0.3) +
  theme_classic() + ylab('') + xlab('') + ggtitle('M2-like TAM') +
  theme(axis.text.x = element_text(color = 'black', size =10),
        strip.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Baseline" = "Baseline\n(Not Mye)",
                              "Treated" = "Treated\n(Mye)")) +
  scale_fill_manual(values = c('Baseline' = '#637b31', 'Treated' = '#cf5e4e'), name = '') +
  stat_compare_means(method = 'wilcox', paired = T, label='p.signif', label.x = 1.4, size = 7, vjust=1); p2

p3 <- score |> select(`Phagocytosis inhibition`, tx_status, patient) |> 
  arrange(patient) |> 
  ggplot(aes(y=`Phagocytosis inhibition`, x = tx_status)) +
  geom_boxplot(aes(fill = tx_status), width=0.5, outliers = F, show.legend = F) + 
  geom_line(aes(group = patient), color = "lightgray",linetype = "dashed", size= 0.4) +
  geom_point(size = 0.3) +
  theme_classic() +
  ylab('') + xlab('') + ggtitle('Phagocytosis Inhibition') + 
  theme(axis.text.x = element_text(color = 'black', size =10),
        strip.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Baseline" = "Baseline\n(Not Mye)",
                              "Treated" = "Treated\n(Mye)")) +
  scale_fill_manual(values = c('Baseline' = '#637b31', 'Treated' = '#cf5e4e'), name = '') +
  stat_compare_means(method = 'wilcox', paired = T, label='p.signif', label.x = 1.4, size = 7, vjust=1); p3

p4 <- score |> select(`Matrix remodeling`, tx_status, patient) |> 
  arrange(patient) |> 
  ggplot(aes(y=`Matrix remodeling`, x = tx_status)) +
  geom_boxplot(aes(fill = tx_status), width=0.5, outliers = F, show.legend = F) + 
  geom_line(aes(group = patient), color = "lightgray",linetype = "dashed", size= 0.4) +
  geom_point(size = 0.3) +
  theme_classic() + 
  ylab('') + xlab('') + ggtitle('Matrix Remodeling') + 
  theme(axis.text.x = element_text(color = 'black', size =10),
        strip.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Baseline" = "Baseline\n(Not Mye)",
                              "Treated" = "Treated\n(Mye)")) +
  scale_fill_manual(values = c('Baseline' = '#637b31', 'Treated' = '#cf5e4e'), name = '') +
  stat_compare_means(method = 'wilcox', paired = T, label='p.signif', label.x = 1.5, size = 7, vjust=1); p4
p1+p2+p3+p4+plot_layout(nrow=1)
ggsave('figures/Dynamics/boxplot_mye_score.pdf', height = 3.5, width = 9)

# Others to Immune-inflamed
ssgseaPar <- ssgseaParam(pseudobulk, gs_list, normalize = T)
score <- gsva(ssgseaPar)|> t() |> data.frame(check.names = F)
pt <- sample_info |>
  filter(paired == "Yes", subset == 'TME') |> 
  dplyr::select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Treated == 'Immune Inflamed') |> pull(patient)
pt_i_bt <- sample_info |>
  filter(patient %in% pt) |> 
  dplyr::select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Baseline  == 'Immune Inflamed', Treated == 'Immune Inflamed') |> pull(patient)
sample_included <- sample_info |> filter(patient %in% setdiff(pt, pt_i_bt)) |> pull(sample) |> str_replace_all('_','-')
score <- score[sample_included,]
score$group <- sample_info$group[match(rownames(score), str_replace_all(sample_info$sample, '_','-'))]
score$tx_status <- sample_info$tx_status[match(rownames(score), str_replace_all(sample_info$sample, '_','-'))]
score$patient <- sample_info$patient[match(rownames(score), str_replace_all(sample_info$sample, '_','-'))]
score$response <- ifelse(score$patient %in% unique(sample_info$patient[sample_info$response == 'R']), 'R', 'NR') |> factor(levels = c('R','NR'))
p1 <- score |> select(`Immune Inflamed`, tx_status, patient, response) |>
  arrange(patient) |> 
  ggplot(aes(y=`Immune Inflamed`, x = tx_status)) +
  geom_boxplot(aes(fill = tx_status), width=0.5, outliers = F, show.legend = F) + 
  geom_line(aes(group = patient), color = "lightgray",linetype = "dashed", size= 0.4) +
  geom_point(size = 0.3) +
  theme_classic() + 
  ylab('') + xlab('') + ggtitle('Immune Inflamed') + 
  theme(axis.text.x = element_text(color = 'black', size =10),
        strip.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Baseline" = "Baseline\n(Not Inf)",
                              "Treated" = "Treated\n(Inf)")) +
  scale_fill_manual(values = c('Baseline' = '#637b31', 'Treated' = '#208CC0'), name = '') +
  facet_wrap(.~ response) +
  stat_compare_means(method = 'wilcox', paired = T, label='p.signif', label.x = 1.5, size = 7, vjust=1); p1

p2 <- score |> select(`Checkpoint molecules`, tx_status, patient, response) |> 
  arrange(patient) |> 
  ggplot(aes(y=`Checkpoint molecules`, x = tx_status)) +
  geom_boxplot(aes(fill = tx_status), width=0.5, outliers = F, show.legend = F) + 
  geom_line(aes(group = patient), color = "lightgray",linetype = "dashed", size= 0.4) +
  geom_point(size = 0.3) +
  theme_classic() + 
  ylab('') + xlab('') + ggtitle('Checkpoint Molecules') + 
  theme(axis.text.x = element_text(color = 'black', size =10),
        strip.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Baseline" = "Baseline\n(Not Inf)",
                              "Treated" = "Treated\n(Inf)")) +
  scale_fill_manual(values = c('Baseline' = '#637b31', 'Treated' = '#208CC0'), name = '') +
  facet_wrap(.~ response) +
  stat_compare_means(method = 'wilcox', paired = T, label='p.signif', label.x = 1.5, size = 7, vjust=1); p2

p3 <- score |> select(`M1-like`, tx_status, patient, response) |> 
  arrange(patient) |> 
  ggplot(aes(y=`M1-like`, x = tx_status)) +
  geom_boxplot(aes(fill = tx_status), width=0.5, outliers = F, show.legend = F) + 
  geom_line(aes(group = patient), color = "lightgray",linetype = "dashed", size= 0.4) +
  geom_point(size = 0.3) +
  theme_classic() + 
  ylab('') + xlab('') + ggtitle('M1-like TAM') + 
  theme(axis.text.x = element_text(color = 'black', size =10),
        strip.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Baseline" = "Baseline\n(Not Inf)",
                              "Treated" = "Treated\n(Inf)")) +
  scale_fill_manual(values = c('Baseline' = '#637b31', 'Treated' = '#208CC0'), name = '') +
  facet_wrap(.~ response) +
  stat_compare_means(method = 'wilcox', paired = T, label='p.signif', label.x = 1.5, size = 7, vjust=1); p3

p4 <- score |> select(`TLS formation`, tx_status, patient, response) |> 
  arrange(patient) |> 
  ggplot(aes(y=`TLS formation`, x = tx_status)) +
  geom_boxplot(aes(fill = tx_status), width=0.5, outliers = F, show.legend = F) + 
  geom_line(aes(group = patient), color = "lightgray",linetype = "dashed", size= 0.4) +
  geom_point(size = 0.3) +
  theme_classic() + 
  ylab('') + xlab('') + ggtitle('TLS formation') + 
  theme(axis.text.x = element_text(color = 'black', size =10),
        strip.text = element_text(size = 12),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Baseline" = "Baseline\n(Not Inf)",
                              "Treated" = "Treated\n(Inf)")) +
  scale_fill_manual(values = c('Baseline' = '#637b31', 'Treated' = '#208CC0'), name = '') +
  facet_wrap(.~ response) +
  stat_compare_means(method = 'wilcox', paired = T, label='p.signif', label.x = 1.5, size = 7, vjust=1); p4
p1+p2+p3+p4+plot_layout(nrow=1)
ggsave('figures/Dynamics/boxplot_inflamed_score.pdf', height = 4, width = 12)

# Dynamics
# subtype_order <- c('Immune Quiescent', 'Immune Inflamed', 'B Cell-enriched', 'Myeloid-enriched')
sample_info <- sample_info |>
  mutate(group = case_when(group == 'Immune Quiescent' ~ 'TIME-Q',
                           group == 'Immune Inflamed' ~ 'TIME-I',
                           group == 'B Cell-enriched' ~ 'TIME-B',
                           group == 'Myeloid-enriched' ~ 'TIME-Mye'))
subtype_order <- c('TIME-Q','TIME-I','TIME-B','TIME-Mye')
color.use <- structure(names = subtype_order, met.brewer("Juarez", length(unique(sample_info$TIME))))

pt <- sample_info |>
  filter(!is.na(tx_status)) |> 
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |> 
  pull(patient) 

ht_res_dyn <- function(mat_r, mat_nr, subtype_order, color.use){
  # Ensure consistent row and column order and presence based on subtype_order
  existing_rows_r <- subtype_order[subtype_order %in% rownames(mat_r)]
  mat_r <- mat_r[existing_rows_r, , drop = FALSE]
  existing_cols_r <- subtype_order[subtype_order %in% colnames(mat_r)]
  mat_r <- mat_r[, existing_cols_r, drop = FALSE]
  
  # Add any missing subtypes from subtype_order with all-zero rows/columns
  # For rows:
  missing_rows_r <- subtype_order[!subtype_order %in% rownames(mat_r)]
  if (length(missing_rows_r) > 0) {
    temp_mat_r_rows <- matrix(0, nrow = length(missing_rows_r), ncol = ncol(mat_r), 
                              dimnames = list(missing_rows_r, colnames(mat_r)))
    mat_r <- rbind(mat_r, temp_mat_r_rows)
    mat_r <- mat_r[subtype_order[subtype_order %in% rownames(mat_r)], , drop = FALSE] # Reorder
  }
  # For columns:
  missing_cols_r <- subtype_order[!subtype_order %in% colnames(mat_r)]
  if (length(missing_cols_r) > 0) {
    temp_mat_r_cols <- matrix(0, nrow = nrow(mat_r), ncol = length(missing_cols_r),
                              dimnames = list(rownames(mat_r), missing_cols_r))
    mat_r <- cbind(mat_r, temp_mat_r_cols)
    mat_r <- mat_r[, subtype_order[subtype_order %in% colnames(mat_r)], drop = FALSE] # Reorder
  }
  # Final reorder to ensure matrix conforms to subtype_order exactly
  mat_r <- mat_r[subtype_order, subtype_order]
  mat_r[is.na(mat_r)] <- 0 # Ensure NAs introduced by full join logic are 0
  
  
  
  # Ensure consistent row and column order and presence
  existing_rows_nr <- subtype_order[subtype_order %in% rownames(mat_nr)]
  mat_nr <- mat_nr[existing_rows_nr, , drop = FALSE]
  existing_cols_nr <- subtype_order[subtype_order %in% colnames(mat_nr)]
  mat_nr <- mat_nr[, existing_cols_nr, drop = FALSE]
  
  # Add any missing subtypes from subtype_order with all-zero rows/columns
  missing_rows_nr <- subtype_order[!subtype_order %in% rownames(mat_nr)]
  if (length(missing_rows_nr) > 0) {
    temp_mat_nr_rows <- matrix(0, nrow = length(missing_rows_nr), ncol = ncol(mat_nr),
                               dimnames = list(missing_rows_nr, colnames(mat_nr)))
    mat_nr <- rbind(mat_nr, temp_mat_nr_rows)
    mat_nr <- mat_nr[subtype_order[subtype_order %in% rownames(mat_nr)], , drop = FALSE]
  }
  missing_cols_nr <- subtype_order[!subtype_order %in% colnames(mat_nr)]
  if (length(missing_cols_nr) > 0) {
    temp_mat_nr_cols <- matrix(0, nrow = nrow(mat_nr), ncol = length(missing_cols_nr),
                               dimnames = list(rownames(mat_nr), missing_cols_nr))
    mat_nr <- cbind(mat_nr, temp_mat_nr_cols)
    mat_nr <- mat_nr[, subtype_order[subtype_order %in% colnames(mat_nr)], drop = FALSE]
  }
  mat_nr <- mat_nr[subtype_order, subtype_order]
  mat_nr[is.na(mat_nr)] <- 0
  
  #  4. Calculate Global Maxima for Barplot Scales 
  row_sums_r <- rowSums(abs(mat_r))
  col_sums_r <- colSums(abs(mat_r))
  
  row_sums_nr <- rowSums(abs(mat_nr))
  col_sums_nr <- colSums(abs(mat_nr))
  
  # Ensure all subtype names are present in sums vectors for max calculation, fill with 0 if not
  all_possible_sums_r_rows <- numeric(length(subtype_order))
  names(all_possible_sums_r_rows) <- subtype_order
  all_possible_sums_r_rows[names(row_sums_r)] <- row_sums_r
  
  all_possible_sums_nr_rows <- numeric(length(subtype_order))
  names(all_possible_sums_nr_rows) <- subtype_order
  all_possible_sums_nr_rows[names(row_sums_nr)] <- row_sums_nr
  
  all_possible_sums_r_cols <- numeric(length(subtype_order))
  names(all_possible_sums_r_cols) <- subtype_order
  all_possible_sums_r_cols[names(col_sums_r)] <- col_sums_r
  
  all_possible_sums_nr_cols <- numeric(length(subtype_order))
  names(all_possible_sums_nr_cols) <- subtype_order
  all_possible_sums_nr_cols[names(col_sums_nr)] <- col_sums_nr
  
  
  global_max_row_sum <- max(c(all_possible_sums_r_rows, all_possible_sums_nr_rows), na.rm = TRUE)
  global_max_col_sum <- max(c(all_possible_sums_r_cols, all_possible_sums_nr_cols), na.rm = TRUE)
  
  # Handle case where all sums are 0 to avoid ylim issue
  if (global_max_row_sum == 0) global_max_row_sum <- 1 
  if (global_max_col_sum == 0) global_max_col_sum <- 1
  
  # Define Annotations for ht1 (Responders) with Consistent Scales 
  df.col_r <- data.frame(group = colnames(mat_r)); rownames(df.col_r) <- colnames(mat_r)
  df.row_r <- data.frame(group = rownames(mat_r)); rownames(df.row_r) <- rownames(mat_r)
  
  col_annotation_r <- HeatmapAnnotation(df = df.col_r, 
                                        col = list(group = color.use),
                                        which = "column",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation_r <- HeatmapAnnotation(df = df.row_r, 
                                        col = list(group = color.use),
                                        which = "row",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))
  
  ha1_r <- rowAnnotation(Strength = anno_barplot(rowSums(abs(mat_r)), 
                                                 ylim = c(0, global_max_row_sum), 
                                                 border = FALSE,
                                                 gp = gpar(fill = color.use[rownames(mat_r)], 
                                                           col = color.use[rownames(mat_r)])), 
                         show_annotation_name = FALSE)
  
  ha2_r <- HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat_r)), 
                                                     ylim = c(0, global_max_col_sum),
                                                     border = FALSE,
                                                     gp = gpar(fill = color.use[colnames(mat_r)], 
                                                               col = color.use[colnames(mat_r)])), 
                             show_annotation_name = FALSE)
  
  # Define ht1 (Responder Heatmap) 
  ht1 <- Heatmap(t(scale(t(mat_r))), # Scaling row-wise as in your original example (t(scale(t(mat))))
                 row_title = 'Baseline', column_title = 'Responder', show_heatmap_legend = FALSE,
                 cluster_rows = FALSE, cluster_columns = FALSE, 
                 col = RColorBrewer::brewer.pal(8,'Reds'),  na_col = 'white',
                 row_names_side = 'left',
                 bottom_annotation = col_annotation_r,
                 left_annotation = row_annotation_r, # This is the colored bar for rows by subtype
                 top_annotation = ha2_r,    # This is the barplot for colSums
                 right_annotation = ha1_r,   # This is the barplot for rowSums
                 row_names_gp = gpar(fontsize = 10),
                 column_names_gp = gpar(fontsize = 10),
                 column_names_rot = 30, # Consider 45 or 90 if names overlap
                 width = ncol(mat_r)*unit(10, "mm"), 
                 height = nrow(mat_r)*unit(10, "mm"),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%d", mat_r[i, j]), x, y, gp = gpar(fontsize = 12))
                 })
  
  # Define Annotations for ht2 (Non-Responders) with Consistent Scales 
  df.col_nr <- data.frame(group = colnames(mat_nr)); rownames(df.col_nr) <- colnames(mat_nr)
  df.row_nr <- data.frame(group = rownames(mat_nr)); rownames(df.row_nr) <- rownames(mat_nr)
  
  col_annotation_nr <- HeatmapAnnotation(df = df.col_nr, 
                                         col = list(group = color.use),
                                         which = "column",
                                         show_legend = FALSE, show_annotation_name = FALSE,
                                         simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation_nr <- HeatmapAnnotation(df = df.row_nr, 
                                         col = list(group = color.use),
                                         which = "row",
                                         show_legend = FALSE, show_annotation_name = FALSE,
                                         simple_anno_size = grid::unit(0.2, "cm"))
  
  ha1_nr <- rowAnnotation(Strength = anno_barplot(rowSums(abs(mat_nr)), 
                                                  ylim = c(0, global_max_row_sum), 
                                                  border = FALSE,
                                                  gp = gpar(fill = color.use[rownames(mat_nr)], 
                                                            col = color.use[rownames(mat_nr)])), 
                          show_annotation_name = FALSE)
  
  ha2_nr <- HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat_nr)), 
                                                      ylim = c(0, global_max_col_sum),
                                                      border = FALSE,
                                                      gp = gpar(fill = color.use[colnames(mat_nr)], 
                                                                col = color.use[colnames(mat_nr)])), 
                              show_annotation_name = FALSE)
  
  # Define ht2 (Non-Responder Heatmap) 
  ht2 <- Heatmap(t(scale(t(mat_nr))), # Scaling row-wise
                 row_title = 'Baseline', column_title = 'Non-responder', show_heatmap_legend = FALSE,
                 cluster_rows = FALSE, cluster_columns = FALSE, 
                 col = RColorBrewer::brewer.pal(8,'Blues'), na_col = 'lightgray',
                 row_names_side = 'left',
                 bottom_annotation = col_annotation_nr,
                 left_annotation = row_annotation_nr,
                 top_annotation = ha2_nr,
                 right_annotation = ha1_nr,
                 row_names_gp = gpar(fontsize = 10),
                 column_names_gp = gpar(fontsize = 10),
                 column_names_rot = 30,
                 width = ncol(mat_nr)*unit(10, "mm"), 
                 height = nrow(mat_nr)*unit(10, "mm"),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%d", mat_nr[i, j]), x, y, gp = gpar(fontsize = 12))
                 })
  return(ht1+ht2)
}
r <- sample_info |>
  filter(paired == "Yes", response == 'R', res_metric != 'T-cell expansion') |>
  select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |>
  group_by(Baseline, Treated) |>
  dplyr::summarise(count = n(), .groups = "drop") |>
  pivot_wider(names_from = 'Treated', values_from = count, values_fill = 0) |>
  column_to_rownames("Baseline") 
nr <- sample_info |>
  filter(paired == "Yes", response == 'NR') |>
  # filter(paired == "Yes", response == 'NR') |>
  select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |>
  group_by(Baseline, Treated) |>
  dplyr::summarise(count = n(), .groups = "drop") |>
  pivot_wider(names_from = 'Treated', values_from = count, values_fill = 0) |>
  column_to_rownames("Baseline") |> as.matrix() 
pdf('figures/Fig4/sc.pdf', height = 4.5)
ht_res_dyn(mat_r = r, mat_nr = nr, subtype_order = subtype_order, 
           color.use = structure(names = subtype_order, met.brewer("Juarez", 4)))
grid.text("Pan-cancer(scRNA-seq)", x = 0.5, y = 0.95, gp = gpar(fontsize = 16))
dev.off()

sample_info_bulk <- read.csv('tables/sample_info_time_bulk.csv') 
subtype_order <- c('Immune Quiescent', 'Immune Inflamed', 'B Cell-enriched', 'Myeloid-enriched')
rownames(sample_info_bulk) <- NULL
sample_info_bulk$patient <- paste0(sample_info_bulk$cohort, '_', sample_info_bulk$patient)

pt <- sample_info_bulk |>
  filter(!is.na(tx_status)) |> 
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |> 
  pull(patient) 

# Response (scRNA)
sample_info <- read.csv('tables/sample_info_time.csv')
sample_info_val <- read.csv('tables/sample_info_validation.csv')
sample_info_val <- filter(sample_info_val, response != 'NE')
sample_info <- rbind(sample_info_val[,c('sample','group','response','tx_status', "patient", 'paired')], sample_info[,c('sample','group','response','tx_status', "patient", 'paired')])

pt_dyn_r <- sample_info |>
  filter(paired == "Yes", !patient %in% pt_duplicated$patient) |> 
  select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Baseline == 'TIME-BE' & Treated == 'TIME-BE'|
           Baseline == 'TIME-BE' & Treated == 'TIME-Teff'|
           Baseline == 'TIME-Teff' & Treated == 'TIME-BE'|
           Baseline == 'TIME-Mye' & Treated == 'TIME-BE'|
           Baseline == 'TIME-IE' & Treated == 'TIME-Teff'|
           Baseline == 'TIME-ISG' & Treated == 'TIME-Teff'
  ) |>
  select(patient) |> pull()

pt_dyn_nr <- sample_info |>
  filter(paired == "Yes", !patient %in% pt_duplicated$patient) |>
  select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Baseline == 'TIME-BE' & Treated == 'TIME-IE'|
           Baseline == 'TIME-BE' & Treated == 'TIME-Mye'|
           Baseline == 'TIME-Mye' & Treated == 'TIME-Mye') |> 
  select(patient) |> pull()
statistic <- sample_info |> 
  filter(response %in% c('R','NR'),
         paired == "Yes", !patient %in% pt_duplicated$patient) |> 
  mutate(response = factor(response, levels = c('R','NR')),
         dynamics = case_when(patient %in% pt_dyn_r ~ 'Fav-Resp',
                              patient %in% pt_dyn_nr ~ 'Unfav-Resp',
                              !patient %in% c(pt_dyn_r, pt_dyn_nr) ~ 'Others')) |> 
  mutate(dynamics = factor(dynamics, levels = c('Fav-Resp','Unfav-Resp','Others'))) |> 
  filter(dynamics %in% c('Fav-Resp','Unfav-Resp')) |> 
  distinct(patient, .keep_all = T) |> 
  tabyl(dynamics, response) |> data.frame() |> column_to_rownames(var = 'dynamics') |> filter(R !=0) |> chisq.test()
pvalue <- statistic$p.value
sample_info |> 
  filter(response %in% c('R','NR'),
         paired == "Yes", !patient %in% pt_duplicated$patient) |> 
  mutate(response = factor(response, levels = c('R','NR')),
         dynamics = case_when(patient %in% pt_dyn_r ~ 'Fav-Resp',
                              patient %in% pt_dyn_nr ~ 'Unfav-Resp',
                              !patient %in% c(pt_dyn_r, pt_dyn_nr) ~ 'Others')) |> 
  mutate(dynamics = factor(dynamics, levels = c('Fav-Resp','Unfav-Resp','Others'))) |> 
  filter(dynamics %in% c('Fav-Resp','Unfav-Resp')) |>
  distinct(patient, .keep_all = T) |> 
  select(patient, dynamics, response) |> 
  group_by(dynamics, response) |>
  tally(name = "n") |> 
  group_by(dynamics) |>
  mutate(prop = n / sum(n)) |>
  ggplot(aes(x = dynamics, y = prop, fill = response)) +
  geom_col(position = "fill", color='black', width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(
    aes(label = n), 
    position = position_stack(vjust = 0.5), 
    size = 3.5, 
    color = "black"
  ) +
  labs(
    x = "",
    y = "",
    fill = "Subtype",
    title = ""
  ) +
  theme_classic() + scale_fill_manual(values = c('R' = '#CC0C00FF','NR' = '#5C88DAFF'), name = 'Response') + ggtitle('') +
  theme(axis.text.y = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 0, hjust = 0.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) + geom_signif(comparison = list(c('Fav-Resp','Unfav-Resp')), y_position = 1.01, annotations=paste0('p<0.001'))
ggsave('figures/Dynamics/barplot_group_sc.pdf', height = 3.5, width = 4)

# Accordance 
sample_info_bulk <- read.csv('tables/sample_info_time_bulk.csv')
pt <- sample_info_bulk |>
  filter(!is.na(tx_status)) |> 
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |> 
  pull(patient) 
pt_duplicated <- sample_info_bulk |>
  filter(paired == "Yes"&response == "NR" | paired == "Yes"&response == "R") |>
  dplyr::count(patient, tx_status) |>
  filter(n > 1)
cutoff <- 2
df_bulk <- sample_info_bulk |> 
  dplyr::filter(paired == 'Yes', !is.na(response), response != 'NE', !patient %in% pt_duplicated) |> 
  dplyr::select(tx_status, response, patient, group) |> 
  pivot_wider(names_from = 'tx_status', values_from = 'group') |> 
  mutate(dynamics = paste0(Baseline, '->', Treated)) |> 
  group_by(dynamics) |>
  summarise(
    R_count = sum(response == 'R'),
    NR_count = sum(response == 'NR'),
    response_rate = R_count / (R_count + NR_count)
  ) |> 
  filter(R_count>=cutoff | NR_count>=cutoff) |> 
  mutate(R_count_norm = R_count / sum(R_count),
         NR_count_norm = NR_count / sum(NR_count))
df_sc <- sample_info |> 
  filter(paired == 'Yes', 
         cohort != 'HCC_Guo',
         # res_metric != 'T-cell expansion',
         response != 'NE') |> 
  dplyr::select(tx_status, response, patient, group) |> 
  pivot_wider(names_from = 'tx_status', values_from = 'group') |> 
  mutate(dynamics = paste0(Baseline, '->', Treated)) |> 
  group_by(dynamics) |>
  summarise(
    R_count = sum(response == 'R'),
    NR_count = sum(response == 'NR'),
    response_rate = R_count / (R_count + NR_count)
  ) |> 
  filter(R_count>=cutoff | NR_count>=cutoff) |> 
  mutate(R_count_norm = R_count / sum(R_count),
         NR_count_norm = NR_count / sum(NR_count))
df <- merge(df_sc, df_bulk, by = 'dynamics')
names(df) <- str_replace_all(names(df), '\\.x', '_sc')
names(df) <- str_replace_all(names(df), '\\.y', '_bulk')
df |> 
  mutate(dynamics = str_replace_all(dynamics, 'B Cell-enriched', 'B'),
         dynamics = str_replace_all(dynamics, 'Immune Inflamed', 'I'),
         dynamics = str_replace_all(dynamics, 'Immune Quiescent', 'Q'),
         dynamics = str_replace_all(dynamics, 'Myeloid-enriched', 'Mye')) |> 
  mutate(resp_group = case_when(dynamics %in% c('B->B', 'I->B','I->Q', 'Q->I', 'Mye->I') ~ 'Resp-fav',
                                dynamics %in% c('Q->Q','Mye->Mye','Q->Mye','B->Q','I->I') ~ 'Resp-unfav',
                                .default = 'Unresolved')) |> 
  filter(resp_group != 'Unresolved') |>
  ggplot(aes(x = response_rate_sc, y=response_rate_bulk)) + 
  geom_smooth(method = "lm", color = "black", se = F) +
  # geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  scale_fill_manual(values = c('Resp-fav' = '#CC0C00FF','Resp-unfav' = '#5C88DAFF'), name = 'Response') +
  theme_bw() + 
  stat_cor(size = 6) + 
  xlim(0,1) + ylim(0,0.8) + 
  ggtitle('') +
  xlab('scRNA-seq') + ylab('Bulk Transcriptomic') + 
  coord_cartesian(clip = "off") +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 15), "pt")) +
  geom_point(aes(fill = resp_group), size = 4, pch =21) +
  ggrepel::geom_text_repel(aes(label = dynamics), size = 5, box.padding = 0.6) +
  theme(legend.position = c(0.85, 0.2))
ggsave('figures/Dynamics/group_res/cor_sc_bulk.pdf', height = 5, width = 5)

# Response (bulk)
sample_info <- read.csv('tables/sample_info_time_bulk.csv') 
rownames(sample_info) <- NULL
sample_info$patient <- paste0(sample_info$cohort, '_', sample_info$patient)
pt <- sample_info |>
  filter(!is.na(tx_status)) |> 
  group_by(patient) |>
  dplyr::summarise(n = n()) |>
  filter(n==2) |> 
  pull(patient) 
sample_info$paired <- ifelse(sample_info$patient %in% pt, 'Yes', 'No')
pt_duplicated <- sample_info |>
  filter(paired == "Yes"&response == "NR" | paired == "Yes"&response == "R") |>
  dplyr::count(patient, tx_status) |>
  filter(n > 1) 
c('Immune Quiescent', 'Immune Inflamed', 'B Cell-enriched', 'Myeloid-enriched')
pt_dyn_r <- sample_info |>
  filter(paired == "Yes", !patient %in% pt_duplicated$patient) |> 
  select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Baseline == 'Immune Inflamed' & Treated == 'B Cell-enriched' |
           Baseline == 'Myeloid-enriched' & Treated %in% c('Immune Inflamed','B Cell-enriched') |
           Baseline == 'Immune Quiescent' & Treated %in% c('Immune Inflamed','B Cell-enriched') |
           Baseline == 'B Cell-enriched' & Treated == 'B Cell-enriched' |
           Baseline == 'B Cell-enriched' & Treated == 'Immune Inflamed') |>
  select(patient) |> pull()

pt_dyn_nr <- sample_info |>
  filter(paired == "Yes", !patient %in% pt_duplicated$patient) |>
  select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Treated == 'Myeloid-enriched'|
           Baseline == 'Immune Quiescent' & Treated=='Immune Quiescent'|
           Baseline == 'B Cell-enriched' & Treated=='Immune Quiescent'
  ) |>
  select(patient) |> pull()

statistic <- sample_info |> 
  filter(response %in% c('R','NR'),
         paired == "Yes", !patient %in% pt_duplicated$patient) |> 
  mutate(response = factor(response, levels = c('R','NR')),
         dynamics = case_when(patient %in% pt_dyn_r ~ 'Fav-Resp',
                              patient %in% pt_dyn_nr ~ 'Unfav-Resp',
                              !patient %in% c(pt_dyn_r, pt_dyn_nr) ~ 'Others')) |> 
  mutate(dynamics = factor(dynamics, levels = c('Fav-Resp','Unfav-Resp','Others'))) |> 
  filter(dynamics %in% c('Fav-Resp','Unfav-Resp')) |> 
  distinct(patient, .keep_all = T) |> 
  tabyl(dynamics, response) |> data.frame() |> column_to_rownames(var = 'dynamics') |> filter(R !=0) |> chisq.test()
pvalue <- statistic$p.value
sample_info |> 
  filter(response %in% c('R','NR'),
         paired == "Yes", !patient %in% pt_duplicated$patient) |> 
  mutate(response = factor(response, levels = c('R','NR')),
         dynamics = case_when(patient %in% pt_dyn_r ~ 'Fav-Resp',
                              patient %in% pt_dyn_nr ~ 'Unfav-Resp',
                              !patient %in% c(pt_dyn_r, pt_dyn_nr) ~ 'Others')) |> 
  mutate(dynamics = factor(dynamics, levels = c('Fav-Resp','Unfav-Resp','Others'))) |> 
  filter(dynamics %in% c('Fav-Resp','Unfav-Resp')) |>
  distinct(patient, .keep_all = T) |> 
  select(patient, dynamics, response) |> 
  group_by(dynamics, response) |>
  tally(name = "n") |> 
  group_by(dynamics) |>
  mutate(prop = n / sum(n)) |>
  ggplot(aes(x = dynamics, y = prop, fill = response)) +
  geom_col(position = "fill", color='black', width = 0.6) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(
    aes(label = n), 
    position = position_stack(vjust = 0.5), 
    size = 3.5, 
    color = "black"
  ) +
  labs(
    x = "",
    y = "",
    fill = "Subtype",
    title = ""
  ) +
  theme_classic() + scale_fill_manual(values = c('R' = '#CC0C00FF','NR' = '#5C88DAFF'), name = 'Response') + ggtitle('Bulk Cohorts') +
  theme(axis.text.y = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 0, hjust = 0.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5)) + geom_signif(comparison = list(c('Fav-Resp','Unfav-Resp')), y_position = 1.01, annotations=paste0('p=', round(pvalue,4)))
ggsave('figures/Dynamics/barplot_group_bulk.pdf', height = 3.5, width = 4)

# Mye2Mye
pt_dyn_m2m <- sample_info |>
  filter(paired == "Yes", !patient %in% pt_duplicated$patient, cancertype=='Melanoma') |> 
  select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Treated == 'Myeloid-enriched') |>
  select(patient) |> pull()
pt_dyn_m2o <- sample_info |>
  filter(paired == "Yes", !patient %in% pt_duplicated$patient) |>
  select(patient, tx_status, group) |>
  pivot_wider(names_from = tx_status, values_from = group) |> 
  filter(Baseline == 'Myeloid-enriched' & Treated != 'Myeloid-enriched') |> 
  select(patient) |> pull()

sample_info_os <- sample_info |> 
  filter(!is.na(os), !is.na(time_os),
         paired == "Yes",
         !patient %in% pt_duplicated$patient) |> 
  mutate(os = case_when(os == 'Alive' ~ 0,
                        os == 'Dead' ~ 1,
                        TRUE ~ as.numeric(os)),
         dynamics = case_when(patient %in% pt_dyn_m2m ~ 'Mye-stable',
                              !patient %in% pt_dyn_m2m ~ 'Others')) |>
  # mutate(dynamics = factor(dynamics, levels = c('Others','Mye-shifted','Mye-stable'))) |> 
  distinct(patient, .keep_all = T) |> filter(dynamics %in% c('Others','Mye-stable'))
fit <- survfit(Surv(time_os, os) ~ dynamics, data = sample_info_os)
p <- ggsurvplot(fit, data = sample_info_os,conf.int = F, pval = T, 
                # legend.labs = c('Others','Mye-shifted','Mye-stable'),
                risk.table = TRUE,        
                risk.table.col = "strata",
                risk.table.height = 0.25,
                ggtheme = theme_bw(),
                # palette = c('#CC0C00FF','#007900','#5C88DAFF'),
                xlab = 'Time(OS)', title = 'Melanoma');p
pdf('figures/Dynamics/km_group_bulk_mye.pdf', height = 5, width = 5)
print(p, newpage = FALSE)
dev.off()










