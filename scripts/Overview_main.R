rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','dittoSeq','RColorBrewer','tibble','ggnewscale','qs','MetBrewer','forcats','lmerTest','grid')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

meta_patient <- read.csv('/bigdata/zlin/PanCancer_ICI/tables/meta_patient.csv')
df <- meta_patient |> 
  dplyr::mutate(cancertype = factor(cancertype, levels = c("SKCM", "BCC", "SCC", "TNBC", "ER+BC", "HER2+BC", "HNSC", "NSCLC", "PCa", "CRC")),
         response = factor(response, levels = c('RE', 'NR'))) |> 
  arrange(cancertype, interval) 

df <- df[nrow(df):1,]
df$y <- 1:nrow(df) 
n <- 5
dfs <- lapply(1:n, function(i) {
  df_rep <- df
  df_rep$x <- i
  return(df_rep)
})
final_df <- do.call(rbind, dfs)
final_df$x_cancertype <- NA
final_df$x_cancertype[final_df$x == 1] <- 0.5
final_df$x_treatment <- NA
final_df$x_treatment[final_df$x == 2] <- 1
final_df$x_prior <- NA
final_df$x_prior[final_df$x == 3] <- 1.5
final_df$x_response <- NA
final_df$x_response[final_df$x == 4] <- 2
final_df$x_res_metric <- NA
final_df$x_res_metric[final_df$x == 5] <- 2.5

p1 <- ggplot(final_df[nrow(final_df):1,]) +
  geom_tile(aes(x_cancertype,y, width=0.4, height=0.4, fill = cancertype)) +
  scale_fill_manual(values = c("PCa" = "#a82203", 
                               "CRC" = "#208cc0",     
                               "NSCLC" = "#f1af3a",  
                               "HNSC" = "#cf5e4e",
                               "SCC" = "#3B7546",     
                               "BCC" = "#0092A5", 
                               "SKCM" = "#000000", 
                               "TNBC" = "#FF0196", 
                               "ER+BC" = "#F073EA",    
                               "HER2+BC" = "#FFB0FF"), name = "Cancer Type") + 
  new_scale_fill() +  # This introduces a new fill scale
  geom_tile(aes(x_treatment,y, width=0.4, height=0.5, fill = treatment)) +
  scale_fill_npg(name="Treatment") +
  new_scale_fill() +  # This introduces a new fill scale
  geom_tile(aes(x_prior,y, width=0.4, height=0.5, fill = prior)) +
  scale_fill_manual(values = c('No' = "#9cc184",
                               'Yes' = "#3c7c3d"), name = "Prior Tx") +
  new_scale_fill() +  # This introduces a new fill scale
  geom_tile(aes(x_response,y, width=0.4, height=0.5, fill = response)) +
  scale_fill_startrek(name = "Response") +
  new_scale_fill() +  # This introduces a new fill scale
  geom_tile(aes(x_res_metric,y, width=0.4, height=0.5, fill = res_metric)) +
  scale_fill_manual(values = met.brewer("Juarez", length(unique(df$res_metric))), name = 'Response Metrix') +
  coord_cartesian(ylim = c(1, length(unique(final_df$patient)))) +
  scale_y_reverse() +
  annotate("text", x = c(0.5, 1, 1.5, 2, 2.5), y = 5, 
           label = c("Cancer Type", "Treatment", "Prior-Tx", "Response", "Response Metric"),
           vjust = 0.5, hjust = 1, angle = 45, size = 3) +
  theme_void() +
  theme(legend.key.size = unit(0.5, 'cm'),  # Adjust legend key size
        legend.text = element_text(size = 8))

df <- df[nrow(df):1,]
df$dataset[df$dataset %in% c('BCC_Yost','SCC_Yost')] <- 'BCC/SCC_Yost'
# df_p2 <- df |>
#   distinct(patient, .keep_all = T)
# df_p2$sequence <- rev(1:nrow(df_p2))
df$sequence <- rev(1:nrow(df))
p2 <- df |>
  ggplot(aes(x = sequence, y = pmin(interval, 42))) +  # Limit y to 40
  geom_segment(aes(x = sequence, xend = sequence, y = 0, yend = pmin(interval, 42)), color = 'black', linewidth = 0.2) +
  scale_color_manual(values = dittoColors()[1:length(unique(df$dataset))], name = 'Cohort') +
  scale_y_continuous(limits = c(0, 42), breaks = c(0, 7, 14, 21, 28, 35, 42), labels = c("0", "1", "2", "3", "4", "5", ">6")) +
  theme_minimal() + 
  geom_point(aes(color = dataset), size = 0.4) +
  geom_point(aes(x = sequence, y = 0, color = dataset), size = 0.4) +  # Add points at y = 0
  labs(x = '', y = 'Time from Treatment (week)') +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(vjust =10, size = 12),
        axis.text.x = element_text(vjust =12),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.2, 'cm')) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  geom_hline(yintercept = 21, linetype = 'dashed') + 
  coord_flip() 
p1+p2+plot_layout(widths = c(0.2, 0.6, 1.2), ncol = 2) 
ggsave('/bigdata/zlin/PanCancer_ICI/figures/dist_1.pdf', height = 12, width = 8)

meta_int <- read.csv('/bigdata/zlin/PanCancer_ICI/tables/meta_int.csv') 
meta_int$celltype_main[str_detect(meta_int$celltype_r2, 'CD4')] <- 'CD4+T'
meta_int$celltype_main[str_detect(meta_int$celltype_r2, 'CD8')] <- 'CD8+T'
meta_int$celltype_main[meta_int$celltype_main == 'T'] <- 'CD8+T'
meta_int$celltype_main[meta_int$celltype_main %in% c('Mono', 'Macro')] <- 'Mo/Mac'
metadata <- meta_int |> group_by(sample) |>
  select(sample, time_point, sample, patient, celltype_main, cancertype, response, interval) |> 
  dplyr::mutate(count_sample = n()) |> 
  group_by(celltype_main, .add = TRUE) |> 
  dplyr::mutate(count_main = n()) |> 
  dplyr::mutate(freq_main = count_main/count_sample) |> 
  distinct(sample, celltype_main, .keep_all = T)
metadata <- meta_int[nrow(meta_int):1, ]
metadata$time_point <- factor(metadata$time_point, levels = c('Pre','On'))

colors <- c('CD4+T' = '#E31A1C',
            'CD8+T' = '#1F78B4',
            'NK' = '#CAB2D6',
            'B' = '#A6CEE3',
            'Plasma' = '#FF7F00',
            'pDC' = '#B15928',
            'Mast' = '#FFFF99',
            'cDC' = '#B2DF8A',
            'Mo/Mac' = '#33A02C',
            'Endo' = '#FB9A99',
            'CAF' = '#FDBF6F')
celltype_order <- c('CD4+T', 'CD8+T', 'NK', 'B', 'Plasma', 'pDC', 'Mast', 'cDC', 'Mo/Mac', 'Endo', 'CAF')
metadata |> 
  dplyr::mutate(cancertype = factor(cancertype, levels = c("SKCM", "BCC", "SCC", "TNBC", "ER+BC", "HER2+BC", "HNSC", "NSCLC", "PCa", "CRC")),
                response = factor(response, levels = c('RE', 'NR')),
                patient = factor(patient, levels = rev(df$patient))) |>
  arrange(patient) |>
  ggplot(aes(x = fct_inorder(patient), y = freq_main, fill = factor(celltype_main, levels = celltype_order))) +
  geom_col(position = "fill", width = 0.4) + coord_flip() +
  scale_fill_manual(values = colors, name = 'Cell type') + facet_wrap(.~time_point) + theme_void() +
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        legend.key.size = unit(0.3, 'cm'))

ggsave('/bigdata/zlin/PanCancer_ICI/figures/dist_2.png', height = 12, width = 6)

# Change at main level
meta_int <- read.csv('/bigdata/zlin/PanCancer_ICI/tables/meta_int.csv')
meta_int$celltype_main[str_detect(meta_int$celltype_r2, 'CD4')] <- 'CD4+T'
meta_int$celltype_main[str_detect(meta_int$celltype_r2, 'CD8')] <- 'CD8+T'
meta_int$celltype_main[meta_int$celltype_main == 'T'] <- 'CD8+T'
# meta_int$celltype_main[meta_int$celltype_main %in% c('Mono', 'Macro')] <- 'Mo/Mac'
meta_int <- meta_int |> group_by(sample) |>
  select(sample, time_point, sample, patient, celltype_main, cancertype, response, interval, dataset, modality) |> 
  dplyr::mutate(count_sample = n()) |> 
  group_by(celltype_main, .add = TRUE) |> 
  dplyr::mutate(count_main = n()) |> 
  dplyr::mutate(freq_main = count_main/count_sample)
meta_int$int_cat <- ifelse(meta_int$interval < 21, '< 21d', '>= 21d')
uni_lmer <- function(meta_int){
  freq_mat <- meta_int |> 
    distinct(celltype_main, sample, .keep_all = T) |> 
    select(freq_main, dataset, response, modality, int_cat, patient, time_point, celltype_main) |> 
    pivot_wider(values_from = freq_main, names_from = time_point, values_fill = 0) |> 
    pivot_longer(cols = c('Pre', 'On'), names_to = 'time_point', values_to = 'freq_main') |> 
    group_by(celltype_main) |> 
    dplyr::mutate(freq_main_scale = scale(freq_main)) |> 
    ungroup()
  subtypes <- unique(freq_mat$celltype_main)
  uni_models <- lapply(subtypes, function(subtype){
    print(subtype)
    freq_mat$timepoint <- ifelse(freq_mat$time_point == 'Pre', 0, 1)
    formula <- as.formula(paste("timepoint ~ freq_main_scale + dataset + modality + int_cat + (1 | patient)"))
    model <- lmer(formula, data = freq_mat[freq_mat$celltype_main == subtype,], REML = FALSE)
    model_summary <- summary(model)
    print('summary')
    coef_table <- coef(summary(model))
    confint_table <- confint(model, level = 0.95)
    subtype_results <- coef_table['freq_main_scale', , drop = FALSE]
    subtype_confint <- confint_table['freq_main_scale', , drop = FALSE]
    combined_results <- data.frame(
      Celltypes = subtype,
      Estimate = subtype_results[1],
      StdError = subtype_results[2],
      tValue = subtype_results[4],
      pValue = coef_table['freq_main_scale', 'Pr(>|t|)'],
      CI_lower = subtype_confint[1],
      CI_upper = subtype_confint[2]
    )
    return(combined_results)
  })
  df <- do.call(rbind, uni_models) |> data.frame()
  df$fdr <- p.adjust(df$pValue, method = 'fdr', n = nrow(df)) 
  # celltype_keep <- filter(df, pValue < 0.05) |> pull(Celltypes) |> unique()
  # df <- filter(df, Celltypes %in% celltype_keep)
  return(df)
}
# Overall
res_lmer <- uni_lmer(meta_int); res_lmer
res_lmer <- res_lmer |> arrange(dplyr::desc(Estimate))
# res_lmer$CI_lower[res_lmer$CI_lower < -1] <- -1
# res_lmer$CI_upper[res_lmer$CI_upper > 0.5] <- 0.5
pdf(paste0('/bigdata/zlin/PanCancer_ICI/figures/Change/uni_lmer_main.pdf'), height = 4, width = 4)
p <- ggplot(res_lmer, aes(x= factor(Celltypes, levels = rev(res_lmer$Celltypes)), y=Estimate, ymin=CI_lower, ymax=CI_upper)) +
  # ggplot(res_lmer, aes(x= factor(Celltypes, levels = rev(res_lmer$Celltypes)), y=Estimate, ymin=CI_lower, ymax=CI_upper, size = -log10(pValue))) +
  geom_linerange(size=1, position=position_dodge(width = 0.5)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(shape=21, fill="#364F77", color='white', stroke = 0.5, position=position_dodge(width = 0.5), size = 4) +
  # geom_point(shape=21, fill="#0047AB", color='white', stroke = 0.5, position=position_dodge(width = 0.5)) +
  # scale_size(breaks = c(-log10(0.1), -log10(0.5)), labels = c('<0.1', '<0.5'), range = c(2,5)) +
  scale_y_continuous(name= "Effect Size", limits = c(-0.09, 0.08)) +
  scale_color_manual(values = c("black", "grey")) +
  coord_flip() + theme_minimal() +
  labs(size = "P value") + xlab("") +
  # guides(size = guide_legend(override.aes = list(fill = "#0047AB"))) +
  theme(plot.margin = unit(c(0.5, 1, 2, 0.5), "lines"),
        axis.title.x = element_text(margin = margin(t = 2, unit = "lines")),
        axis.text.y = element_text(size = 9, colour = "black"))
grid.draw(p)
y_text = 0.18
x1 = 0.25
x2 = 0.9
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
