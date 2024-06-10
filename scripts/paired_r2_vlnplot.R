rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggplot2','rstatix','ggpubr')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

meta_int <- read.csv('/bigdata/zlin/Melanoma_meta/tables/meta_int.csv') 
# modality
df_add <- meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 == 'Mono_CD16') |> 
  select(freq_r2_comp, time_point, patient, celltype_r2, modality) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(direction = ifelse(Post > Pre*1.1, 'Up', ifelse(Post < Pre*0.9, 'Down', 'not'))) |>
  group_by(celltype_r2, modality) |> 
  mutate(pt_count = n()) |> 
  group_by(direction, .add = T) |> 
  mutate(dir_count = n()) |> 
  distinct(celltype_r2, direction, pt_count, dir_count, modality) |> 
  filter(!direction == 'not') |> 
  pivot_wider(values_from = 'dir_count', names_from = 'direction')

stat.test <- meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% c('Mono_CD16')) |> 
  select(freq_r2_comp, time_point, modality,patient) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  pivot_longer(cols = c('Pre', 'Post'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  group_by(modality) |> 
  t_test(freq_r2_comp~time_point, paired = T) |> 
  add_significance("p") |> 
  merge(df_add, by = 'modality') |> 
  mutate(strip_label = paste0('Up: ', Up, '/', pt_count, '\n Down: ', Down, '/', pt_count, '\n p=', p)) 

meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% c('Mono_CD16')) |> 
  select(freq_r2_comp, time_point, modality,patient) |> 
  merge(stat.test, by = 'modality') |> 
  ggplot(aes(x = factor(time_point, levels = c('Pre', 'Post')), y = freq_r2_comp)) +
  geom_violin(aes(fill = time_point), alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0.2) +
  geom_line(aes(group = patient), color = "gray",linetype = "dashed") +
  geom_point(size = 0.5) +
  scale_fill_manual(values = c('#154999', '#CF0034')) +
  facet_wrap(.~ modality + strip_label, scales = 'free') +
  xlab("") + ylab("Relative Frequency") + ggtitle('Mono_CD16') +
  theme_classic2() + 
  theme(axis.text.x = element_text(size = 8, colour = 'black'),
        strip.text = element_text(size = 8, margin = margin(2,1,2,1))) +
  guides(fill="none") 
ggsave('/bigdata/zlin/Melanoma_meta/figures/Abundance/boxplot_MonoCD16.pdf', height = 3, width = 4.5)

# response
df_add <- meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 == 'CD4_Treg_TNFRSF9') |> 
  select(freq_r2_comp, time_point, patient, celltype_r2, response) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |> 
  mutate(direction = ifelse(Post > Pre*1.1, 'Up', ifelse(Post < Pre*0.9, 'Down', 'not'))) |>
  group_by(celltype_r2, response) |> 
  mutate(pt_count = n()) |> 
  group_by(direction, .add = T) |> 
  mutate(dir_count = n()) |> 
  distinct(celltype_r2, direction, pt_count, dir_count, response) |> 
  filter(!direction == 'not') |> 
  pivot_wider(values_from = 'dir_count', names_from = 'direction')

stat.test <- meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% 'CD4_Treg_TNFRSF9') |> 
  select(freq_r2_comp, time_point, response, patient, celltype_r2) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = time_point, values_fill = 0) |>
  pivot_longer(cols = c('Pre', 'Post'), names_to = 'time_point', values_to = 'freq_r2_comp') |> 
  group_by(response) |> 
  t_test(freq_r2_comp~time_point, paired = T) |> 
  add_significance("p") |> 
  merge(df_add, by = c('response')) |> 
  mutate(strip_label = paste0('Up: ', Up, '/', pt_count, '\n Down: ', Down, '/', pt_count, '\n p=', p)) 

meta_int |> 
  distinct(sample, celltype_r2, .keep_all = T) |> 
  filter(celltype_r2 %in% c('CD4_Treg_TNFRSF9')) |> 
  select(freq_r2_comp, time_point, response,patient) |> 
  merge(stat.test, by = 'response') |> 
  ggplot(aes(x = factor(time_point, levels = c('Pre', 'Post')), y = freq_r2_comp)) +
  geom_violin(aes(fill = time_point), alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0.2) +
  geom_line(aes(group = patient), color = "gray",linetype = "dashed") +
  geom_point(size = 0.5) +
  scale_fill_manual(values = c('#154999', '#CF0034')) +
  facet_wrap(.~ response + strip_label, scales = 'free') + ggtitle('CD4_Treg_TNFRSF9') +
  xlab("") + ylab("Relative Frequency") + 
  theme_classic2() + 
  theme(axis.text.x = element_text(size = 8, colour = 'black'),
        strip.text = element_text(size = 8, margin = margin(2,1,2,1))) +
  guides(fill="none") 
ggsave('/bigdata/zlin/Melanoma_meta/figures/Abundance/boxplot_treg_tnfrsf9.pdf', height = 3, width = 4.5)








