rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','tibble','qs2','janitor','RColorBrewer')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(max.print = 10000)
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 'SCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'NSCLC_Liu',
              'PCa_Hawley')
metadata_list <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_r2.qs2')) 
  metadata <- seu@meta.data |> 
    select(celltype_r2, celltype_main, patient, time_point, sample, cohort) |> 
    filter(celltype_r2 != 'NK_CD56hiCD16hi')
  return(metadata)
})

metadata <- do.call(rbind, metadata_list)

label_malignant <- c('Malignant(CNA+)','Epithelial(CNA-)','Melanocytes(CNA-)')
metadata <- metadata |> 
  group_by(sample) |> 
  mutate(cell_count=n()) |> 
  mutate(malignant_count = sum(celltype_r2 %in% label_malignant),
         non_malignant_count = sum(!celltype_r2 %in% label_malignant)) |> 
  mutate(malignant_pct = malignant_count/cell_count,
         non_malignant_pct = non_malignant_count/cell_count,
         keep = case_when(non_malignant_count >= 200 ~ 'Yes',
                          non_malignant_count < 200 ~ 'No')) |> 
  ungroup() |> 
  filter(keep == 'Yes')

# filtering pt with paired samples 
pt <- metadata |>
  distinct(sample, .keep_all = TRUE) |>
  group_by(patient) |>
  summarise(n = n()) |>
  filter(n==2) |>
  pull(patient) |> length()




