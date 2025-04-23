library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(janitor)
library(survminer)
library(survival)
library(ggplot2)
dirs <- c(paste0('data/', 
                 c('bulk_datasets/Melanoma/GSE91061/', 
                   'bulk_datasets/TIGER/Melanoma-Nathanson_2017/',
                   'bulk_datasets/TIGER/Melanoma-phs000452/',
                   'bulk_datasets/TIGER/Melanoma-GSE78220/',
                   'bulk_datasets/TIGER/Melanoma-GSE115821/',
                   'bulk_datasets/Melanoma/Liu/',
                   'bulk_datasets/Melanoma/Ribas/',
                   'bulk_datasets/TIGER/Melanoma-PRJEB23709/',
                   'carroll_etal_2023/processed_data/',
                   'bulk_datasets/IMvigor210/',
                   'bulk_datasets/IMmotion150/',
                   'bulk_datasets/RCC_Braun/',
                   'bulk_datasets/NSCLC_Kang/'
                   ), 'comb.csv'))
genelist <- read.csv('tables/cosg_genes.csv', check.names = F) |> as.list()
genelist$`Melanocytes(CNA+)` <- NULL
genelist$`Epithelial(CNA+)` <- NULL
genelist$Cycling <- NULL
colorder <- c('sample', 'patient', 'response', 'os', 'time_os', 'tx', 'tx_status', 'cancertype', 'cohort', names(genelist))
comb_list <- lapply(dirs, function(x){
  print(x)
  comb <- read.csv(x, check.names = F)
  comb <- comb[,colorder]
  return(comb)
})
comb <- do.call(rbind, comb_list)

df <- comb |> 
  drop_na(os) |> 
  filter(tx_status == 'Baseline') |>
  mutate(os = ifelse(os %in% c('Alive', 0), 0, 1)) 
colnames(df) <- str_replace(colnames(df),'-','_')
names(genelist) <- str_replace(names(genelist),'-','_')
cox_results <- lapply(names(genelist), function(celltype) {
  print(celltype)
  formula_obj <- as.formula(paste("Surv(time_os, os) ~ ", celltype, "+ cohort"))
  cox_model <- coxph(formula_obj, df)
  summary_cox <- summary(cox_model)
  # Extract HR and correctly calculate CI
  HR <- exp(cox_model$coefficients[celltype])  # Hazard Ratio
  HR_CI <- exp(confint(cox_model)[celltype, ])  # Exponentiate CI values
  HR_CI_Lower <- HR_CI[1]  # Lower bound
  HR_CI_Upper <- HR_CI[2]  # Upper bound
  p_value <- summary_cox$coefficients[celltype, "Pr(>|z|)"]
  return(data.frame(CellType = celltype, HR, HR_CI_Lower, HR_CI_Upper, p_value))
})
cox_results_df <- bind_rows(cox_results)

# Create a forest plot using ggplot2
ggplot(cox_results_df, aes(x = reorder(CellType, HR), y = HR)) +
  geom_point(color = "blue", size = 3) +  # Point for HR
  geom_errorbar(aes(ymin = HR_CI_Lower, ymax = HR_CI_Upper), width = 0.2, color = "black") +  # Error bars for CI
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Reference line at HR = 1
  coord_flip() +  # Flip coordinates for better readability
  labs(x = "Cell Type", y = "Hazard Ratio (HR)", title = "Forest Plot of Univariate Cox Regression") +
  theme_minimal(base_size = 14)

df <- comb |> 
  filter(tx_status == 'Baseline') |>
  filter(!response %in% c(NA, 'NE', 'MR')) |> 
  mutate(binary_response = ifelse(response %in% c(1, 'CR', 'CRPR', 'PR', 'R', 'Responder'), 1, 0))
colnames(df) <- str_replace(colnames(df),'-','_')
names(genelist) <- str_replace(names(genelist),'-','_')
logistic_results <- lapply(names(genelist), function(celltype) {
  print(celltype)
  df[, celltype] <- scale(df[, celltype])  # Optional: scale cell type abundances
  formula_obj <- as.formula(paste("binary_response ~ ", celltype, "+ cohort"))
  logit_model <- glm(formula_obj, data = df, family = binomial())
  # Extract odds ratio (OR) and CI
  OR <- exp(coef(logit_model)[celltype])
  OR_CI <- exp(confint(logit_model)[celltype, ])
  OR_CI_Lower <- OR_CI[1]
  OR_CI_Upper <- OR_CI[2]
  p_value <- coef(summary(logit_model))[celltype, "Pr(>|z|)"]
  return(data.frame(
    CellType = celltype,
    OR = OR,
    OR_CI_Lower = OR_CI_Lower,
    OR_CI_Upper = OR_CI_Upper,
    p_value = p_value
  ))
})
logistic_results <- bind_rows(logistic_results)

ggplot(logistic_results, aes(x = reorder(CellType, OR), y = OR)) +
  geom_point(color = "blue", size = 3) +  # Point for HR
  geom_errorbar(aes(ymin = OR_CI_Lower, ymax = OR_CI_Upper), width = 0.2, color = "black") +  # Error bars for CI
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Reference line at HR = 1
  coord_flip() +  # Flip coordinates for better readability
  labs(x = "Cell Type", y = "Hazard Ratio (HR)", title = "Forest Plot of Univariate Logistic Regression") +
  theme_minimal(base_size = 14)






