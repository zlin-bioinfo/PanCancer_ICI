rm(list=ls())
pkgs <- c('Seurat','tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','qs2', 'lmerTest','grid','msigdbr','ggplotify','pheatmap','msigdbr','MetBrewer','tibble','ComplexHeatmap','UCell')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
setwd("/home/zlin/workspace/PanCancer_ICI")
metadata <- read.csv('tables/meta_all.csv')
filter_sample <- metadata |> 
  distinct(sample, .keep_all = T) |> 
  pull(sample)
datasets <- c('SKCM_Becker','SKCM_Plozniak', 
              'BCC_Yost', 'SCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 
              'NSCLC_Yan', 'NSCLC_Liu',
              'PCa_Hawley','RCC_Bi','HCC_Guo','HCC_Ma')
# datasets <- c('SKCM_Plozniak', 
#               'BRCA_Bassez1', 
#               'HCC_Guo')

# # CD4T
# gs_list <- readxl::read_xlsx('tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 7, skip = 1) |> 
#   lapply(function(x){return(x[!is.na(x)])})
# gs_list$Exhaustion <- NULL
# names(gs_list)[names(gs_list  ) == 'OXPHOS'] = 'Oxidative phosphorylation'
# names(gs_list)[names(gs_list) == 'Lipid metabolism'] = 'Fatty acid metabolism'
# gs_list$`Fatty acid metabolism`[gs_list$`Fatty acid metabolism` == 'PLA2G16'] <- 'PLAAT3'
# gs_list$`Fatty acid metabolism`[gs_list$`Fatty acid metabolism` == 'ARNTL'] <- 'BMAL1'
# gs_list$`Oxidative phosphorylation` <- gsub('\\.', '-', gs_list$`Oxidative phosphorylation`)
# Differentiation <- c("Naïve", "Activation/Effector function")
# Function <- c("TCR signaling", "Cytotoxicity", "Cytokine/Cytokine receptor",
#               "Chemokine/Chemokine receptor", "Stress response", "Adhesion",
#               "IFN response", "Treg signature", "Costimulatory molecules")
# Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
# Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")

# Gene sets for CD4+ T cells
gs_list_cd4 <- list(
  Naïve = c("IL7R", "CCR7", "SELL", "FOXP1", "KLF2", "KLF3", "LEF1", "TCF7", "ACTN1", "BTG1", "BTG2", "TOB1"),
  `Activation/Effector Function` = c("FAS", "CD44", "CD69", "CD38", "NKG7", "KLRB1", "KLRD1", "KLRG1", "CX3CR1", "CD300A", "FGFBP2", "ID2", "ID3", "PRDM1", "RUNX3", "TBX21", 
                                     "ZEB2", "BATF", "NR4A1", "NR4A2", "HOPX", "FOS", "FOSB", "FOSL2", "JUN", "JUNB", "JUND", "STAT1", "STAT3", "EOMES", "AHR"),
  `TCR Signaling` = c("CALM1", "CALM2", "CALM3", "CD4", "CAST", "CD247", "CD3D", "CD3E", "CD3G", "CSK", "DOK2", "FYN", "LCK", "NFATC1", 
                      "NFATC2", "PLEK", "PTPN11", "PTPN13", "PTPN2", "PTPN22", "PTPN4", "PTPN6", "PTPN7", "PTPRC", "PTPRCAP", "S100A10", "S100A11", "S100A4", 
                      "S100A6", "ZAP70", "DUSP1", "DUSP2", "DUSP4", "DUSP16", "LAT", "FOS", "FOSB", "FOSL2", "JUN", "JUNB", "JUND", "NR4A1", 
                      "NR4A2", "BATF", "IRF1", "SH2D1A", "SH2D2A", "MAP2K3", "MAP3K4", "MAP3K8", "MAP4K1", "NFKB2", "NFKBIA", "NFKBIZ", "REL", "RELB"),
  Cytotoxicity = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "GNLY", "PRF1", "IFNG", "TNF", "SERPINB9", "CTSA", "CTSB", "CTSC", 
                   "CTSD", "CTSH", "CTSW", "CST7", "CAPN2", "PLEK"),
  `Cytokine/Cytokine Receptor` = c("CSF1", "CSF2", "ADAM19", "ADAM8", "ADAM12", "CD70", "IL12RB2", "IL17A", "IL17F", "IL10RA", "IL1R1", "IL1R2", "IL21R", "IL21", 
                                   "IL26", "IL2RA", "IL2RB", "IL2RG", "IL32", "IL6R", "TGFB1", "TGFBR2", "TGFBR3"),
  `Chemokine/Chemokine Receptor` = c("CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CXCR3", "CXCR4", "CXCR5", "CXCR6", "CCL3", "CCL4", "CCL5", "CCL20", "CXCL13", "CXCL8", "XCL1"),
  Adhesion = c("ITGA6", "ITGA4", "ITGAE", "ITGAL", "ITGAM", "ITGB1", "ITGB2", "ITGB7", "SELL", "SELPLG", "S1PR1", "ICAM2"),
  `IFN Response` = c("STAT1", "STAT3", "MX1", "IRF1", "ISG15", "ISG20", "IFITM1", "IFITM2", "IFITM3", "OAS1", "OAS2", "OASL", "SOCS1", "SOCS3", "TRIM22", "APOL6", 
                     "IFNAR2", "IFNGR1", "GBP1", "GBP2", "GBP4", "GBP5", "BST2", "IFI16", "IFI35", "IFI44L", "IFI6", "PARP8", "PARP9"),
  `Treg Signature` = c("FOXP3", "IKZF2", "IKZF4", "IL2RA", "ENTPD1", "CCR4", "ICOS", "IL10RA", "TGFB1", "TIGIT", "CTLA4", "LAG3", "HAVCR2", "PDCD1"),
  `Costimulatory Molecules` = c("TNFRSF25", "TNFRSF1B", "TNFRSF4", "TNFRSF9", "TNFRSF18", "CD27", "CD28", "CD44", "CD48", "ICOS", "CD2", "SLAMF1", "CD40LG", "CD84"),
  Autophagy = c("BECN1", "ULK1", "AMBRA1", "UVRAG", "ATG3", "ATG5", "ATG7", "ATG9A", "ATG9B", "ATG12", "ATG13", "ATG16L1", 
                "MAP1LC3A", "MAP1LC3B", "MAP1LC3C", "WIPI1", "WIPI2", "ATG4A", "ATG4B", "ATG4C", "ATG4D", "ATG14"),
  `MAPK Signaling` = c("MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK10", "MAPK11", "MAPK14", 
                       "MAP3K5", "MAP2K3", "MAP2K4", "MAP2K6", "MAP2K7"),
  `Oxidative Phosphorylation` = c("ATP1B1", "ATP1B3", "ATP2B1", "ATP2B4", "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5E", "ATP5F1", "ATP5F1A", "ATP5F1B", "ATP5F1C", 
                                  "ATP5F1D", "ATP5F1E", "ATP5G1", "ATP5G2", "ATP5G3", "ATP5H", "ATP5I", "ATP5IF1", "ATP5J", "ATP5J2", "ATP5MG", "ATP5MC2", "ATP5MC3", 
                                  "ATP5MD", "ATP5MF", "ATP5MG", "ATP5MPL", "ATP5O", "ATP5PB", "ATP5PD", "ATP5PF", "ATP6V0B", "ATP6V0C", "ATPIF1", "COX16", "COX17", 
                                  "COX5A", "COX6A1", "COX6C", "COX8A", "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", 
                                  "MT-ND4L", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "NDUFA13", "NDUFA6", "NDUFB2", 
                                  "NDUFC1", "NDUFV2", "UQCR10", "UQCR11", "UQCR11-1", "UQCRB", "UQCRC2", "UQCRHL", "SURF4", "CYCS", "ISCU", "MTCH1"),
  Glycolysis = c("SLC2A3", "SLC16A1", "PFKFB3", "ALDOA", "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA", "LDHB", "GYG1", "GPI", "UGP2", "UBB", "UBC"),
  `Fatty Acid Metabolism` = c("ABHD3", "ACADVL", "ACOT9", "AHR", "ALOX5AP", "AP2M1", "AP2S1", "APOBEC3C", "APOBEC3G", "APOBEC3H", "APOL6", "ARF1", "ARF5", 
                              "BMAL1", "B3GALT2", "B4GALT1", "CAV1", "CDIPT", "CERK", "CLTA", "CPNE7", "ELOVL1", "ELOVL5", "DBI", "EPHX2", "FABP5", 
                              "FURIN", "GDE1", "GK", "GPX1", "GPX4", "HDLBP", "IDH2", "IDI1", "HPGD", "INPP4B", "INPP5D", "INPP5F", "INSIG1", 
                              "KDSR", "LDLRAP1", "LPCAT1", "MFSD2A", "NCOA3", "NDUFAB1", "NEU1", "NPC2", "NUDT7", "P4HB", "PLAAT3", "PLIN2", "PPP1CA", 
                              "PPP1CB", "PPP1CC", "PTGES3", "SAR1B", "SARDH", "SC5D", "SGMS1", "SH3KBP1", "SREBF2", "STARD7", "TNFAIP8"),
  `Pro-Apoptosis` = c("BAX", "BAG3", "CASP1", "CASP4", "CYCS", "BCL2L11"),
  `Anti-Apoptosis` = c("BIN1", "BIN2", "BIRC3", "BCL2", "BCL2L1", "MCL1"),
  `Oxidative Stress` = c("SOD1", "SOD2", "SOD3", "GPX1", "GPX2", "GPX3", "GPX5", "GPX6", "GPX7", "GPX8", 
                         "PRDX1", "PRDX2", "PRDX3", "PRDX5", "PRDX6", "TXN", "TXN2", "TXNRD1", "TXNRD2", 
                         "CAT", "GSR", "GSTP1", "NOX4", "NOX5", "NCF1", "NCF2", "NCF4"),
  `Heat Shock Response` = c("HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "HSPA7", 
                            "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14", 
                            "HSP90AA1", "HSP90AB1", "DNAJB1", "DNAJB6", "DNAJC2", "DNAJC7", "HSPH1"),
  `DNA Damage Response` = c("ATM", "ATR", "MRE11", "NBN", "RAD50", "TP53", "RB1", 
                            "CDK2", "CDK4", "CDK6", "CCNA1", "CCNA2", "CCNE1", "CCNE2", 
                            "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", 
                            "TERF1", "TERF2", "POT1", "TINF2", "TERF2IP", "ACD", 
                            "HIRA", "ASF1A", "KAT5", "EP400", "UBN1", "CABIN1"),
  `Hypoxia Response` = c("HIF1A", "HIF1AN", "HIF3A", "EPAS1", "EGLN1", "EGLN2", "EGLN3", 
                         "VEGFA", "VHL", "MTOR", "RPTOR", "TSC1", "TSC2", "PRKAA1", "PRKAA2", "PRKAG1", "PRKAG2", "PRKAG3"),
  `Immune Stress Response` = c("IL1A", "IL6", "CXCL8", "IFNB1", "NFKB1", "RELA", "MDM2", "MDM4"),
  `ER Stress Response` = c("HSPA5", "P4HB", "ERO1A", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G")
)

# Categories for CD4+ T cells
Differentiation <- c("Naïve", "Activation/Effector Function")
Function <- c("TCR Signaling", "Cytotoxicity", "Cytokine/Cytokine Receptor", "Chemokine/Chemokine Receptor", "Adhesion", "IFN Response", "Treg Signature", "Costimulatory Molecules","MAPK Signaling","Autophagy")
Metabolism <- c("Oxidative Phosphorylation", "Glycolysis", "Fatty Acid Metabolism")
Apoptosis <- c("Pro-Apoptosis", "Anti-Apoptosis")
StressResponse <- c("Oxidative Stress", "Heat Shock Response", "DNA Damage Response", "Hypoxia Response", "Immune Stress Response", "ER Stress Response")

list_cd4 <- lapply(datasets, function(dataset) {
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = celltype_main == 'CD4+T') |> 
    subset(subset = sample %in% filter_sample) |> 
    NormalizeData() |> 
    AverageExpression(group.by = 'celltype_r2', layer = 'data', return.seurat = T) |> 
    AddModuleScore_UCell(features = gs_list_cd4)
  for(i in 1:length(gs_list_cd4)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list_cd4)[i], "_UCell")] <- names(gs_list_cd4)[i]
  }
  df_score <- seu@meta.data[,c('celltype_r2', names(gs_list_cd4))]
  df_score$dataset <- dataset
  return(df_score)
})

df_score <- do.call(rbind, list_cd4)
write.csv(df_score, 'tables/cd4_func.csv')
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis, StressResponse)
col_order <- c('CD4-T-naive','CD4-Tcm','CD4-Tctl','CD4-T-ISG','CD4-Tfh','CD4-Th17','CD4-Treg','CD4-Tstr')
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(df_score$celltype_r2)),
                              nrow = length(MarkerNameVector))
colnames(FunctionScoreMatrix) <- unique(df_score$celltype_r2)
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)) {
  for(ri in 1:nrow(FunctionScoreMatrix)) {
    FunctionVec <- as_tibble(df_score) |> pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[df_score$celltype_r2 == unique(df_score$celltype_r2)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, scales::rescale, to=c(0, 1)))
FunctionScoreMatrix <- FunctionScoreMatrix[, col_order]
cols <- colorRampPalette(brewer.pal('YlGnBu', n=8))(10)
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis)),
  rep("Stress Response", length(StressResponse))))
signatureType_row$Signature.type <- factor(signatureType_row$Signature.type, levels = c('Differentiation', 'Function','Metabolism','Apoptosis','Stress Response'))
pdf('figures/Functional_score_T/ht_func_cd4.pdf', height = 7, width = 7)
Heatmap(FunctionScoreMatrix, name = 'Signature score \n(scaled)',
        column_title = 'CD4+ T cells',
        cluster_rows = F, cluster_columns = F,
        column_names_rot = 45,
        col = cols, 
        row_names_gp = gpar(fontsize = 10), 
        column_names_gp = gpar(fontsize = 9), 
        width = ncol(FunctionScoreMatrix)*unit(6, "mm"), 
        height = nrow(FunctionScoreMatrix)*unit(6, "mm"),
        row_title_gp = gpar(fontsize = 7.5, fontface = 'bold'),
        rect_gp = gpar(col = "white", lwd = 1),
        heatmap_legend_param = list(
          legend_direction = "horizontal",
          legend_width = unit(2, "cm"), at = c(0, 1),labels = c("Min",  "Max"),
          legend_side = 'bottom',
          title_position = "topcenter"),
        row_split = signatureType_row$Signature.type)
dev.off()

# CD8T
# gs_list <- readxl::read_xlsx('tables/41591_2023_2371_MOESM3_ESM.xlsx', sheet = 5, skip = 1) |> 
#   lapply(function(x){return(x[!is.na(x)])})
# names(gs_list)[names(gs_list) == 'Activation:Effector function'] <- 'Activation/Effector function'
# names(gs_list)[names(gs_list) == 'TCR Signaling'] = 'TCR signaling'
# names(gs_list)[names(gs_list) == 'IFN Response'] = 'IFN response'
# gs_list$`IFN response`[gs_list$`IFN response` == 'DDX58'] <- 'RIGI'
# gs_list$`Oxidative phosphorylation` <- gsub('\\.', '-', gs_list$`Oxidative phosphorylation`)
# gs_list <- gs_list[-which(names(gs_list)=='Senescence')]

# Gene sets for CD8+ T cells
gs_list_cd8 <- list(
  Naïve = c("IL7R", "CCR7", "SELL", "FOXO1", "KLF2", "KLF3", "LEF1", "TCF7", "ACTN1", "FOXP1"),
  `Activation/Effector Function` = c("FAS", "FASLG", "CD44", "CD69", "CD38", "NKG7", "KLRB1", "KLRD1", "KLRF1", "KLRG1", "KLRK1", "FCGR3A", "CX3CR1", "CD300A", "FGFBP2", "ID2", 
                                     "ID3", "PRDM1", "RUNX3", "TBX21", "ZEB2", "BATF", "IRF4", "NR4A1", "NR4A2", "NR4A3", "PBX3", "ZNF683", "HOPX", "FOS", "FOSB", "JUN", 
                                     "JUNB", "JUND", "STAT1", "STAT2", "STAT5A", "STAT6", "STAT4", "EOMES"),
  Exhaustion = c("PDCD1", "LAYN", "HAVCR2", "LAG3", "CD244", "CTLA4", "LILRB1", "TIGIT", "TOX", "VSIR", "BTLA", "ENTPD1", "CD160", "LAIR1"),
  `TCR Signaling` = c("CALM1", "CALM2", "CALM3", "CAST", "CD247", "CD3D", "CD3E", "CD3G", "CSK", "DOK1", "DOK2", "FYN", "LCK", "NFATC2", 
                      "NFATC1", "NFATC4", "NFATC3", "PLEK", "PAG1", "PTPN11", "PTPN2", "PTPN22", "PTPN4", "PTPN6", "PTPN7", "PTPRC", "PTPRCAP", "S100A10", 
                      "S100A11", "S100A13", "S100A4", "S100A6", "ZAP70", "DUSP1", "DUSP2", "DUSP4", "DUSP5", "DUSP10", "LAT", "PLCG1", "PLCG2", "PPP3CA", 
                      "PPP3CC", "FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "NR4A1", "NR4A2", "NR4A3", "BATF", "IRF4", "SH2D2A"),
  Cytotoxicity = c("GZMA", "GZMB", "GZMH", "GZMK", "GZMH", "GNLY", "PRF1", "IFNG", "TNF", "SERPINB1", "SERPINB6", "SERPINB9", "CTSA", 
                   "CTSB", "CTSC", "CTSD", "CTSW", "CST3", "CST7", "CSTB", "LAMP1", "LAMP3", "CAPN2"),
  `Cytokine/Cytokine Receptor` = c("CSF1", "IL10RA", "IL16", "IL17RA", "IL18RAP", "IL21R", "IL2RB", "IL2RG", "IL32", "IL9R", "ADAM10", "ADAM8", "METRNL", "CD70"),
  `Chemokine/Chemokine Receptor` = c("CCR4", "CCR5", "CCR7", "CXCR3", "CXCR4", "CXCR5", "CXCR6", "CCL3", "CCL4", "CCL4L1", "CCL4L2", "CCL5", "CXCL13", "CXCL8", "XCL1", "XCL2"),
  Anergy = c("NT5E", "IZUMO1R", "LAG3", "NRP1", "DGKA", "CBLB", "RNF128", "ITCH", "NFATC2", "EGR2", "EGR3", "NR4A1", "TOB1"),
  `NFKB Signaling` = c("NFKB1", "NFKB2", "NFKBIA", "NFKBIB", "NFKBIZ", "CHUK", "IKBKB", "IKBKG", "REL", "RELA", "RELB"),
  Adhesion = c("ITGA1", "ITGA4", "ITGAE", "ITGAL", "ITGAM", "ITGB1", "ITGB2", "ITGB7", "SELL", "SELPLG", "S1PR1", "VCAM1", "ICAM2", "ICAM3"),
  `IFN Response` = c("IFIT1", "IFIT2", "IFIT3", "IFIT5", "STAT1", "STAT2", "MX1", "IRF1", "IRF4", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "IFITM1", "IFITM2", 
                     "IFITM3", "OAS1", "OAS2", "OAS3", "JAK1", "JAK2", "SOCS1", "SOCS3", "TRIM14", "TRIM21", "TRIM22", "APOL1", "APOL2", "APOL6", "IFNGR1", "GBP1", 
                     "GBP2", "GBP4", "GBP5", "GBP3", "BST2", "CMPK2", "RIGI", "DDX60", "DDX60L", "IFI30", "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1", "PARP10", 
                     "PARP12", "PARP14"),
  Autophagy = c("BECN1", "ULK1", "AMBRA1", "UVRAG", "ATG3", "ATG5", "ATG7", "ATG9A", "ATG9B", "ATG12", "ATG13", "ATG16L1", 
                "MAP1LC3A", "MAP1LC3B", "MAP1LC3C", "WIPI1", "WIPI2", "ATG4A", "ATG4B", "ATG4C", "ATG4D", "ATG14"),
  `MAPK Signaling` = c("MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK10", "MAPK11", "MAPK14", 
                       "MAP3K5", "MAP2K3", "MAP2K4", "MAP2K6", "MAP2K7"),
  `Oxidative Phosphorylation` = c("ATP1B3", "ATP2A3", "ATP2B1", "ATP2B4", "ATP5A1", "ATP5B", "ATP5C1", "ATP5D", "ATP5E", "ATP5EP2", "ATP5F1", "ATP5F1A", "ATP5F1B", 
                                  "ATP5F1C", "ATP5F1D", "ATP5F1E", "ATP5G2", "ATP5G3", "ATP5H", "ATP5I", "ATP5IF1", "ATP5J", "ATP5J2", "ATP5L", "ATP5MC1", "ATP5MC2", 
                                  "ATP5MC3", "ATP5MD", "ATP5ME", "ATP5MF", "ATP5MG", "ATP5MPL", "ATP5O", "ATP5PB", "ATP5PD", "ATP5PF", "ATP6AP1", "ATP6AP2", "ATP6V0B", 
                                  "ATP6V0C", "ATP6V0D1", "ATP6V0E1", "ATP6V0E2", "ATP6V1D", "ATP6V1F", "ATP6V1G1", "ATP8B2", "ATP8B4", "COX16", "COX17", "COX4I1", 
                                  "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7B", "COX7C", "COX8A", "CYC1", "MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", 
                                  "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "NDUFA1", "NDUFA10", "NDUFA11", "NDUFA12", 
                                  "NDUFA13", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFAB1", "NDUFAF3", "NDUFAF4", "NDUFAF8", "NDUFB10", "NDUFB11", 
                                  "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFC1", "NDUFC2", "NDUFS2", "NDUFS3", "NDUFS5", "NDUFS6", 
                                  "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "SDHA", "SDHB", "TCIRG1", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ", "SURF1", 
                                  "SURF4"),
  Glycolysis = c("SLC2A3", "SLC2A8", "PFKFB3", "ALDOA", "ENO1", "GAPDH", "GPI", "PGAM1", "PGK1", "PKM", "TPI1", "GYG1", "MDH1", "MDH2", "LDHA", "LDHB", "ADH5", "IDH2"),
  `Fatty Acid Metabolism` = c("ACAA1", "ACADM", "ACADVL", "ACSL5", "ALDH16A1", "ALDH9A1", "ECH1", "ECHS1", "HADHA", "HSD17B10", "SLC27A2", "DECR1", "PEX16", "SCP2"),
  `Pro-Apoptosis` = c("BCL2L11", "BID", "BIN1", "BIN2", "BAX", "BAK1", "CASP8", "CASP3"),
  `Anti-Apoptosis` = c("BCL11B", "BCL2L1", "BAG3", "BAG4", "BIRC3", "GADD45B"),
  `Oxidative Stress` = c("SOD1", "PRDX1", "PRDX2", "PRDX3", "PRDX5", "PRDX6", "TXN", "TXN2", "TXNRD1", "TXNRD2", "CAT", "GSR", "GSTP1", "NOX4", "NOX5", "NCF1", "NCF2", "NCF4"),
  `Heat Shock Response` = c("HSPA1A", "HSPA1B", "HSPA4", "HSPA5", "HSPA6", "HSPA8", "HSP90AA1", "HSP90AB1", "DNAJB1", "DNAJB6", "HSPH1"),
  `DNA Damage Response` = c("ATM", "ATR", "TP53", "RB1", "CDK4", "CDK6", "CDKN1A", "CDKN2C", "CDKN2D", "TERF1", "TERF2IP", "TINF2", "HIRA", "ASF1A", "KAT5", "EP400", "UBN1", "CABIN1"),
  `Hypoxia Response` = c("HIF1A", "VHL", "MTOR", "RPTOR", "TSC1", "TSC2", "PRKAA1", "PRKAA2", "PRKAG1", "PRKAG2", "PRKAG3"),
  `Immune Stress Response` = c("IL6", "CXCL8", "NFKB1", "RELA", "MDM4"),
  `ER Stress Response` = c("HSPA5", "P4HB", "ERO1A", "CAMK2G")
)

# Categories for CD8+ T cells
Differentiation <- c("Naïve", "Activation/Effector Function", "Exhaustion")
Function <- c("TCR Signaling", "Cytotoxicity", "Cytokine/Cytokine Receptor", "Chemokine/Chemokine Receptor", "Adhesion", "IFN Response", "Anergy", "NFKB Signaling", "Autophagy", "MAPK Signaling")
Metabolism <- c("Oxidative Phosphorylation", "Glycolysis", "Fatty Acid Metabolism")
Apoptosis <- c("Pro-Apoptosis", "Anti-Apoptosis")
StressResponse <- c("Oxidative Stress", "Heat Shock Response", "DNA Damage Response", "Hypoxia Response", "Immune Stress Response", "ER Stress Response")

cd8t <- c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',"CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',"CD8_NK-like")
list_cd8 <- lapply(datasets, function(dataset) {
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = celltype_r2 %in% cd8t) |> 
    subset(subset = sample %in% filter_sample) |> 
    NormalizeData() |> 
    AverageExpression(group.by = 'celltype_r2', layer = 'data', return.seurat = T) |> 
    AddModuleScore_UCell(features = gs_list_cd8)
  for(i in 1:length(gs_list_cd8)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list_cd8)[i], "_UCell")] <- names(gs_list_cd8)[i]
  }
  df_score <- seu@meta.data[,c('celltype_r2', names(gs_list_cd8))]
  df_score$dataset <- dataset
  return(df_score)
})

df_score <- do.call(rbind, list_cd8)
write.csv(df_score, 'tables/cd8_func.csv')
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis, StressResponse)
col_order <- c('CD8-T-naive','CD8-Tm','CD8-Trm','CD8-Tem-early','CD8-Tem','CD8-Tpex',"CD8-Tex-CXCL13", "CD8-Tex-GZMK",'CD8-Temra','CD8-T-ISG','CD8-Tstr',"CD8-NK-like")
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(df_score$celltype_r2)),
                              nrow = length(MarkerNameVector))
colnames(FunctionScoreMatrix) <- unique(df_score$celltype_r2)
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)) {
  for(ri in 1:nrow(FunctionScoreMatrix)) {
    FunctionVec <- as_tibble(df_score) |> pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[df_score$celltype_r2 == unique(df_score$celltype_r2)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, scales::rescale, to=c(0, 1)))
FunctionScoreMatrix <- FunctionScoreMatrix[, col_order]
cols <- colorRampPalette(brewer.pal('YlGnBu', n=8))(10)
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis)),
  rep("Stress Response", length(StressResponse))))
signatureType_row$Signature.type <- factor(signatureType_row$Signature.type, levels = c('Differentiation', 'Function','Metabolism','Apoptosis','Stress Response'))
pdf('figures/Functional_score_T/ht_func_cd8.pdf', height = 7, width = 7)
Heatmap(FunctionScoreMatrix, name = 'Signature score \n(scaled)',
        column_title = 'CD8+ T cells',
        cluster_rows = F, cluster_columns = F,
        column_names_rot = 45,
        col = cols, 
        row_names_gp = gpar(fontsize = 10), 
        column_names_gp = gpar(fontsize = 9), 
        width = ncol(FunctionScoreMatrix)*unit(6, "mm"), 
        height = nrow(FunctionScoreMatrix)*unit(6, "mm"),
        row_title_gp = gpar(fontsize = 7.5, fontface = 'bold'),
        rect_gp = gpar(col = "white", lwd = 1),
        heatmap_legend_param = list(
          legend_direction = "horizontal",
          legend_width = unit(2, "cm"), at = c(0, 1),labels = c("Min",  "Max"),
          legend_side = 'bottom',
          title_position = "topcenter"),
        row_split = signatureType_row$Signature.type)
dev.off()

# Dataset names
datasets <- c('SKCM_Becker','SKCM_Plozniak','BCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2', 'TNBC_Zhang', 'TNBC_Shiao',
              'HNSC_Franken', 'HNSC_vanderLeun', 'HNSC_Luoma', 
              'CRC_Li', 'CRC_Chen', 'NSCLC_Yan', 
              'PCa_Hawley','RCC_Bi','HCC_Guo')
# Macrophages
ifna <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_INTERFERON_ALPHA_RESPONSE') |>
  pull(gene_symbol) |>
  unique()
ifng <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_INTERFERON_GAMMA_RESPONSE') |>
  pull(gene_symbol) |>
  unique()
ifng[ifng == 'DDX58'] <- 'RIGI'
mhc1 <- msigdbr(species = "Homo sapiens", category = 'C2') |> 
  filter(gs_name  == 'REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION') |> 
  pull(gene_symbol) |> 
  unique()
mhc2 <- msigdbr(species = "Homo sapiens", category = 'C2') |> 
  filter(gs_name  == 'REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION') |> 
  pull(gene_symbol) |> 
  unique()
complement <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_COMPLEMENT') |> 
  pull(gene_symbol) |> 
  unique()
tnfa <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_TNFA_SIGNALING_VIA_NFKB') |> 
  pull(gene_symbol) |> 
  unique()
tnfa[tnfa == 'DDX58'] <- 'RIGI'
tgfb <- msigdbr(species = "Homo sapiens", category = 'H') |> 
  filter(gs_name  == 'HALLMARK_TGF_BETA_SIGNALING') |> 
  pull(gene_symbol) |> 
  unique()
infla <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_INFLAMMATORY_RESPONSE') |>
  pull(gene_symbol) |>
  unique()
mtorc1 <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_MTORC1_SIGNALING') |>
  pull(gene_symbol) |>
  unique()
ros <- msigdbr(species = "Homo sapiens", category = 'H') |>
  filter(gs_name  == 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY') |>
  pull(gene_symbol) |>
  unique()
ecm <- msigdbr(species = "Homo sapiens", category = 'C2') |>
  filter(gs_name  == 'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION') |>
  pull(gene_symbol) |>
  unique()
pd1 <- genesets <- msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME") |> 
  filter(gs_name  == 'REACTOME_PD_1_SIGNALING') |>
  pull(gene_symbol) |>
  unique()
sig_mac <- readxl::read_xlsx('tables/1-s2.0-S0092867421000106-mmc5.xlsx', skip = 1)
sig_mac$M1[sig_mac$M1 == 'IL23'] <- 'IL23A'
sig_mac$M2[sig_mac$M2 == 'FASL'] <- 'FASLG'
gs_list <- list(`M1 Signature` = sig_mac$M1[!is.na(sig_mac$M1)],
                `M2 Signature` = sig_mac$M2[!is.na(sig_mac$M2)],
                Angiogenesis = sig_mac$Angiogenesis[!is.na(sig_mac$Angiogenesis)],
                Phagocytosis = sig_mac$Phagocytosis[!is.na(sig_mac$Phagocytosis)],
                Inflammatory = infla,
                `TNF-alpha Signaling` = tnfa,
                `TGF-beta Signaling` = tgfb,
                `IFN-alpha Response` = ifna,
                `IFN-gamma Response` = ifng,
                Complement = complement,
                `MHC I Antigen Proc.&Pres.` = mhc1,
                `MHC II Antigen Presentation` = mhc2,
                `mTORC1 Signaling` = mtorc1,
                `Reactive Oxigen Species` = ros,
                `Extracellular matrix organization` = ecm,
                `PD1 Signaling` = pd1)

list_mac <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = celltype_main %in% c('Mono/macro', 'Mono', 'Macro')) |> 
    subset(subset = sample %in% filter_sample) |> 
    NormalizeData() |> 
    AverageExpression(group.by = 'celltype_r2', layer = 'data', return.seurat = T) |> 
    AddModuleScore_UCell(features = gs_list)
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  df_score <- seu@meta.data[,c('celltype_r2', names(gs_list))]
  df_score$dataset <- dataset
  return(df_score)
})

df_score <- do.call(rbind, list_mac)
write.csv(df_score, 'tables/monomac_func.csv')
# Heatmap
pdf('figures/Functional_score_Mac/ht_func_mac.pdf', height = 6, width = 6)
df_score |> 
  group_by(celltype_r2) |> 
  summarise(across(where(is.numeric), median, na.rm = TRUE)) |> 
  column_to_rownames(var = 'celltype_r2') |> 
  apply(2, scales::rescale, to=c(0, 1)) |> t() |>
  Heatmap(column_names_rot = 45, name = 'Signature score \n(scaled)', column_title = 'Monocyte/macrophages',
          heatmap_legend_param = list(
            legend_direction = "horizontal",
            legend_width = unit(2, "cm"), at = c(0, 1),labels = c('Min','Max'),
            title_position = "topcenter"),
          row_names_gp = gpar(fontsize = 10), 
          rect_gp = gpar(col = "white", lwd = 0.5),
          column_names_gp = gpar(fontsize = 9), 
          col = colorRampPalette(brewer.pal(10, "YlGnBu"))(10),
          width = 9*unit(6, "mm"), 
          height = 12*unit(6, "mm"))
dev.off()

# Dataset names
caf_pathways <- list(
  # Core shared pathways with annotations
  Shared = c(
    "GOBP_TISSUE_DEVELOPMENT",                         # Pro-tumor: stromal support
    "GOBP_VASCULATURE_DEVELOPMENT",                    # Pro-tumor: angiogenesis
    "GOBP_CELL_MIGRATION",                             # Pro-tumor: invasion
    "GOBP_REGULATION_OF_CELL_POPULATION_PROLIFERATION", # Pro-tumor: growth
    "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX" # Pro-tumor: ECM; Anti-immune: barrier
  ),
  
  # Distinct pathways with annotations
  iCAF_MMP1 = c(
    "GOBP_CYTOKINE_PRODUCTION",               # Pro-tumor, Pro-immune, Anti-immune
    "GOBP_INNATE_IMMUNE_RESPONSE"            # Pro-tumor, Pro-immune, Anti-immune
  ),
  Myofibroblasts = c(
    "GOBP_BLOOD_VESSEL_MORPHOGENESIS",                # Pro-tumor: vascularization
    "GOBP_REGULATION_OF_CELLULAR_COMPONENT_MOVEMENT"   # Pro-tumor: migration regulation
  ),
  iCAF_IL6 = c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",       # Pro-tumor, Pro-immune, Anti-immune
    "GOBP_RESPONSE_TO_LIPID"                          # Pro-tumor: lipid signaling
  ),
  CAF_desmo = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",      # Pro-tumor: EMT
    "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT" # Pro-tumor, Anti-immune
  ),
  CAF_SFRP2 = c(
    "GOBP_ENZYME_LINKED_RECEPTOR_PROTEIN_SIGNALING_PATHWAY", # Pro-tumor: kinase signaling
    "GOBP_ANIMAL_ORGAN_MORPHOGENESIS"                 # Pro-tumor: structural support
  ),
  CAF_ap = c(
    "GOBP_DEFENSE_RESPONSE",                  # Pro-tumor, Pro-immune, Anti-immune
    "GOBP_REGULATION_OF_IMMUNE_SYSTEM_PROCESS" # Pro-tumor, Pro-immune, Anti-immune
  ),
  CAF_prog = c(
    "GOBP_POSITIVE_REGULATION_OF_DEVELOPMENTAL_PROCESS",   # Pro-tumor: developmental support
    "GOBP_NEGATIVE_REGULATION_OF_RESPONSE_TO_STIMULUS" # Pro-tumor, Anti-immune
  )
)
pathways <- unlist(tme_caf_pathways) |> unique()
datasets <- c('SKCM_Becker','SKCM_Plozniak','BCC_Yost',
              'BRCA_Bassez1', 'BRCA_Bassez2',  'TNBC_Shiao',
              'HNSC_Franken', 
              'CRC_Li', 'CRC_Chen', 'NSCLC_Yan', 
              'PCa_Hawley','RCC_Bi','HCC_Ma')
go <- msigdbr(species = "Homo sapiens", category = 'C5') |> filter(gs_name %in% pathways) |> select(gs_name, gene_symbol)
go <- split(go$gene_symbol, go$gs_name)
hmk <- msigdbr(species = "Homo sapiens", category = 'H') |> filter(gs_name %in% pathways) |> select(gs_name, gene_symbol)
hmk <- split(hmk$gene_symbol, hmk$gs_name)
gs_list <- c(go, hmk)
list_caf <- lapply(datasets, function(dataset){
  print(dataset)
  seu <- qs_read(paste0('data/', dataset, '/seu_final.qs2')) |> 
    subset(subset = celltype_main == 'CAF') |> 
    NormalizeData() |> 
    AverageExpression(group.by = 'celltype_r2', layer = 'data', return.seurat = T) |> 
    AddModuleScore_UCell(features = gs_list, maxRank = 3000)
  seu$dataset <- dataset
  for(i in 1:length(gs_list)) {
    colnames(seu@meta.data)[colnames(seu@meta.data) == paste0(names(gs_list)[i], "_UCell")] <- names(gs_list)[i]
  }
  df_score <- seu@meta.data[,c('celltype_r2', names(gs_list))]
  df_score$dataset <- dataset
  return(df_score)
})

df_score <- do.call(rbind, list_caf)
write.csv(df_score, 'tables/caf_func.csv')
pdf('figures/Functional_score_Mac/ht_func_caf.pdf', height = 5, width = 10)
df_score |> 
  group_by(celltype_r2) |> 
  summarise(across(where(is.numeric), median, na.rm = TRUE)) |> 
  column_to_rownames(var = 'celltype_r2') |> 
  apply(2, scales::rescale, to=c(0, 1))  |>
  ComplexHeatmap::Heatmap(column_names_rot = 45, name = 'Signature score \n(scaled)', column_title = 'CAF', 
                          row_names_side = 'left',show_row_dend = F,
                          heatmap_legend_param = list(
                            # legend_direction = "horizontal",
                            legend_width = unit(2, "cm"), at = c(0, 1),labels = c('Min','Max'),
                            title_position = "topcenter"),
                          row_names_gp = gpar(fontsize = 10), 
                          rect_gp = gpar(col = "white", lwd = 0.5),
                          column_names_gp = gpar(fontsize = 8), 
                          col = colorRampPalette(brewer.pal(10, "YlGnBu"))(10),
                          width = 18*unit(6, "mm"), 
                          height = 6*unit(6, "mm"))
dev.off()

