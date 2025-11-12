t_nk <- c('CD4_T-naive','CD4_Tcm','CD4_Treg','CD4_T-ISG','CD4_Tfh','CD4_Tstr','CD4_Tctl','CD4_Th17',
          'CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
          "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
          "CD8_NK-like",'CD8_Tc17','MAIT','gdT','NK_CD56loCD16hi','NK_CD56hiCD16lo','Cycling T','Cycling NK')
bplasma <- c("B-naive", "B-ISG", "B-HSP", "B_MT2A", 
             "ACB_EGR1", "ACB_NR4A2", "ACB_CCR7", "B-memory", "B-AtM", 
             "GCB-pre", "GCB-DZ_SUGCT", "GCB-LZ_LMO2",
             "GCB-cycling", "PC-cycling","PC-trans",
             "PC-early_RGS13", "PC_IGHG", "PC_IGHA")
mye <- c('Mast','pDC','cDC1', 
         'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'MoDC', 'mregDC', 'cDC2', 
         'Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
         'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro-ISG', 
         'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2','Cycling myeloids')
nonimmune <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein",
               "Pericytes","SMC", "Myofibroblasts", "CAF_SFRP2", 
               "CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap", 'Cycling non-immune')
immune <- c(t_nk, bplasma, mye)

cd4t <- c('CD4_T-naive','CD4_Tcm','CD4_Treg','CD4_T-ISG','CD4_Tfh','CD4_Tstr','CD4_Tctl','CD4_Th17')
cd8t <- c('CD8_T-naive','CD8_Tm','CD8_Trm','CD8_Tem-early','CD8_Tem','CD8_Tpex',
          "CD8_Tex_CXCL13", "CD8_Tex_GZMK",'CD8_Temra','CD8_T-ISG','CD8_Tstr',
          "CD8_NK-like",'CD8_Tc17')
nk <- c('NK_CD56loCD16hi','NK_CD56hiCD16lo')
cdc <- c('cDC1', 'cDC2_CD1C', 'cDC2_IL1B','cDC2-ISG', 'cDC2_CXCL9', 'DC_LC-like', 'MoDC', 'mregDC', 'cDC2')
monomac <- c('Mono_CD14', 'Mono_CD14CD16', 'Mono_CD16',
             'Macro_IL1B', 'Macro_INHBA', 'Macro_SPP1', 'Macro_FN1', 'Macro-ISG', 
             'Macro_TNF', 'Macro_LYVE1', 'Macro_C1QC', 'Macro_TREM2')
mural <- c("Pericytes","SMC")
caf <- c("CAF_SFRP2", "CAF-prog", "CAF-desmo", "iCAF_MMP1", "iCAF_IL6", "CAF-ap","Myofibroblasts")
endo <- c("Endo-lymphatic", "Endo-artery", "Endo-capillary", "Endo-tip", "Endo-vein")
plasma <- c("PC-cycling","PC-trans","PC-early_RGS13", "PC_IGHG", "PC_IGHA")
b <- c("B-naive", "B-ISG", "B-HSP", "B_MT2A", 
       "ACB_EGR1", "ACB_NR4A2", "ACB_CCR7", "B-memory", "B-AtM", 
       "GCB-pre", "GCB-DZ_SUGCT", "GCB-LZ_LMO2",
       "GCB-cycling")
putative_malignant <- c('Epithelial(CNA-)','Epithelial(CNA+)','Melanocytes(CNA-)','Melanocytes(CNA+)','Malignant(CNA+)')



