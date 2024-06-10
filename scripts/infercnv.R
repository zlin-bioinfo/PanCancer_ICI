#!/usr/bin/env Rscript
rm(list=ls())
pkgs <- c('Seurat','infercnv','magrittr','qs','dittoSeq','ComplexHeatmap','RColorBrewer','colorRamp2')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
options(scipen = 100)

# seu_list <- readRDS('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/processed.rds') %>% SplitObject(split.by = "patient")
# # https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt
# gene_order <- read.table('/bigdata/zlin/Melanoma_meta/data/hg38_gencode_v27.txt', header = F,row.names = 1)
# for (i in 1:length(seu_list)){
#   seu <- seu_list[[i]]
#   patient <- unique(seu$patient)
#   print(paste0('Patient ', i))
#   # counts matrix
#   counts_matrix <- as.data.frame(seu@assays$RNA@counts)
#   # classify immune/non-immune cells
#   anno <- data.frame(celltype = seu$celltype_major)
#   anno$cell <- 'Others'
#   anno$cell[anno$celltype == 'Melanoma' ] <- 'Melanoma'
#   anno$celltype <- NULL
#   # Create infercnv object
#   infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
#                                        annotations_file = anno,
#                                        gene_order_file = gene_order,
#                                        ref_group_names = 'Others')
#   # Run infercnv
#   path_dir <- paste0('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/inferCNV/inferCNV_subcluster_', patient)
#   infercnv_obj <- suppressWarnings(infercnv::run(infercnv_obj,
#                                                  cutoff = 0.1,
#                                                  out_dir = path_dir,
#                                                  cluster_by_groups = F,
#                                                  denoise = T,
#                                                  analysis_mode = "subclusters",
#                                                  HMM = T,
#                                                  output_format = 'pdf',
#                                                  num_threads = 30))
#   # Identify malignant cells(cut-off >0.1)
#   seu <- add_to_seurat(seu, infercnv_output_path = path_dir)
#   cnv_cols <- grep('proportion_scaled_cnv_chr', names(seu@meta.data), value = T)
#   cnvs <- seu@meta.data[, cnv_cols]
#   seu$cnv_avg <- rowMeans(cnvs)
#   seu$malignant <- ifelse(seu$cnv_avg > 0.1, 'malignant', 'non-malignant')
#   # Add CNV metrics
#   cnv_cols <- grep('proportion_scaled_cnv_chr', names(seu@meta.data), value = T)
#   cnvs <- seu@meta.data[, cnv_cols]
#   seu$proportion_scaled_cnv_avg <- rowMeans(cnvs)
#   
#   cnv_cols <- grep('proportion_cnv_chr', names(seu@meta.data), value = T)
#   cnvs <- seu@meta.data[, cnv_cols]
#   seu$proportion_cnv_avg <- rowMeans(cnvs)
#   
#   cnv_cols <- grep('has_cnv_chr', names(seu@meta.data), value = T)
#   cnvs <- seu@meta.data[, cnv_cols]
#   seu$has_cnv_avg <- rowMeans(cnvs)
#   seu_list[[i]] <- seu
# }
# seu <- merge(x = seu_list[[1]], y = seu_list[(2:length(seu_list))])
# saveRDS(seu, file = '/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/processed.rds')
# 
# # SKCM_Becker  HNSC_Franken BRCA_Bassez1 BRCA_Bassez2 BCC_Yost CRC_Li TNBC_Li
# dataset <- 'SKCM_Becker'
# seu <- dior::read_h5(file = paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/anno_r2.h5'), target.object = 'seurat')
# # seu <- dior::read_h5(file = paste0('/bigdata/zlin/Melanoma_meta/data/', dataset, '/seu.h5'), target.object = 'seurat')
# seu <- CreateSeuratObject(counts = seu@assays$RNA@counts, meta.data = seu@meta.data)
# Idents(seu) <- seu$celltype_major
# seu <- subset(seu, subset = celltype_major %in% c('Melanoma','Endothelium','Fibroblasts','Myeloids'))
# seu$anno <- as.character(seu$celltype_major)
# seu$patient <- as.character(seu$patient)
# seu$anno[seu$celltype_major =='Melanoma'] <- seu$patient[seu$celltype_major =='Melanoma']
# Idents(seu) <- seu$anno

output_dir_full = '/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/infercnv'
# https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt
gene_order <- read.table('/bigdata/zlin/Melanoma_meta/data/hg38_gencode_v27.txt', header = F,row.names = 1)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=GetAssayData(seu,assay = "RNA", slot ='counts'),
                                         annotations_file=as.matrix(seu@active.ident),
                                         delim="\t",
                                         gene_order_file=gene_order,
                                         ref_group_names=c("Endothelium" ,"Fibroblasts", "Myeloids"))

infercnv_obj = suppressWarnings(infercnv::run(infercnv_obj,
                                              cutoff=0.1, 
                                              out_dir=output_dir_full, 
                                              cluster_by_groups=T,
                                              cluster_references = F,
                                              leiden_resolution = 0.001,
                                              HMM=T,
                                              per_chr_hmm_subclusters=F,
                                              denoise=T,
                                              num_threads = 60,
                                              no_plot=T))

infercnv_obj <- readRDS('/bigdata/zlin/Melanoma_meta/data/SKCM_Becker/infercnv/preliminary.infercnv_obj')
expr <- infercnv_obj@expr.data[, unlist(infercnv_obj@observation_grouped_cell_indices)] |> t()
if(T){
  tmp1 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-1`]
  
  dim(tmp1)
  head(tmp1)[,1:9]
  
  tmp2 = expr[,infer_CNV_obj@reference_grouped_cell_indices$`ref-2`]
  tmp= cbind(tmp1,tmp2)
  
  dim(tmp)
  head(tmp)[,1:9]
  
  down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))
  up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))
  oneCopy=up-down
  oneCopy
  
  a1= down- 2*oneCopy
  a2= down- 1*oneCopy
  down;up
  a3= up +  1*oneCopy
  a4= up + 2*oneCopy 
  
  cnv_score_table<-infer_CNV_obj@expr.data
  cnv_score_table[1:4,1:4]
  dim(cnv_score_table)
  
  cnv_score_mat <- as.matrix(cnv_score_table)
  # Scoring
  cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A" #complete loss. 2pts
  cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B" #loss of one copy. 1pts
  cnv_score_table[cnv_score_mat >= down & cnv_score_mat <  up ] <- "C" #Neutral. 0pts
  cnv_score_table[cnv_score_mat >= up  & cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
  cnv_score_table[cnv_score_mat > a3  & cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
  cnv_score_table[cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts
  
  # Check
  table(cnv_score_table[,1])
  # Replace with score 
  cnv_score_table_pts <- cnv_score_table
  rm(cnv_score_mat)
  # 
  cnv_score_table_pts[cnv_score_table == "A"] <- 2
  cnv_score_table_pts[cnv_score_table == "B"] <- 1
  cnv_score_table_pts[cnv_score_table == "C"] <- 0
  cnv_score_table_pts[cnv_score_table == "D"] <- 1
  cnv_score_table_pts[cnv_score_table == "E"] <- 2
  cnv_score_table_pts[cnv_score_table == "F"] <- 2
  
  # Scores are stored in "cnv_score_table_pts". Use colSums to add up scores for each cell and store as vector 
  cnv_score_table_pts[1:4,1:4]
  tmp  = apply(cnv_score_table_pts, 1, as.numeric)
  tmp[1:4,1:4]
  rownames(tmp)  = colnames(cnv_score_table_pts) 
  cell_scores_CNV <- as.data.frame(rowSums(tmp))
  colnames(cell_scores_CNV) <- "cnv_score"
  head(cell_scores_CNV)
  #write.csv(x = cell_scores_CNV, file = "cnv_scores.csv")
  
}







# reproducing cnv heatmap (not working because the row limit of Complexheatmap)
gene_order <- infercnv_obj@gene_order
color_mapping <- colorRamp2(c(0.8, 1,1.2), c("#2166AC","white","#B2182B"))
patient <- seu@meta.data[unlist(infercnv_obj@observation_grouped_cell_indices), 'patient']
timepoint <- seu@meta.data[unlist(infercnv_obj@observation_grouped_cell_indices), 'time_point']
response <- seu@meta.data[unlist(infercnv_obj@observation_grouped_cell_indices), 'response']
response <- ifelse(response %in% c('CR','PR'), 'RE', 
                   ifelse(response %in% c('SD','PD'), 'NR', 'NE'))
t_vec <- brewer.pal(9, "Set1")[1:length(unique(timepoint))]
names(t_vec) <- unique(timepoint)
r_vec <- brewer.pal(9, "Set1")[3:(length(unique(timepoint))-2)]
names(r_vec) <- unique(response)
p_vec <- dittoColors()[1:length(unique(patient))]
names(p_vec) <- unique(patient)
row_ha = rowAnnotation(Time_point = timepoint, Response = response, Patient = patient,col = list(Time_point = t_vec, Response = r_vec, Patient = p_vec))
pdf("/bigdata/zlin/Melanoma_meta/figures/heatmap_Becker_cnv.pdf", width = 8, height = 8)
ht <- Heatmap(expr, 
        name = 'InferCNV',
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        column_split = gene_order$chr,
        column_title_side = c("bottom"),
        col = color_mapping,
        left_annotation = row_ha,
        column_title_gp = grid::gpar(fontsize = 6),
        use_raster = TRUE
        )
draw(ht)
dev.off()




