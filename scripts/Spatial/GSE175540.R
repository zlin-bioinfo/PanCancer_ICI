# #!/bin/bash
# 
# # Define the directory path
# dir_path="data/GSE175540/GSE175540_RAW"
# 
# # Iterate through all files in the directory
# for file in "$dir_path"/*; do
# # Extract the base filename (without path)
# base_name=$(basename "$file")
# 
# # Remove the GSM prefix and extract the sample ID (e.g., c_2, a_1)
# # Assumes the format is GSM[0-9]+_(ffpe|frozen)_sampleid_*
# sample_id=$(echo "$base_name" | sed -E 's/^GSM[0-9]+_(ffpe|frozen)_//' | cut -d'_' -f1-2)
# 
# # Create a new directory for the sample if it doesn't exist
# new_dir="$dir_path/$sample_id"
# mkdir -p "$new_dir"
# 
# # Move the file to the new directory
# mv "$file" "$new_dir/$base_name"
# 
# # Unzip the file if it ends with .gz
# if [[ "$base_name" == *.gz ]]; then
# gunzip "$new_dir/$base_name"
# fi
# done

# base_dir="/bigdata/zlin/PanCancer_ICI/data/GSE175540/GSE175540_RAW"
# 
# for subdir in "$base_dir"/*; do
# if [ -d "$subdir" ]; then
# mkdir -p "$subdir/spatial"
# mv "$subdir/aligned_fiducials.jpg" \
# "$subdir/detected_tissue_image.jpg" \
# "$subdir/scalefactors_json.json" \
# "$subdir/tissue_hires_image.png" \
# "$subdir/tissue_lowres_image.png" \
# "$subdir/tissue_positions_list.csv" \
# "$subdir/spatial/" 2>/dev/null
# fi
# done

library(Seurat)
library(ggplot2)
library(patchwork)
library(janitor)
library(dplyr)
imaged <- Read10X_Image(image.dir = 'data/GSE175540/GSE175540_RAW/a_1/spatial',
                        image.name = 'tissue_lowres_image.png',
                        slice = 'a_1')
seu1 <- Load10X_Spatial('data/GSE175540/GSE175540_RAW/a_1/',
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial",
                       slice = "a_1", image = imaged)
seu1$sample <- 'a_1'
imaged <- Read10X_Image(image.dir = 'data/GSE175540/GSE175540_RAW/a_15/',
                        image.name = 'tissue_lowres_image.png',
                        slice = 'a_15')
seu2 <- Load10X_Spatial('data/GSE175540/GSE175540_RAW/a_15/',
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial",
                       slice = "a_15", image = imaged)
seu2$sample <- 'a_15'
seu <- merge(seu1, seu2)
Idents(seu) <- seu$sample
VlnPlot(seu, features = c('nCount_Spatial', 'nFeature_Spatial'))
seu <- seu |> NormalizeData()
seu <- JoinLayers(seu)
mye_sig <- read.csv('tables/marker_myeloids.csv') 
macc1qc <- mye_sig$names.Macro_C1QC
seu <- AddModuleScore(seu, list(macc1qc), name = 'Macro_C1QC')
cd4t_sig <- read.csv('tables/marker_cd4t.csv') 
cd4t <- cd4t_sig$names.CD4_Treg
seu <- AddModuleScore(seu, list(cd4t), name = 'CD4_Treg')
SpatialFeaturePlot(seu, features = c('CD4_Treg1','Macro_C1QC1','CD4_T-naive1'), max.cutoff = 1, images = 'a_15', image.alpha = 0.5, min.cutoff = c(0,0.1,0.1))







