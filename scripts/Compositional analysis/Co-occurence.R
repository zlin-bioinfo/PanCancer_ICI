rm(list=ls())
pkgs <- c('tidyr','plyr','dplyr','stringr','ggsci','patchwork','ggplot2','RColorBrewer','tibble','qs2','MetBrewer','grid','corrplot','ComplexHeatmap','colorRamp2','corrr','igraph','ggraph','tidygraph','graphlayouts','ggforce','ggnetwork','ggnewscale','janitor')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)
metadata <- read.csv('tables/meta_all.csv') 
source('scripts/Celltype_classification.R')
# unwanted_celltypes <- c('Melanocytes(CNA-)', 'Melanocytes(CNA+)', 'Epithelial(CNA+)', 'Epithelial(CNA-)', 'Malignant(CNA+)')
unwanted_celltypes <- c('Melanocytes(CNA-)', 'Epithelial(CNA-)','Cycling T','Cycling NK', "GCB-cycling", "PC-cycling", 'Cycling non-immune', 'Cycling myeloids')
metadata$freq_r2_comp[metadata$celltype_r2 == 'Malignant(CNA+)'] <- metadata$freq_r2[metadata$celltype_r2 == 'Malignant(CNA+)']
# Overall
mat <- metadata |> 
  filter(!celltype_r2 %in% unwanted_celltypes, time_point == c('ICI_exposed','On'), subset == 'All TME') |> 
  select(sample, celltype_r2, freq_r2_comp) |> 
  distinct(sample, celltype_r2, freq_r2_comp, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = sample, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
cor.res <- cor(t(mat))
row.order <- hclust(dist(cor.res))$order
cor.res <- cor.res[row.order, row.order]
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cor.test.res <- testRes$p
cor.test.res <- cor.test.res[row.order, row.order]
diag(cor.test.res) <- NA
# pdf('figures/Co-occurence/pre.pdf', height = 12, width = 12)
ht1 <- Heatmap(cor.res, col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), 
        cluster_rows = F, cluster_columns = F, column_names_rot = 45, 
        show_row_names = T, show_column_names = T,
        row_names_side = "right", column_names_side = 'top',
        show_column_dend = F, show_row_dend = F,
        heatmap_legend_param = list(title = "Pearson's \nCorrelation", 
                                    at = c(-0.8, 0, 0.8), 
                                    labels = c("0.8", "0", "0.8"),
                                    legend_direction = "horizontal", 
                                    legend_width = unit(1.5, "cm"), 
                                    legend_side = 'bottom',
                                    title_position = "topcenter"),
        rect_gp = gpar(type = "none",col ="grey", lwd = 0.5), 
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 8),
        width = ncol(cor.res)*unit(2.5, "mm"),
        height = nrow(cor.res)*unit(2.5, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i <= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
          if(!is.na(cor.test.res[i, j])) {
            if(cor.test.res[i, j] < 0.05 & i < j) {
              gb = textGrob("*", gp = gpar(fontsize=14))
              gb_h = convertHeight(grobHeight(gb), "mm")
              gb_w = convertWidth(grobWidth(gb), "mm")
              grid.text('*', x, y - gb_h*0.5 + gb_w*0.4, gp=gpar(fontsize=14))
            }
          }
        })
# dev.off()

row.name <- rownames(cor.res)
mat <- metadata |> 
  filter(!celltype_r2 %in% unwanted_celltypes, time_point == c('No_ICI','Pre'), subset == 'All TME') |> 
  select(sample, celltype_r2, freq_r2_comp) |> 
  distinct(sample, celltype_r2, freq_r2_comp, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = sample, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
# meta_matrix <- apply(mat, 2, function(col) scale(log(col + 1e-6)))
# rownames(meta_matrix) <- rownames(mat)
cor.res <- cor(t(mat))
cor.res <- cor.res[row.name, row.name]
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cor.test.res <- testRes$p
cor.test.res <- cor.test.res[row.name, row.name]
diag(cor.test.res) <- NA
# pdf('figures/Co-occurence/on.pdf', height = 12, width = 12)
ht2 <- Heatmap(cor.res, col = circlize::colorRamp2(c(-0.8, 0, 0.8), c("#154999", "white", "#CF0034")), 
        cluster_rows = F, cluster_columns = F, column_names_rot = 90, row_names_side = "left", column_names_side = 'bottom',
        show_column_dend = F, show_row_dend = F,
        heatmap_legend_param = list(title = "Pearson's \nCorrelation", 
                                    at = c(-0.8, 0, 0.8), 
                                    labels = c("0.8", "0", "0.8"),
                                    legend_direction = "horizontal", 
                                    legend_width = unit(1.5, "cm"), 
                                    legend_side = 'bottom',
                                    title_position = "topcenter"),
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        width = ncol(cor.res)*unit(2.5, "mm"),
        height = nrow(cor.res)*unit(2.5, "mm"),
        rect_gp = gpar(type = "none",col ="grey", lwd = 0.5), 
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
          if(!is.na(cor.test.res[i, j])) {
            if(cor.test.res[i, j] < 0.05 & i > j) {
              gb = textGrob("*", gp = gpar(fontsize=14))
              gb_h = convertHeight(grobHeight(gb), "mm")
              gb_w = convertWidth(grobWidth(gb), "mm")
              grid.text('*', x, y - gb_h*0.5 + gb_w*0.4, gp=gpar(fontsize=14))
            }
          }
        })
# dev.off()
pdf('figures/Co-occurence/pre_on.pdf', height = 12, width = 12)
draw(ht2 + ht1, ht_gap = unit(-195, "mm"))
dev.off()

mat <- metadata |> 
  filter(!celltype_r2 %in% unwanted_celltypes, time_point == c('ICI_exposed','On'), subset == 'All TME') |> 
  select(sample, celltype_r2, freq_r2_comp) |> 
  distinct(sample, celltype_r2, freq_r2_comp, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = sample, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
# meta_matrix <- apply(mat, 2, function(col) scale(log(col + 1e-6)))
# rownames(meta_matrix) <- rownames(mat)
cor.res <- cor(t(mat))
row.order <- hclust(dist(cor.res))$order
cor.res <- cor.res[row.order, row.order]
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cor.test.res <- testRes$p
cor.test.res <- cor.test.res[row.order, row.order]


row.name <- rownames(cor.res)
mat <- metadata |> 
  filter(!celltype_r2 %in% unwanted_celltypes, time_point == c('No_ICI','Pre'), subset == 'All TME') |> 
  select(sample, celltype_r2, freq_r2_comp) |> 
  distinct(sample, celltype_r2, freq_r2_comp, .keep_all = T) |> 
  pivot_wider(values_from = freq_r2_comp, names_from = sample, values_fill = 0) |>
  column_to_rownames(var = 'celltype_r2')
# meta_matrix <- apply(mat, 2, function(col) scale(log(col + 1e-6)))
# rownames(meta_matrix) <- rownames(mat)
cor.res1 <- cor(t(mat))
cor.res1 <- cor.res1[row.name, row.name]
testRes <- cor.mtest(t(mat), conf.level = 0.95)
cor.test.res <- testRes$p
cor.test.res <- cor.test.res[row.name, row.name]
pdf('figures/Co-occurence/ht_delta.pdf', height = 12, width = 12)
Heatmap((cor.res1-cor.res), col = circlize::colorRamp2(c(-0.5, 0, 0.5), c("#154999", "white", "#CF0034")), 
         column_names_rot = 90, row_names_side = "left", column_names_side = 'bottom',
        show_column_dend = F, show_row_dend = F,
        heatmap_legend_param = list(title = "Pearson's \nCorrelation(delta)", 
                                    at = c(-0.5, 0, 0.5),
                                    legend_width = unit(1.5, "cm")),
        column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8),
        rect_gp = gpar(col ="grey", lwd = 0.5), 
        width = ncol(cor.res)*unit(2.5, "mm"),
        height = nrow(cor.res)*unit(2.5, "mm")
        )
dev.off()
Heatmap(cor.res)






