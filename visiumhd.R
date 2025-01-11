library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

object <- Load10X_Spatial(data.dir = "/bigdata/zlin/data/1325_1_XS-VHD/",
                          bin.size = c(8, 16))

# Setting default assay changes between 8um and 16um binning
Assays(object)
DefaultAssay(object) <- "Spatial.008um"

vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot



function (data.dir, filename = "filtered_feature_bc_matrix.h5", 
          assay = "Spatial", slice = "slice1", bin.size = NULL, filter.matrix = TRUE, 
          to.upper = FALSE, image = NULL, ...) 
{
  if (length(x = data.dir) > 1) {
    data.dir <- data.dir[1]
    warning(paste0("`data.dir` expects a single value but recieved multiple - ", 
                   "continuing using the first: '", data.dir, "'."), 
            immediate. = TRUE, )
  }
  if (!file.exists(data.dir)) {
    stop(paste0("No such file or directory: ", "'", data.dir, 
                "'"))
  }
  if (is.null(bin.size) & file.exists(paste0(data.dir, "/binned_outputs"))) {
    bin.size <- c(16, 8)
  }
  if (!is.null(bin.size)) {
    bin.size.pretty <- paste0(sprintf("%03d", bin.size), 
                              "um")
    data.dirs <- paste0(data.dir, "/binned_outputs/", "square_", 
                        bin.size.pretty)
    assay.names <- paste0(assay, ".", bin.size.pretty)
    slice.names <- paste0(slice, ".", bin.size.pretty)
  }
  else {
    data.dirs <- data.dir
    assay.names <- assay
    slice.names <- slice
  }
  counts.paths <- lapply(data.dirs, file.path, filename)
  counts.list <- lapply(counts.paths, Read10X_h5, ...)
  if (to.upper) {
    rownames(counts) <- lapply(rownames(counts), toupper)
  }
  if (is.null(image)) {
    image.list <- mapply(Read10X_Image, file.path(data.dirs, 
                                                  "spatial"), assay = assay.names, slice = slice.names, 
                         MoreArgs = list(filter.matrix = filter.matrix))
  }
  else {
    image.list <- c(image)
  }
  if (length(image.list) != length(counts.list)) {
    stop(paste0("The number of images does not match the number of counts matrices. ", 
                "Ensure each spatial dataset has a corresponding image."))
  }
  object.list <- mapply(CreateSeuratObject, counts.list, assay = assay.names)
  object.list <- mapply(function(.object, .image, .assay, .slice) {
    .image <- .image[Cells(.object)]
    .object[[.slice]] <- .image
    return(.object)
  }, object.list, image.list, assay.names, slice.names)
  object <- merge(object.list[[1]], y = object.list[-1])
  return(object)
}