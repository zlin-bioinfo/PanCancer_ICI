rm(list=ls())
pkgs <- c('GEOquery')
unlist(lapply(pkgs, function(x) require(package = x,  character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
options(warn = -1)

# GSE115821 Melanoma
GEO <- "GSE115821"
gds <- getGEO("GSE115821", destdir = "/bigdata/zlin/Melanoma_meta/data/bulk_datasets/Melanoma/GSE115821")
meta <- rbind(pData(gds[[1]]), pData(gds[[2]]))
# getGEOSuppFiles(GEO, baseDir = "/bigdata/zlin/Melanoma_meta/data/bulk_datasets/Melanoma/GSE115821",
#                 fetch_files = TRUE, filter_regex = NULL)








