library(Matrix)
library(data.table)
setwd("~/Desktop/STATS254/opMultiome10x")

# Binerize and compress MTX data
matFile <- "opmultiome_GEXATAC_s2.quantitative.mtx"
mat <- Matrix::readMM(matFile)
#mat <- as.matrix((mat > 0) + 0)
mat <- as(mat, "nsparseMatrix")
head(mat)
Matrix::writeMM(mat, file=paste0('~/Desktop/STATS254/opMultiome10x/opmultiome_GEXATAC_s2.binarized','.mtx.gz'))
rm(mat)
