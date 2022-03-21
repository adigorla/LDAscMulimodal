library(Matrix)
library(data.table)
setwd("~/Desktop/STATS254/scIB")

# Binerize and compress MTX data
matFile <- "scBI_mousebrain_GEXATAC.quantitative.mtx"
mat <- Matrix::readMM(matFile)
#mat <- as.matrix((mat > 0) + 0)
mat <- as(mat, "nsparseMatrix")
head(mat)
Matrix::writeMM(mat, file=paste0('./scBI_mousebrain_GEXATAC.binarized','.mtx.gz'))
rm(mat)
