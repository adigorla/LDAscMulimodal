library(devtools)
devtools::install_github("aertslab/cisTopic", args = c('--library="~/Desktop/STATS254/repo"'))
library(cisTopic)

#cisTopicObject <- readRDS('/Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/library/cisTopic/examples/cisTopicObject_melanoma.Rds')

#epimat <- Matrix::readMM('~/Desktop/STATS254/mouseData/activity_scores.quantitative.mtx.gz')
epimat <- as.matrix(readRDS('~/Desktop/STATS254/mouseData/activity_scores.binarized.rds'))
epimat <- epimat[, 1:18000]
gene.names = read.delim('~/Desktop/STATS254/mouseData/activity_scores.binarized.genes.txt',
                           header = FALSE,
                           stringsAsFactors = FALSE)
cell.names = read.delim('~/Desktop/STATS254/mouseData/activity_scores.binarized.cells.txt',
                        header = FALSE,
                        stringsAsFactors = FALSE)
rownames(epimat) <- gene.names$V1
colnames(epimat) <- cell.names$V1
head(epimat)
cisTopicObject <- createcisTopicObjectFrom10Xmatrix(epimat, project.name='mouseATLAS_epi')
