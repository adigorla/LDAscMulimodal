library(Matrix)
library(data.table)
library(ggplot2)
library(tidyverse)
setwd("~/Desktop/STATS254/mouseData")


# Find Shared Genes
genes.atac = read.delim('~/Desktop/STATS254/mouseData/activity_scores.binarized.genes.txt',
                        header = FALSE,
                        stringsAsFactors = FALSE)
print(dim(genes.atac))

genes.rna = read.delim('~/Desktop/STATS254/mouseData/droplet/Marrow-10X_P7_3/genes.tsv',
                        header = FALSE,
                        sep = '\t',
                        stringsAsFactors = FALSE)
print(dim(genes.rna))

genes.shared <- intersect(
  as.list(genes.atac[,1]),
  lapply(as.list(genes.rna[,2]), toupper)
)
genes.shared <- as.vector(genes.shared)
print(length(genes.shared))


# Process RNA annot file
annot.rna = read.delim('~/Desktop/STATS254/mouseData/annotations_droplet.csv',
                       header = TRUE,
                       sep = ',')

as.data.frame(table(annot.rna$cell_ontology_class)) %>% 
ggplot(aes(x = Var1, y = Freq)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_bar(stat="identity")
as.data.frame(table(annot.rna$tissue)) %>% 
  ggplot(aes(x = Var1, y = Freq)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_bar(stat="identity")

print(attributes(table(annot.rna$cell_ontology_class))$dimnames[[1]])


# Process ATAC annot file
annot.atac = read.delim('~/Desktop/STATS254/mouseData/cell_metadata.tissue_freq_filtered.txt',
                       header = TRUE,
                       sep = '\t')

as.data.frame(table(annot.atac$cell_label)) %>% 
  ggplot(aes(x = Var1, y = Freq)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_bar(stat="identity")
as.data.frame(table(annot.atac$tissue)) %>% 
  ggplot(aes(x = Var1, y = Freq)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_bar(stat="identity")

print(attributes(table(annot.atac$cell_label))$dimnames[[1]])


# Renaming cell types for greater concordance 
annot.rna$cell_ontology_class[annot.rna$cell_ontology_class =='cardiac muscle cell'] <- 'Cardiomyocytes'
annot.rna$cell_ontology_class[(annot.rna$cell_ontology_class == 'early pro-B cell') | (annot.rna$cell_ontology_class == 'Fraction A pre-pro B cell') |
             (annot.rna$cell_ontology_class == 'immature B cell') | (annot.rna$cell_ontology_class == 'late pro-B cell')] <- 'B cell'
annot.rna$cell_ontology_class[(annot.rna$cell_ontology_class == 'immature T cell')] <- 'T cell'
annot.rna$cell_ontology_class[(annot.rna$cell_ontology_class == 'alveolar macrophage')] <- 'macrophage'
annot.rna$cell_ontology_class[(annot.rna$cell_ontology_class == 'natural killer cell')] <- 'NK cell'
annot.rna$cell_ontology_class[(annot.rna$cell_ontology_class == 'endothelial cell of hepatic sinusoid')] <- 'endothelial cell'
annot.rna$cell_ontology_class[(annot.rna$cell_ontology_class == 'hematopoietic precursor cell')] <- 'hematopoietic progenitor'
annot.rna$cell_ontology_class[(annot.rna$cell_ontology_class == 'non-classical monocyte')] <- 'monocyte'
annot.rna$cell_ontology_class[(annot.rna$cell_ontology_class == 'blood cell')] <- 'Enterocytes'

annot.atac$cell_label[(annot.atac$cell_label == 'Regulatory T cell') | (annot.atac$cell_label == 'T cells')] <- 'T cell'
annot.atac$cell_label[(annot.atac$cell_label == 'Immature B cells') | (annot.atac$cell_label == 'Activated B cells') | (annot.atac$cell_label == 'B cells')] <- 'B cell'
annot.atac$cell_label[(annot.atac$cell_label == 'Ex. neurons CPN') | (annot.atac$cell_label == 'Ex. neurons CThPN') |
             (annot.atac$cell_label == 'Ex. neurons SCPN') | (annot.atac$cell_label == 'Inhibitory neurons') ] <- 'Neuron'
annot.atac$cell_label[(annot.atac$cell_label == 'Alveolar macrophages') | (annot.atac$cell_label == 'Macrophages')] <- 'macrophage'
annot.atac$cell_label[(annot.atac$cell_label == 'NK cells')] <- 'NK cell'
annot.atac$cell_label[(annot.atac$cell_label == 'Dendritic cells')] <- 'Dendritic cell'
annot.atac$cell_label[(annot.atac$cell_label == 'Endothelial I (glomerular)') | (annot.atac$cell_label == 'Endothelial II cells') |
             (annot.atac$cell_label == 'Endothelial I cells') ] <- 'endothelial cell'
annot.atac$cell_label[(annot.atac$cell_label == 'Erythroblasts')] <- 'erythroblast'
annot.atac$cell_label[(annot.atac$cell_label == 'Hematopoietic progenitors')] <- 'hematopoietic progenitor'
annot.atac$cell_label[(annot.atac$cell_label == 'Hepatocytes')] <- 'hepatocyte'
annot.atac$cell_label[(annot.atac$cell_label == 'Monocytes')] <- 'monocyte'
annot.atac$cell_label[(annot.atac$cell_label == 'Proximal tubule S3')] <- 'Proximal tubule'
annot.atac$cell_label[ (annot.atac$cell_label == 'Type II pneumocytes')] <- 'type II pneumocyte'


annot.rna$tissue[annot.rna$tissue == 'Heart_and_Aorta'] <- 'Heart'
annot.rna$tissue[annot.rna$tissue == 'Marrow'] <- 'BoneMarrow'



as.data.frame(table(annot.rna$cell_ontology_class)) %>% 
  ggplot(aes(x = Var1, y = Freq)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_bar(stat="identity")
as.data.frame(table(annot.rna$tissue)) %>% 
  ggplot(aes(x = Var1, y = Freq)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_bar(stat="identity")
as.data.frame(table(annot.atac$cell_label)) %>% 
  ggplot(aes(x = Var1, y = Freq)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_bar(stat="identity")
as.data.frame(table(annot.atac$tissue)) %>% 
  ggplot(aes(x = Var1, y = Freq)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_bar(stat="identity")

print(attributes(table(annot.rna$cell_ontology_class))$dimnames[[1]])
print(attributes(table(annot.atac$cell_label))$dimnames[[1]])
shared.cellTypes<- intersect(attributes(table(annot.rna$cell_ontology_class))$dimnames[[1]], 
                 attributes(table(annot.atac$cell_label))$dimnames[[1]])
print(length(shared.cellTypes))
shared.tissues<- intersect(attributes(table(annot.rna$tissue))$dimnames[[1]], 
                             attributes(table(annot.atac$tissue))$dimnames[[1]])
print(length(shared.tissues))


# Compile a single matrix for RNA, and a simple annot file 
# annot file has: cellID, cell_label, tissue, source
scrnaFiles <- list.files("~/Desktop/STATS254/mouseData/droplet", pattern="matrix.mtx", recursive=TRUE, full.names=TRUE)

NEWannot.rna <- data.frame(cellID = character(),
                           cellLabel = character(),
                           tissue = character(),
                           source = character())
cellIDs <- list()
batchIDs <- list()
tissIDs <- c()
cLabs <- c()
rnaMAT <- NULL
ct = 1
for (matFile in scrnaFiles){
  
  subd<- strsplit(matFile,"/")[[1]]
  batch <- subd[length(subd) - 1]
  batchid <- strsplit(batch, '-')[[1]]
  
  bcFile <- read.delim(gsub("matrix.mtx", "barcodes.tsv", matFile),
                                   header = FALSE,
                                   sep = '\t',
                                   stringsAsFactors = FALSE)
  idx = 1
  print(dim(bcFile))
  dropcell <- c()
  for (cid in bcFile[,1]) {
    tmp<- strsplit(cid, '-')[[1]][1]
    bcFile[idx,1] <- paste(batchid[2], tmp, sep='_')
    row <- annot.rna[annot.rna$cell == bcFile[idx,1], ]
    tissIDs <- c(tissIDs, row$tissue)
    cLabs <- c(cLabs, row$cell_ontology_class)
    if (length(row$tissue) != 1 | length(row$cell_ontology_class) != 1){
      dropcell <- c(dropcell, bcFile[idx,1])
    }
    idx = idx + 1
  }
  print(batchid)
  batchIDs[[ct]] <- batchid[2]
  gFile <- read.delim(gsub("matrix.mtx", "genes.tsv", matFile),
                      header = FALSE,
                      sep = '\t',
                      stringsAsFactors = FALSE)
  mat <- as.matrix(Matrix::readMM(matFile))
  rownames(mat) <- lapply(genes.rna[,2], toupper)
  colnames(mat) <- bcFile[,1]
  print(length(colnames(mat)))
  mat <- mat[rownames(mat) %in% genes.shared,]
  if (length(dropcell) > 0 ){
    mat <- mat[, !colnames(mat) %in% dropcell]
  }
  cellIDs[[ct]] <- colnames(mat)
  mat <- as.matrix((mat > 0) + 0)
  mat <- Matrix::t(as(mat, "dgCMatrix"))
  print("mat binerized")
  if (!is.null(rnaMAT)){
    rnaMAT <- rbind(rnaMAT, mat)
  } else{
    rnaMAT <- mat
    print("firstmat")
  }
  #Matrix::writeMM(mat, file=paste0('~/Desktop/STATS254/mouseData/droplet/',batchid[2],'.mtx.gz'))
  #rnaMAT[[ct]] <- mat
  ct = ct + 1
}
print(dim(rnaMAT))
Matrix::writeMM(rnaMAT, file=paste0('~/Desktop/STATS254/mouseData/scrna10X.binarized','.mtx.gz'))
NEWannot.rna <- data.frame(cellID = unlist(cellIDs))
NEWannot.rna$cellLabel <- cLabs
NEWannot.rna$tissue <- tissIDs
NEWannot.rna$source <- rep("scRNA", length(NEWannot.rna$cellID))
write.table(NEWannot.rna, file=paste0('~/Desktop/STATS254/mouseData/scrna10X.binarized.annot','.tsv'), quote=FALSE, sep='\t', row.names = FALSE)
write.table(as.data.frame(colnames(rnaMAT)), file=paste0('~/Desktop/STATS254/mouseData/scrna10X.binarized.genes','.tsv'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
rm(rnaMAT, mat, cLabs, tissIDs, cellIDs, ct, idx, batch, batchid)


# Compile a single matrix for ATAC, and a simple annot file 
# annot file has: cellID, cell_label, tissue, source

atacMAT <- Matrix::t(Matrix::readMM('~/Desktop/STATS254/mouseData/activity_scores.binarized.mtx.gz'))*1
gene.names = read.delim('~/Desktop/STATS254/mouseData/activity_scores.binarized.genes.txt',
                        header = FALSE,
                        stringsAsFactors = FALSE)
cell.names = read.delim('~/Desktop/STATS254/mouseData/activity_scores.binarized.cells.txt',
                        header = FALSE,
                        stringsAsFactors = FALSE)
colnames(atacMAT) <- gene.names$V1
rownames(atacMAT) <- cell.names$V1
atacMAT <- atacMAT[, colnames(atacMAT) %in% genes.shared]

# Drop Testes, Cerebellum and unkown cells
dropcell <- annot.atac$cell[annot.atac$tissue == "Testes"]
dropcell <- c(dropcell, annot.atac$cell[annot.atac$cell_label == "Unknown"])
dropcell <- c(dropcell, annot.atac$cell[annot.atac$tissue == "Cerebellum"])
length(dropcell)
atacMAT <- atacMAT[!rownames(atacMAT) %in% dropcell,]

cellIDs <- rownames(atacMAT)
tissIDs <- c()
cLabs <- c()
dropcell <- c()
for (cid in cellIDs) {
  row <- annot.atac[annot.atac$cell == cid, ]
  tissIDs <- c(tissIDs, row$tissue)
  cLabs <- c(cLabs, row$cell_label)
  if (length(row$tissue) != 1 | length(row$cell_label) != 1){
    dropcell <- c(dropcell, cid)
  }
}
atacMAT <- atacMAT[!rownames(atacMAT) %in% dropcell,]
cellIDs <- rownames(atacMAT)
dim(atacMAT)

Matrix::writeMM(atacMAT, file=paste0('~/Desktop/STATS254/mouseData/scatacAS.binarized','.mtx.gz'))
NEWannot.atac <- data.frame(cellID = cellIDs)
NEWannot.atac$cellLabel <- cLabs
NEWannot.atac$tissue <- tissIDs
NEWannot.atac$source <- rep("scATAC", length(NEWannot.atac$cellID))
write.table(NEWannot.atac, file=paste0('~/Desktop/STATS254/mouseData/scatacAS.binarized.annot','.tsv'), quote=FALSE, sep='\t', row.names = FALSE)
write.table(as.data.frame(colnames(atacMAT)), file=paste0('~/Desktop/STATS254/mouseData/scatacAS.binarized.genes','.tsv'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
rm(atacMAT, cLabs, tissIDs, cellIDs)

# Combine RNA+ATAC into a single matrix .....
rnaMAT <- Matrix::readMM('~/Desktop/STATS254/mouseData/scrna10X.binarized.mtx.gz')
atacCols <- read.delim('./scatacAS.binarized.genes.tsv',
                       header = FALSE,
                       sep = '\t',
                       stringsAsFactors = FALSE)
rnaCols <- read.delim('./scrna10X.binarized.genes.tsv',
                       header = FALSE,
                       sep = '\t',
                       stringsAsFactors = FALSE)
colnames(rnaMAT) <- rnaCols$V1
rnaMAT <- rnaMAT[, atacCols$V1]
atacMAT <- Matrix::readMM('~/Desktop/STATS254/mouseData/scatacAS.binarized.mtx.gz')
colnames(atacMAT) <- atacCols$V1

combMAT <- rbind(atacMAT, rnaMAT)
dim(combMAT)
head(combMAT)

Matrix::writeMM(combMAT, file=paste0('./scATACRNA.binarized','.mtx.gz'))
write.table(as.data.frame(colnames(combMAT)), file=paste0('./scATACRNA.binarized.genes','.tsv'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
rm(rnaMAT, atacMAT, combMAT)

atacAnnot <- read.delim('./scatacAS.binarized.annot.tsv',
                        sep = '\t',
                        stringsAsFactors = FALSE)
rnaAnnot <- read.delim('./scrna10X.binarized.annot.tsv',
                       sep = '\t',
                       stringsAsFactors = FALSE)
combAnnot <- rbind(atacAnnot, rnaAnnot)
dim(combAnnot)
head(combAnnot)
tail(combAnnot)
write.table(combAnnot, file=paste0('./scATACRNA.binarized.annot','.tsv'), quote=FALSE, sep='\t', row.names = FALSE)
rm(combAnnot, atacAnnot, rnaAnnot, rnaCols, atacCols)




