library(Matrix)
library(text2vec)
library(data.table)
library(ggplot2)
library(tidyverse)
library(fastcluster)
library(devtools)
library(ComplexHeatmap)
library(scales)
library(Rtsne)
library(umap)
library(destiny)

setwd("~/Desktop/STATS254/opMultiome10x")

# Load GEX data & build model
combMAT <- Matrix::readMM(paste0('./opmultiome_GEXATAC.binarized','.mtx.gz'))*1
gene.names = read.delim(paste0('./opmultiome_GEXATAC.quantitative.features','.tsv'),
                        header = TRUE,
                        stringsAsFactors = FALSE)
gex.locs <- c(1:length(gene.names$featType))[gene.names$featType == "GEX"]
gexMAT <- combMAT[,gex.locs]
dim(gexMAT)
head(gexMAT)
rm(combMAT)
gene.names <- gene.names[gex.locs,]
annots = read.delim(paste0('./opmultiome_GEXATAC.quantitative.annot','.tsv'),
                    header = TRUE,
                    stringsAsFactors = FALSE)
colnames(gexMAT) <- gene.names$featID
rownames(gexMAT) <- annots$cellLabel

topics <- 25
lda_model <- LDA$new(n_topics = topics, doc_topic_prior = 50/topics, topic_word_prior = 0.1)
doc_topic_distr <- lda_model$fit_transform(x = gexMAT, 
                                           n_iter = 150,
                                           convergence_tol = 0.0001, 
                                           n_check_convergence = 4,
                                           progressbar = TRUE, normalize='none')
model <- list()
model$topics <- lda_model$components
model$topic_sums <- as.matrix(rowSums(lda_model$components))
#model$document_sums <- t(doc_topic_distr)
model$topic_doc_distr <- t(doc_topic_distr)
model$topic_word_distr <- lda_model$topic_word_distribution
model$log.likelihoods <- t(attributes(doc_topic_distr)$likelihood)
model$perplexity <- perplexity(gexMAT, lda_model$topic_word_distribution, doc_topic_distr)

# Plot Results
print(model$perplexity)
plot(model$log.likelihoods[1,], model$log.likelihoods[2,])
# TxC heatmap
dat <- model$topic_doc_distr
rownames(dat) <- paste0("Topic",seq(1,topics))
colnames(dat) <- annots$cellID
dat %>% 
  as.data.frame() %>%
  rownames_to_column("Topic") %>%
  pivot_longer(-c(Topic), names_to = "Sample", values_to = "porb") %>%
  ggplot(aes(x=Topic, y=Sample, fill=porb)) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_raster() +
  scale_fill_viridis_c()

# TxG heatmap
dat <- model$topic_word_distr
rownames(dat) <- paste0("Topic",seq(1,topics))
colnames(dat) <- gene.names$featID
dat %>% 
  as.data.frame() %>%
  rownames_to_column("Topic") %>%
  pivot_longer(-c(Topic), names_to = "Gene", values_to = "porb") %>%
  ggplot(aes(x=Topic, y=Gene, fill=porb)) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_raster() +
  scale_fill_viridis_c()


# Complex TxC heatmap
col.low = "floralwhite"
col.mid = "pink"
col.high = "red"
distinctColorPalette <-function(k) {
  ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=85)(2e3))))
  km <- kmeans(ColorSpace, k, iter.max=20)
  colors <- rgb(round(km$centers), maxColorValue=255)
  return(colors)
}

dat <- as.matrix(model$topic_doc_distr)
rownames(dat) <- paste0("Topic",seq(1,topics))
colnames(dat) <- annots$cellLabel

cl.cells <- fastcluster::hclust.vector(t(dat), method="ward", metric="euclidean")
dd.cells <- as.dendrogram(cl.cells)
colorPal <- colorRampPalette(c(col.low, col.mid, col.high))

colorBy <- c( 'cellLabel', 'batch', 'phase')
colVars <- list()
for (variable in colorBy){
  colVars[[variable]] <- setNames(distinctColorPalette(length(unique(annots[variable])[,1])), as.vector(sort(unique(annots[variable])[,1])))
  #cellColor <- setNames(colVars[[variable]][annots[variable]], rownames(dat))
}

annotation <- ComplexHeatmap::HeatmapAnnotation(df = annots[,colorBy,drop=FALSE] , col = colVars, which='column')
heatmap <- ComplexHeatmap::Heatmap(dat, col=colorPal(20), cluster_columns=dd.cells, name='Probability', use_raster = TRUE,
                                   show_column_names=FALSE, show_row_names = TRUE, top_annotation = annotation, 
                                   column_order = sort(colnames(dat)),
                                   heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"), 
                                                               column_title = "Topic contribution per cell", column_title_gp = gpar(fontface = 'bold')))
pdf("./cplxHeatmap_gexonly.pdf", width = 16, height = 8)
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()

# Run UMAP
Umap <- umap::umap(t(dat))
# rownames(Umap$layout) <- annots$cellID
# colnames(Umap$layout) <- paste0('UMAP', 1:ncol(Umap$layout))

ggplot(data=data.frame(UMAP1 = Umap$layout[,1], UMAP2 = Umap$layout[,2], cell = annots$cellLabel), 
       aes(x=UMAP1, y=UMAP2, color=cell)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

ggplot(data=data.frame(UMAP1 = Umap$layout[,1], UMAP2 = Umap$layout[,2], SRC = annots$batch), 
       aes(x=UMAP1, y=UMAP2, color=SRC)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

ggplot(data=data.frame(UMAP1 = Umap$layout[,1], UMAP2 = Umap$layout[,2], PH = annots$phase), 
       aes(x=UMAP1, y=UMAP2, color=PH)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))



# Run tSNE
tSNE <- Rtsne::Rtsne(t(dat))
#rownames(tSNE$Y) <- annots$cellID
#colnames(tSNE$Y) <- paste0('tSNE', 1:ncol(tSNE$Y))

ggplot(data=data.frame(x = tSNE$Y[,1], y = tSNE$Y[,2], col = annots$cellLabel), 
       aes(x=x, y=y, color=col)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

ggplot(data=data.frame(x = tSNE$Y[,1], y = tSNE$Y[,2], col = annots$batch), 
       aes(x=x, y=y, color=col)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

ggplot(data=data.frame(x = tSNE$Y[,1], y = tSNE$Y[,2], col = annots$phase), 
       aes(x=x, y=y, color=col)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

# Diffusion maps
dm <- destiny::DiffusionMap(t(dat))
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Timepoint = annots$phase,
                  Batch = annots$batch)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Batch)) +
  geom_point() +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  cell = annots$cellLabel)
ggplot(tmp, aes(x = DC1, y = DC2, colour = cell)) +
  geom_point()+
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

ggplot(tmp, aes(x = DC1, y = cell, colour = cell)) +
    geom_point()+
  xlab("Diffusion component 1") + 
    ylab("Cell Type") +
    theme_classic()
plot(eigenvalues(dm), ylim = 0:1, pch = 20, xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')
# shit eigval dropoff so results not very trsutworthy! 
dpt <- DPT(dm)
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  cell = annots$cellLabel,
                  pseudotime_dpt = rank(dpt$dpt))
ggplot(tmp, aes(x = pseudotime_dpt, y = cell, colour = cell)) +
  geom_point() +
  xlab("Diffusion component 1") + 
  ylab("Cell Type") +
  theme_classic()






##############################################################################################################################################################################################################################################
rm(list = ls())



# Load ATAC data & build model
combMAT <- Matrix::readMM(paste0('./opmultiome_GEXATAC.binarized','.mtx.gz'))*1
gene.names = read.delim(paste0('./opmultiome_GEXATAC.quantitative.features','.tsv'),
                        header = TRUE,
                        stringsAsFactors = FALSE)
atac.locs <- c(1:length(gene.names$featType))[gene.names$featType == "ATAC"]
atacMAT <- combMAT[,atac.locs]
dim(atacMAT)
head(atacMAT)
rm(combMAT)
gene.names <- gene.names[atac.locs,]
annots = read.delim(paste0('./opmultiome_GEXATAC.quantitative.annot','.tsv'),
                    header = TRUE,
                    stringsAsFactors = FALSE)
colnames(atacMAT) <- gene.names$featID
rownames(atacMAT) <- annots$cellLabel

topics <- 25
lda_model <- LDA$new(n_topics = topics, doc_topic_prior = 50/topics, topic_word_prior = 0.1)
doc_topic_distr <- lda_model$fit_transform(x = atacMAT, 
                                           n_iter = 150,
                                           convergence_tol = 0.0001, 
                                           n_check_convergence = 4,
                                           progressbar = TRUE, normalize='none')
model <- list()
model$topics <- lda_model$components
model$topic_sums <- as.matrix(rowSums(lda_model$components))
#model$document_sums <- t(doc_topic_distr)
model$topic_doc_distr <- t(doc_topic_distr)
model$topic_word_distr <- lda_model$topic_word_distribution
model$log.likelihoods <- t(attributes(doc_topic_distr)$likelihood)
model$perplexity <- perplexity(atacMAT, lda_model$topic_word_distribution, doc_topic_distr)

# Plot Results
print(model$perplexity)
plot(model$log.likelihoods[1,], model$log.likelihoods[2,])
# TxC heatmap
dat <- model$topic_doc_distr
rownames(dat) <- paste0("Topic",seq(1,topics))
colnames(dat) <- annots$cellID
dat %>% 
  as.data.frame() %>%
  rownames_to_column("Topic") %>%
  pivot_longer(-c(Topic), names_to = "Sample", values_to = "porb") %>%
  ggplot(aes(x=Topic, y=Sample, fill=porb)) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_raster() +
  scale_fill_viridis_c()

# TxG heatmap
dat <- model$topic_word_distr
rownames(dat) <- paste0("Topic",seq(1,topics))
colnames(dat) <- gene.names$featID
dat %>% 
  as.data.frame() %>%
  rownames_to_column("Topic") %>%
  pivot_longer(-c(Topic), names_to = "Gene", values_to = "porb") %>%
  ggplot(aes(x=Topic, y=Gene, fill=porb)) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_raster() +
  scale_fill_viridis_c()


# Complex TxC heatmap
col.low = "floralwhite"
col.mid = "pink"
col.high = "red"
distinctColorPalette <-function(k) {
  ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=85)(2e3))))
  km <- kmeans(ColorSpace, k, iter.max=20)
  colors <- rgb(round(km$centers), maxColorValue=255)
  return(colors)
}

dat <- as.matrix(model$topic_doc_distr)
rownames(dat) <- paste0("Topic",seq(1,topics))
colnames(dat) <- annots$cellLabel

cl.cells <- fastcluster::hclust.vector(t(dat), method="ward", metric="euclidean")
dd.cells <- as.dendrogram(cl.cells)
colorPal <- colorRampPalette(c(col.low, col.mid, col.high))

colorBy <- c( 'cellLabel', 'batch', 'phase')
colVars <- list()
for (variable in colorBy){
  colVars[[variable]] <- setNames(distinctColorPalette(length(unique(annots[variable])[,1])), as.vector(sort(unique(annots[variable])[,1])))
  #cellColor <- setNames(colVars[[variable]][annots[variable]], rownames(dat))
}

annotation <- ComplexHeatmap::HeatmapAnnotation(df = annots[,colorBy,drop=FALSE] , col = colVars, which='column')
heatmap <- ComplexHeatmap::Heatmap(dat, col=colorPal(20), cluster_columns=dd.cells, name='Probability', use_raster = TRUE,
                                   show_column_names=FALSE, show_row_names = TRUE, top_annotation = annotation, 
                                   column_order = sort(colnames(dat)),
                                   heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(5, "cm"), 
                                                               column_title = "Topic contribution per cell", column_title_gp = gpar(fontface = 'bold')))
pdf("./cplxHeatmap_ataconly.pdf", width = 16, height = 8)
ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()

# Run UMAP
Umap <- umap::umap(t(dat))
# rownames(Umap$layout) <- annots$cellID
# colnames(Umap$layout) <- paste0('UMAP', 1:ncol(Umap$layout))

ggplot(data=data.frame(UMAP1 = Umap$layout[,1], UMAP2 = Umap$layout[,2], cell = annots$cellLabel), 
       aes(x=UMAP1, y=UMAP2, color=cell)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

ggplot(data=data.frame(UMAP1 = Umap$layout[,1], UMAP2 = Umap$layout[,2], SRC = annots$batch), 
       aes(x=UMAP1, y=UMAP2, color=SRC)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

ggplot(data=data.frame(UMAP1 = Umap$layout[,1], UMAP2 = Umap$layout[,2], PH = annots$phase), 
       aes(x=UMAP1, y=UMAP2, color=PH)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))



# Run tSNE
tSNE <- Rtsne::Rtsne(t(dat))
#rownames(tSNE$Y) <- annots$cellID
#colnames(tSNE$Y) <- paste0('tSNE', 1:ncol(tSNE$Y))

ggplot(data=data.frame(x = tSNE$Y[,1], y = tSNE$Y[,2], col = annots$cellLabel), 
       aes(x=x, y=y, color=col)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

ggplot(data=data.frame(x = tSNE$Y[,1], y = tSNE$Y[,2], col = annots$batch), 
       aes(x=x, y=y, color=col)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

ggplot(data=data.frame(x = tSNE$Y[,1], y = tSNE$Y[,2], col = annots$phase), 
       aes(x=x, y=y, color=col)) +
  geom_point() +
  guides(fill=guide_legend(ncol=10)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "cm"))

# Diffusion maps
dm <- destiny::DiffusionMap(t(dat))
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  Timepoint = annots$phase,
                  Batch = annots$batch)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Batch)) +
  geom_point() +
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  cell = annots$cellLabel)
ggplot(tmp, aes(x = DC1, y = DC2, colour = cell)) +
  geom_point()+
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

ggplot(tmp, aes(x = DC1, y = cell, colour = cell)) +
  geom_point()+
  xlab("Diffusion component 1") + 
  ylab("Cell Type") +
  theme_classic()
plot(eigenvalues(dm), ylim = 0:1, pch = 20, xlab = 'Diffusion component (DC)', ylab = 'Eigenvalue')
# shit eigval dropoff so results not very trsutworthy! 
dpt <- DPT(dm)
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  cell = annots$cellLabel,
                  pseudotime_dpt = rank(dpt$dpt))
ggplot(tmp, aes(x = pseudotime_dpt, y = cell, colour = cell)) +
  geom_point() +
  xlab("Diffusion component 1") + 
  ylab("Cell Type") +
  theme_classic()