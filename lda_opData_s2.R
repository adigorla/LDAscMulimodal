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

# Load data & build model
combMAT <- Matrix::readMM(paste0('./opmultiome_GEXATAC_s2.binarized','.mtx.gz'))*1
head(combMAT)
gene.names = read.delim(paste0('./opmultiome_GEXATAC_s2.quantitative.features','.tsv'),
                        header = TRUE,
                        stringsAsFactors = FALSE)
annots = read.delim(paste0('./opmultiome_GEXATAC_s2.quantitative.annot','.tsv'),
                    header = TRUE,
                    stringsAsFactors = FALSE)
colnames(combMAT) <- gene.names$featID
rownames(combMAT) <- annots$cellLabel

topics <- 25
lda_model <- LDA$new(n_topics = topics, doc_topic_prior = 50/topics, topic_word_prior = 0.1)
doc_topic_distr <- lda_model$fit_transform(x = combMAT, 
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
model$perplexity <- perplexity(combMAT, lda_model$topic_word_distribution, doc_topic_distr)

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
pdf("./s2cplxHeatmap.pdf", width = 16, height = 8)
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
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint, shape = factor(Batch))) +
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




# Topic-CT proportions
dat <- model$topic_doc_distr
rownames(dat) <- paste0("Topic_",seq(1,topics))
colnames(dat) <- annots$cellID
dat <- as.data.frame(t(dat))
dat$Type <- annots$cellLabel
dat[is.na(dat)] = 0
# dat <- dat %>% 
#   group_by(Type) %>% 
#   summarise(across(everything(), sum))
dat %>% 
  group_by(Type) %>% 
  summarise(across(everything(), sum)) %>% 
  pivot_longer(!Type, names_to = c("Topic"), names_prefix = "Topic_",
               values_to = "porb", values_ptypes = list(prob = double())) %>%
  ggplot(aes(x=Type, y=porb, fill=Topic)) + 
  geom_bar(position="fill", stat="identity") + 
  labs(y = "probability Prop" ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6))


# Topic-Mode proportions
dat <- model$topic_word_distr
rownames(dat) <- paste0("Topic_",seq(1,topics))
dat <- as.data.frame(t(dat))
dat$Type <- gene.names$featType
dat[is.na(dat)] = 0
# dat <- dat %>% 
#   group_by(Type) %>% 
#   summarise(across(everything(), sum))

dat %>% 
  group_by(Type) %>% 
  summarise(across(everything(), sum)) %>% 
  pivot_longer(!Type, names_to = c("Topic"), names_prefix = "Topic_",
               values_to = "porb", values_ptypes = list(prob = double())) %>%
  ggplot(aes(x=Topic, y=porb, fill=Type)) + 
  geom_bar(position="fill", stat="identity") + 
  labs(y = "probability Prop" ) +
  theme(axis.text.x = element_text(vjust = 0.6))

# Weighted props ??

dat1 <- model$topic_doc_distr
rownames(dat1) <- paste0("Topic_",seq(1,topics))
colnames(dat1) <- annots$cellID
dat1 <- as.data.frame(t(dat1))
dat1$Type <- annots$cellLabel
dat1 <- dat1 %>% group_by(Type) %>% 
  summarise(across(everything(), sum))

dat2 <- model$topic_word_distr
rownames(dat2) <- paste0("Topic_",seq(1,topics))
dat2 <- as.data.frame(t(dat2))
dat2$Type <- gene.names$featType
dat2[is.na(dat)] = 0
dat2 <- dat2 %>% group_by(Type) %>% 
  summarise(across(everything(), sum))

dat <- as.matrix(dat1[, !names(dat1) %in% c("Type")]) %*% t(as.matrix(dat2[, !names(dat1) %in% c("Type")]))
colnames(dat) <- dat2$Type
dat <- as.data.frame(dat)
dat$Cell <- dat1$Type
dat %>% 
  pivot_longer(!Cell, names_to = c("SRC"),
               values_to = "porb", values_ptypes = list(prob = double())) %>%
  ggplot(aes(x=Cell, y=porb, fill=SRC)) + 
  geom_bar(position="fill", stat="identity") + 
  labs(y = "Weight" ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6))
#only 10.33% of the features are GEX