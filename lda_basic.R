# wrapLDA docs: https://cran.r-project.org/web/packages/text2vec/text2vec.pdf
library(Matrix)
library(text2vec)
library(data.table)
library(ggplot2)
library(tidyverse)

# Load Data
epimat_bin <- Matrix::t(Matrix::readMM('~/Desktop/STATS254/mouseData/activity_scores.binarized.mtx.gz'))*1 #multiply by 1 to convert ngTMatrix --> dgCMatrix
gene.names = read.delim('~/Desktop/STATS254/mouseData/activity_scores.binarized.genes.txt',
                        header = FALSE,
                        stringsAsFactors = FALSE)
cell.names = read.delim('~/Desktop/STATS254/mouseData/activity_scores.binarized.cells.txt',
                        header = FALSE,
                        stringsAsFactors = FALSE)
colnames(epimat_bin) <- gene.names$V1
rownames(epimat_bin) <- cell.names$V1
head(epimat_bin)


# Build model
topics <- 10
lda_model <- LDA$new(n_topics = topics, doc_topic_prior = 50/topics, topic_word_prior = 0.1)
doc_topic_distr <- lda_model$fit_transform(x = epimat_bin, 
                                           n_iter = 100,
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
model$perplexity <- perplexity(epimat_bin, lda_model$topic_word_distribution, doc_topic_distr)


# Plot Results
print(model$perplexity)
plot(model$log.likelihoods[1,], model$log.likelihoods[2,])
# TxC heatmap
dat <- model$topic_doc_distr
rownames(dat) <- paste0("Topic",seq(1,10))
colnames(dat) <- paste0(seq(1,81173))
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
rownames(dat) <- paste0("Topic",seq(1,10))
colnames(dat) <- paste0(seq(1,20783))
dat %>% 
  as.data.frame() %>%
  rownames_to_column("Topic") %>%
  pivot_longer(-c(Topic), names_to = "Gene", values_to = "porb") %>%
  ggplot(aes(x=Topic, y=Gene, fill=porb)) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_raster() +
  scale_fill_viridis_c()
# TxG count heatmap
dat <- model$topics
rownames(dat) <- paste0("Topic",seq(1,10))
colnames(dat) <- paste0(seq(1,20783))
dat %>% 
  as.data.frame() %>%
  rownames_to_column("Topic") %>%
  pivot_longer(-c(Topic), names_to = "Gene", values_to = "count") %>%
  ggplot(aes(x=Topic, y=Gene, fill=porb)) + 
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_raster() +
  scale_fill_viridis_c()
dat %>% 
  as.data.frame() %>%
  rownames_to_column("Topic") %>%
  pivot_longer(-c(Topic), names_to = "Gene", values_to = "count") %>%
  ggplot(aes(x=Topic, y=count)) +
  geom_bar(stat="identity") 
  
  