# Load libraries ----------------------------------------------------------
library("tidyverse")
library("scPred")
library("here")
library("cowplot")
library("Seurat")
library("viridis")

output <- file.path("results", "2018-05-19_pancreas_prediction_align")

# Read datasets -----------------------------------------------------------

# Read training
training <- readRDS(here("results", "2018-05-12_pancreas_processing", 
                         "pancreas_cpm_train.RDS"))
training_metadata <- readRDS(here("results", "2018-05-12_pancreas_processing", 
                                  "pancreas_metadata_train.RDS"))

# Read testing
test <- readRDS(here("results", "2018-05-12_pancreas_processing",
                     "baron_cpm_test.RDS"))
test_metadata <- readRDS(here("results", "2018-05-12_pancreas_processing",
                              "baron_metadata_test.RDS"))


# Pre-process training datasets -------------------------------------------

normalize_data <- function(x, meta, label){
  x <- CreateSeuratObject(raw.data = x, 
                     meta.data = meta)
  x <- NormalizeData(x)
  x <- FindVariableGenes(x, do.plot = F, display.progress = F)
  x@meta.data$batch <- label
  x@meta.data$cell_type1 <- x@meta.data$x.cell_type.i.
  x
}


training_seurat <- mapply(normalize_data, training, training_metadata, names(training))

assignData <- function(x){
  ScaleData(x)
}

training_seurat %>% 
  lapply(assignData) -> training_seurat


# Pre-process test dataset ------------------------------------------------

baron <- CreateSeuratObject(raw.data = test,
                          meta.data = test_metadata)
baron <- NormalizeData(baron)
baron <- FindVariableGenes(baron, do.plot = F, display.progress = F)
baron@meta.data$batch <- "baron"


# Determine genes to use for PCA, must be highly variable in at least 2 datasets

genes.use <- c()
for (i in 1:length(training_seurat)) {
  genes.use <- c(genes.use, head(rownames(training_seurat[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(training_seurat)) {
  genes.use <- genes.use[genes.use %in% rownames(training_seurat[[i]]@scale.data)]
}



# Run PCA -----------------------------------------------------------------

pancreas.integrated <- MergeSeurat(training_seurat$muraro, training_seurat$segerstolpe)
pancreas.integrated <- MergeSeurat(pancreas.integrated, training_seurat$xin)
pancreas.integrated <- ScaleData(pancreas.integrated)

pancreas.integrated <- RunPCA(pancreas.integrated, 
                              pc.genes = genes.use, 
                              pcs.compute = 30)


# Test dataset projection -------------------------------------------------

# Get loadings gene expression values of variable genes in training dataset
loadings <- pancreas.integrated@dr$pca@gene.loadings
genes <- pancreas.integrated@scale.data[rownames(pancreas.integrated@scale.data) %in% genes.use, ]

# Center and scale test dataset 
new.center <- rowMeans(genes)
new.scale <- apply(genes, 1, sd)

i <- match(rownames(loadings), names(new.center))
new.center <- new.center[i]

i <- match(rownames(loadings), names(new.scale))
new.scale <- new.scale[i]

share_genes <- intersect(rownames(baron@raw.data), rownames(loadings))

baron_sub <- baron@data[rownames(baron@data) %in% share_genes, ]
loadings_sub <- loadings[rownames(loadings) %in% share_genes,]
new.center <- new.center[rownames(loadings) %in% share_genes]
new.scale <- new.scale[rownames(loadings) %in% share_genes]
all(rownames(baron_sub) == rownames(loadings_sub))

# Normalization and projection
projection <- scale(t(baron_sub), center = new.center, scale = new.scale) %*% loadings_sub


# Merge training eigenspace and projection
test_projection <- rbind(pancreas.integrated@dr$pca@cell.embeddings, projection)

# Assign projection to training object
pancreas.integrated@dr$pca@cell.embeddings <- as.matrix(test_projection)

# Add test metadata to training object
trainInfo <- pancreas.integrated@meta.data[, c("batch", "cell_type1")]
testInfo <- data.frame(batch = "baron", 
                       cell_type1 = baron@meta.data$cell_type1, 
                       row.names = rownames(baron@meta.data))

pancreas.integrated@meta.data <- rbind(trainInfo, testInfo)

pancreas.integrated@meta.data$batch <- as.factor(pancreas.integrated@meta.data$batch)
all(rownames(pancreas.integrated@meta.data)  == rownames(pancreas.integrated@dr$pca@cell.embeddings))


unaligned_test_train <- as.data.frame(pancreas.integrated@dr$pca@cell.embeddings)

var_exp <- (pancreas.integrated@dr$pca@sdev**2/sum(pancreas.integrated@dr$pca@sdev**2))*100


cell_labels <- ifelse(pancreas.integrated@meta.data$cell_type1 %in% c("alpha", "beta", "delta", "gamma"),
                      as.character(pancreas.integrated@meta.data$cell_type1),
                      "other")

unaligned_test_train$cellType <- factor(cell_labels, 
                                        levels = c("alpha", "beta", "delta", "gamma", "other"),
                                        labels = c("Alpha", "Beta", "Delta", "Gamma", "Other"))
unaligned_test_train$batch <- factor(pancreas.integrated@meta.data$batch, 
                                     levels = c("muraro", "segerstolpe", "xin", "baron"),
                                     labels = c("Muraro", "Segerstolpe", "Xin", "Baron"))

unaligned_test_train %>% 
  ggplot(aes(x = PC1, y = PC2, col = batch)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, direction = -1) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  xlab(paste0("PC1 ", "(", round(var_exp[1], 2), "% exp. var.)" ,collapse = "")) + 
  ylab(paste0("PC2 ", "(", round(var_exp[2], 2), "% exp. var.)" ,collapse = "")) -> p1

unaligned_test_train %>% 
  ggplot(aes(x = PC1, y = PC2, col = cellType)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.title=element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  xlab(paste0("PC1 ", "(", round(var_exp[1], 2), "% exp. var.)" ,collapse = "")) + 
  ylab(paste0("PC2 ", "(", round(var_exp[2], 2), "% exp. var.)" ,collapse = "")) -> p2





# Merge training and test datasets in a single Seurat object --------------

pancreas.train.test <- MergeSeurat(pancreas.integrated, baron)
pancreas.train.test@meta.data <- pancreas.integrated@meta.data
pancreas.train.test@dr$pca <- pancreas.integrated@dr$pca



pancreas.train.test <- AlignSubspace(pancreas.train.test,
                                     reduction.type = "pca",
                                     grouping.var = "batch",
                                     dims.align = 1:30)


saveRDS(pancreas.train.test, file = here(output, "train_test_alignment_30pcs.RDS"))
# pancreas.train.test <- readRDS(file = here(output, "train_test_alignment_30pcs.RDS"))


pancreas.aligned <- as.data.frame(pancreas.train.test@dr$pca.aligned@cell.embeddings)

cell_labels <- ifelse(pancreas.train.test@meta.data$cell_type1 %in% c("alpha", "beta", "delta", "gamma"),
                      as.character(pancreas.train.test@meta.data$cell_type1),
                      "other")
pancreas.aligned$cellType <- factor(cell_labels, levels = c("alpha", "beta", "delta", "gamma", "other"))

pancreas.aligned$batch <- factor(pancreas.train.test@meta.data$batch, 
       levels = c("muraro", "segerstolpe", "xin", "baron"),
       labels = c("Muraro", "Segerstolpe", "Xin", "Baron"))


pancreas.aligned %>% 
  ggplot(aes(x = APC1, y = APC2, col = batch)) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, direction = -1) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3))) -> p3


pancreas.aligned %>% 
  ggplot(aes(x = APC1, y = APC2, col = cellType)) +
  geom_point(alpha = 0.2) +
  geom_point(size = 0.6, alpha = 0.8) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.title=element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3))) -> p4



# Get row 1
batch_legend <- get_legend(p1)
row1 <- plot_grid(p1 + theme(legend.position="none"),
          p3 + theme(legend.position="none"))
row1 <- plot_grid(row1, batch_legend,  rel_widths = c(2, .5))


# Get row 2

cellType_legend <- get_legend(p2)
row2 <- plot_grid(p2 + theme(legend.position="none"),
                  p4 + theme(legend.position="none"))
row2 <- plot_grid(row2, cellType_legend,  rel_widths = c(2, .5))

# Bind rowa

pca_results <- plot_grid(row1, row2, 
                         nrow = 2, 
                         align = "v", 
                         labels = c("a", "b"), label_size = 14)

ggsave(here(output, "pancreas_pca_results.png"), plot = pca_results, 
       width = 7, height = 5, dpi = 250)



all(pancreas.train.test@dr$pca@cell.embeddings == test_projection)



# Make predictions --------------------------------------------------------

# Create training and test datasets
i <- pancreas.train.test@meta.data$batch == "baron"
train_eigen <- pancreas.train.test@dr$pca.aligned@cell.embeddings[!i, ]
test_eigen <- pancreas.train.test@dr$pca.aligned@cell.embeddings[i, ]


train_metadata <- pancreas.train.test@meta.data[!i, ]
train_metadata$cell_type1 <- factor(train_metadata$cell_type1, levels = unique(train_metadata$cell_type1))

test_metadata <- pancreas.train.test@meta.data[i, ]


# Get informative principal components ------------------------------------
inf_pcs <- getInformativePCs(as.data.frame(train_eigen), 
                             classes = train_metadata$cell_type1, 
                             "cell_type1", 
                             sig = 0.05)



scpred <- list(prcomp = train_eigen, 
               metadata = train_metadata,
               pVar = "cell_type1",
               features = inf_pcs)



# Train model -------------------------------------------------------------
scpred <- trainModel(object = scpred, seed = 66, method = "svmRadial", savePredictions = TRUE)
predictions <- eigenPredict(scpred, train_eigen, threshold = 0.9)
# predictions_train <- eigenPredict(scpred, train_eigen, threshold = 0.9)


saveRDS(scpred, file = here(output, "trained_model_glm_30PCs.RDS"))
# scpred <- readRDS(file = here(output, "trained_model_svmRadial_30PCs.RDS"))


# Classify cells in new dataset -------------------------------------------
predictions <- eigenPredict(scpred, test_eigen, threshold = 0.9)

write.csv(predictions, 
          file = here(output, "pancreas_predictions.csv"), 
          quote = FALSE)


predictions$true <- test_metadata$cell_type1


getMisClass <- function(x){
predictions %>% 
  filter(true == x, class != x) %>% 
  pull(class) %>% 
  table()
}

lapply(c("alpha", "beta", "delta", "gamma"), getMisClass)


predictions %>% 
  filter(!true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  pull(class) %>% 
  table()



# Get accuracy results ----------------------------------------------------

getAccuracy <- function(predictions, metadata){
  res <- data.frame(pred = predictions$class, true = test_metadata$cell_type1, 
                    row.names = rownames(predictions))
  
  res %>% 
    mutate(true = if_else(as.character(true) %in% c("alpha", "beta", "delta", "gamma"), 
                          as.character(true), "other" )) %>% 
    mutate(label = if_else(pred == as.character(true), as.character(true),
                           if_else(pred ==  "Unassigned" & true == "other", 
                                   "other", "Misclassified"))) %>% 
    group_by(true, label) %>% 
    summarise(n = n()) %>% 
    mutate(accuracy = (n / sum(n))*100) %>% 
    filter(label != "Misclassified") 
}


getAccuracy(predictions, test_metadata)


predictions$true <- test_metadata$cell_type1
# predictions_train$true <- train_metadata$cell_type1


predictions %>% 
  mutate(true = ifelse(true %in% c("alpha", "beta", "delta", "gamma"),
                                    as.character(true), "other")) %>% 
  mutate(true = factor(true, 
                       levels = c("alpha", "beta", "delta", "gamma", "other"), 
                       labels = c("Alpha", "Beta", "Delta", "Gamma", "Other") )) %>% 
  mutate(class = factor(class, 
                        levels = c("alpha", "beta", "delta", "gamma", "Unassigned"), 
                        labels = c("Alpha", "Beta", "Delta", "Gamma", "Unassigned"))) %>% 
  ggplot() +
  aes(x = probability, fill = class) +
  geom_histogram(color = "black", alpha = 0.8) +
  facet_wrap(~true, scale = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  xlab("Probability") +
  ylab("Counts") +
  theme_bw() +
  theme(strip.text.x = element_text(size = 14), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.title = element_text(size = 16), 
        strip.background = element_rect(fill = "grey95")) +
  guides(fill = guide_legend(title = "Predicted class")) -> p_probs


ggsave(here(output, "pancreas_probs.png"), 
       plot = p_probs, width = 10, height = 4.5, dpi = 350)


## Plot roc curves for training dataset

devtools::install_github("sachsmc/plotROC")
library(plotROC)

scpred_beta <- scpred
scpred_beta$features[c("alpha", "delta", "gamma")] <- NULL
scpred_beta$train[c("alpha", "delta", "gamma")] <- NULL
scpred_beta$metadata$cell_type1 <- factor(ifelse(scpred_beta$metadata$cell_type1 == "beta", "beta", "other"),
                                          levels = c("beta", "other"))
predictions_train <- eigenPredict(scpred_beta, train_eigen, threshold = 0.9)
predictions_train$true <- scpred_beta$metadata$cell_type1

predictions_train %>% 
  filter(true == "beta") %>% 
  ggplot() +
  aes(x = probability, fill = true) +
  geom_histogram() +
  theme_bw() +
  geom_vline(xintercept = 0.9)

predictions_train %>% 
  mutate(true = ifelse(true == "beta", 1, 0)) %>% 
 ggplot(aes(d = true, m = probability, fill = true)) + 
  geom_roc() +
   style_roc()



plot(roc(predictions_train$true, 
         predictions_train$probability), 
     col = predictions_train$true, print.thres = c(1.4))



plot.roc(scpred$train$beta$pred$obs,
         scpred$train$beta$pred$pred)


library(precrec)

# Load a test dataset
data(P10N10)

# Calculate ROC and Precision-Recall curves
sscurves <- evalmod(scores = P10N10$scores, labels = P10N10$labels)


sscurves <- evalmod(scores = predictions_train$probability, labels = predictions_train$true)

plot(sscurves)
