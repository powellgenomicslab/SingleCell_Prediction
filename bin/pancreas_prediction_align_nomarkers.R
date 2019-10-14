# Load libraries ----------------------------------------------------------
library("tidyverse")
library("scPred")
library("here")
library("cowplot")
library("Seurat")
library("viridis")
library("ggridges")

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

marker_genes <- c("INS", "GCG", "SST" ,"PPY")

removeMarkers <- function(x){
i <- which(row.names(x) %in% marker_genes)
x[-i,]
}

training <- lapply(training, removeMarkers)



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


saveRDS(pancreas.train.test, file = here(output, "train_test_alignment_30pcs_nomarkers.RDS"))
# pancreas.train.test <- readRDS(file = here("train_test_alignment.RDS"))


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
scpred <- trainModel(object = scpred, seed = 66, method = "svmRadial")

saveRDS(scpred, file = here(output, "trained_model_svmRadial_30PCs_nomarkers.RDS"))


# Classify cells in new dataset -------------------------------------------
predictions <- eigenPredict(scpred, test_eigen, threshold = 0.9)


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
y <- predictions
y$true <- test_metadata$cell_type1

library(tidyverse)

y %>% 
  filter(probability > 0.9) %>% 
  ggplot(aes(y = class, x = probability, fill = true)) +
  geom_density_ridges()



  geom_boxplot()

getAccuracy <- function(predictions, metadata){
  res <- data.frame(pred = predictions$class, true = metadata$cell_type1, 
                    row.names = rownames(predictions))
  
  res %>% 
    mutate(label = if_else(pred == as.character(true), as.character(true),
                           if_else(pred ==  "Unassigned" & true == "other", 
                                   "other", "Misclassified")))
  
}


table(x[, c("true", "pred")])



getAccuracy(muraro_predictions, muraro_metadata)

table(x[as.character(x$true) != as.character(x$pred),"pred"])

table(x$label)



res <- data.frame(pred = muraro_predictions$class, true = muraro_metadata$cell_type1, 
                  row.names = rownames(muraro_predictions))


res %>% 
  mutate(true = if_else(as.character(true) %in% c("alpha", "beta", "delta", "gamma"), 
                        as.character(true), "other")) %>% 
  mutate(label = if_else(pred == as.character(true), as.character(true),
                         if_else(pred ==  "Unassigned" & true == "other", 
                                 "other", "Misclassified"))) -> x
group_by(true, label) %>% 
  summarise(n = n()) %>% 
  mutate(accuracy = (n / sum(n))*100)


table(x[x$true == "other","pred"])









res %>% 
  mutate(true = if_else(as.character(true) %in% c("alpha", "beta", "delta", "gamma"), 
                        as.character(true), "other" )) %>% 
  mutate(label = if_else(pred == as.character(true), as.character(true),
                         if_else(pred ==  "Unassigned" & true == "other", 
                                 "other", "Misclassified"))) %>% 
  group_by(true, label) %>% 
  summarise(n = n()) %>% 
  mutate(accuracy = (n / sum(n))*100) %>% 
  filter(label != "Misclassified") %>% # accuracy per group
  pull(n) %>% 
  sum() -> x


x / nrow(res) # overall accuracy



res %>% 
  mutate(true = if_else(as.character(true) %in% c("alpha", "beta", "delta", "gamma"), 
                        as.character(true), "other" )) %>% 
  group_by(true, pred) %>% 
  summarise(n = n()) %>% 
  mutate(accuracy = (n / sum(n))*100)



# Get number of pcs and variance explained -------------------------------

pancreas@features %>% 
  Reduce(rbind, .) %>% 
  pull(PC) %>% 
  unique() %>% 
  as.character() -> unique_pcs
length(unique_pcs)
sum(pancreas@expVar[names(pancreas@expVar) %in% unique_pcs])
pancreas@features %>% 
  lapply(nrow)

# Get number of support vectors -------------------------------------------

pancreas@train %>% 
  lapply("[", "finalModel") %>% 
  unlist() %>% 
  lapply("slot", "nSV")


# Intersection of support cells
pancreas@train %>% 
  lapply("[", "finalModel") %>% 
  unlist() %>% 
  lapply("slot", "SVindex") %>% 
  unlist() %>% 
  unique() %>% 
  length()







res$human <- as.factor(baron_metadata$human)

res$human %>% table()
res %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  group_by(true, human) %>% 
  summarise(n = n())


res %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  group_by(pred, human) %>% 
  summarise(n = n())


res$lib <- baron_metadata$lib

res %>% 
  mutate(label = if_else(pred == as.character(true), as.character(true), "Misclassified")) %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  group_by(true, label, human, lib) %>% 
  summarise(prediction = n()) %>% 
  ungroup() %>% 
  group_by(true, human, lib) %>% 
  mutate(total = sum(prediction)) %>% 
  mutate(accuracy = prediction/total) %>% 
  filter(label != "Misclassified") %>% 
  ungroup() %>% 
  select(-true) %>% 
  View()




x <- cbind(res, baron_predictions)


x %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  ggplot(aes(y = probability, x = true, fill = human)) +
  geom_boxplot(alpha = 0.3) + 
  geom_hline(yintercept = 0.70) +
  facet_wrap(~lib)

x %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  ggplot(aes(y = probability, x = true, fill = lib)) +
  geom_boxplot(alpha = 0.3) + 
  geom_hline(yintercept = 0.70) +
  facet_wrap(~human)



#########

# Get gamma cells from reference
i <- reference_metadata$cellType == "gamma"
ref_gamma_metadata <- reference_metadata[i,,drop = FALSE]
ref_gamma <- reference[,i]
all(colnames(ref_gamma) == rownames(ref_gamma_metadata))

# Get gamma cells from test
i <- baron_metadata$cell_type1 == "gamma"
test_gamma_metadata <- as.data.frame(baron_metadata[i,])
test_gamma <- baron_cpm[,i]
all(colnames(test_gamma) == rownames(test_gamma_metadata))

ref_gamma <- ref_gamma[apply(ref_gamma, 1, var) != 0, ]
ref_gamma_cor <- cor(t(ref_gamma))
ref_gamma_cor <- abs(ref_gamma_cor)


test_gamma <- test_gamma[apply(test_gamma, 1, var) != 0, ]
test_gamma_cor <- cor(t(test_gamma))
test_gamma_cor <- abs(test_gamma_cor)


x <- apply(ref_gamma_cor, 1, 
           function(x){ 
             unlist(lapply(x, function(y){if(y >= 0.8){ 1}else{0}}))
           }
)

y <- apply(test_gamma_cor, 1, 
           function(x){ 
             unlist(lapply(x, function(y){if(y >= 0.8){ 1}else{0}}))
           }
)



<- x %*% diag(0,  nrow(x), ncol(x))

