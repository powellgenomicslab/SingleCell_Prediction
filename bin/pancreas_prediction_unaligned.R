# Load libraries ----------------------------------------------------------
library("tidyverse")
library("scPred")
library("here")
library("cowplot")
library("Seurat")
library("viridis")

output <- file.path("results", "2019-03-03_pancreas_prediction_unalign")

if(!dir.exists(output)){
  dir.create(here(output))
}

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



PCAPlot(pancreas.integrated, group = "batch")

# Get informative principal components ------------------------------------
model <- getFeatureSpace(pancreas.integrated, pVar = "cell_type1")


# Train model -------------------------------------------------------------
model <- trainModel(object = model)

saveRDS(model, file = here(output, "model.RDS"))


# Classify cells in new dataset -------------------------------------------
model <- scPredict(model, baron)
model@predMeta <- as.data.frame(baron@meta.data)
model %>% getPredictions()

orderLab <- c("alpha", "beta",  "delta", "gamma",  "acinar", "activated_stellate",  "ductal", "endothelial",  "epsilon", "macrophage",  "mast", "quiescent_stellate",  "schwann", "t_cell")

accuracy <- crossTab(model, true = "cell_type1")[,orderLab]


## Overall sensitivity
mean(c(0.88, 0.7, 0.64, 0.43))

## Overall specificity
accuracy[5,] %>% as.numeric() %>% mean()

