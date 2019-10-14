# Script information ------------------------------------------------------
# title: Classification of dendritic cells and monocytes
# author: José Alquicira Hernández
# date: 2018-10-12
# description: None


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary
library("Seurat")
library("scPred")

# Set output --------------------------------------------------------------


output_dir_name <- "10x_immune_dcs_vs_monocytes" # <------ Output directory
date <- "2018-10-13_" # <------ Fix date

output <- file.path("results", paste0(date, output_dir_name))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input <- file.path("results", "2018-10-12_10x_immune_myeloid_vs_lymphoid_prediction") # <------ Input directory

# Read files

filename <- "predictions.RDS" # <------ Input file
predictions <- readRDS(here(input, filename))


input <- file.path("results", "2018-10-12_10x_immune_tree") # <------ Input directory


filename <- "train_myeloid.RDS" # <------ Input file
train_data <- readRDS(here(input, filename))


filename <- "pred_data.RDS" # <------ Input file
pred_data <- readRDS(here(input, filename))

# Subset cells

pred_data@meta.data$pred1 <- as.factor(predictions$predClass)
j <- rownames(pred_data@meta.data)[pred_data@meta.data$pred1 %in% c("Myeloid_cell")]
pred_data <- SubsetData(pred_data, cells.use = j, do.clean = TRUE)
pred_data@meta.data$pred1 <- factor(pred_data@meta.data$pred1)


# Process data ------------------------------------------------------------


rawData <- train_data@raw.data
i <- names(which(rowSums(rawData) != 0))
rawData <- rawData[i, ]

train_data <- CreateSeuratObject(rawData, 
                                 meta.data = train_data@meta.data,
                                 min.cells = round(ncol(rawData)*0.01), 
                                 min.genes = 200)

train_data <- train_data %>% 
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>% 
  # FindVariableGenes(mean.function = ExpMean, 
  #                   dispersion.function = LogVMR, 
  #                   x.low.cutoff = 0.02, 
  #                   x.high.cutoff = 5, 
  #                   y.cutoff = 0.5) %>% 
  ScaleData() %>% 
  RunPCA(pc.genes = rownames(.@data), 
         do.print = TRUE, 
         pcs.print = 1:5, 
         genes.print = 5)


# Mnocytes
FeaturePlot(train_data, features.plot = c("FTL"), reduction.use = "pca", cols.use = c("grey", "red"))

# Dendiritic cells
FeaturePlot(train_data, features.plot = c("CLEC9A"), reduction.use = "pca", cols.use = c("grey", "red"))
FeaturePlot(train_data, features.plot = c("CD1C"), reduction.use = "pca", cols.use = c("grey", "red"))

# Markers: Monocytes
FeaturePlot(train_data, features.plot = c("CD14"), reduction.use = "pca", cols.use = c("grey", "red"))


# Get feature space -------------------------------------------------------

sc_pred <- getFeatureSpace(train_data, pVar = "cellType1")
rm(train_data)


# Train models ------------------------------------------------------------

library(doParallel)
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
sc_pred <- trainModel(sc_pred)
stopCluster(cl)


# Make predictions --------------------------------------------------------
sc_pred <- scPredict(sc_pred, pred_data)
predictions <- getPredictions(sc_pred)

predictions$true <- factor(pred_data@meta.data$cellType1)

tidy_table <- function(x, left, top, fill = 0, prop = FALSE, digits = 2){
  x %>% 
    group_by_(left, top) %>% 
    summarise(n = n()) %>% 
    spread(key = top, value = "n", fill = fill) %>% 
    as.data.frame() %>% 
    column_to_rownames(left) -> x
  
  if(prop){
    row_names <- rownames(x)
    x <- mapply(function(x,d){x/d}, x, colSums(x))
    rownames(x) <- row_names
    x %>% 
      round(digits) %>% 
      as.data.frame() -> x
    
  }
  
  x
  
}

predictions <- predictions %>% 
  mutate(predClass = if_else(Dendritic_cell > 0.9, "DC", 
                             if_else(Dendritic_cell < 0.1, "Monocyte", "Unassigned")))

tidy_table(predictions, "predClass", "true") %>% 
  select(-B_cell, -Blood_progenitor, -T_cell, -Natural_killer)

tidy_table(predictions, "predClass", "true", prop = TRUE) %>% 
  select(-B_cell, -Blood_progenitor, -T_cell, -Natural_killer)



# Feature transformation --------------------------------------------------
sc_pred_ss <- sc_pred

featureList <- as.character(sc_pred@features$Dendritic_cell$PC)
featureSpace <- scPred:::subsetMatrix(sc_pred@svd$x, s = featureList)

svd_ss <- spatialSign(featureSpace)
sc_pred_ss@svd$x[,colnames(sc_pred_ss@svd$x) %in% colnames(svd_ss)] <- svd_ss



sc_pred_ss@projection <- matrix(nrow = 0, ncol = 0)
plotEigen(sc_pred_ss, group = "cellType1")

sc_pred_ss <- getFeatureSpace(sc_pred_ss, pVar = "cellType1")

cl <- makePSOCKcluster(2)
registerDoParallel(cl)
sc_pred_ss <- trainModel(sc_pred_ss, model = "earth")
stopCluster(cl)


proj <- projectNewData(sc_pred_ss, as.matrix(pred_data@raw.data))
proj_ss <- spatialSign(proj)
sc_pred_ss@projection <- proj_ss


sc_pred_ss <- scPredict(sc_pred_ss, useProj = TRUE)
predictions_ss <- getPredictions(sc_pred_ss)

predictions_ss$true <- factor(pred_data@meta.data$cellType1)
predictions_ss <- predictions_ss %>% 
  mutate(predClass = if_else(Dendritic_cell > 0.9, "DC", 
                             if_else(Dendritic_cell < 0.1, "Monocyte", "Unassigned")))



tidy_table(predictions_ss, "predClass", "true") %>% 
  select(-B_cell, -Blood_progenitor, -T_cell, -Natural_killer)

tidy_table(predictions_ss, "predClass", "true", prop = TRUE) %>% 
  select(-B_cell, -Blood_progenitor, -T_cell, -Natural_killer)


# scPred version: d7661592

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
