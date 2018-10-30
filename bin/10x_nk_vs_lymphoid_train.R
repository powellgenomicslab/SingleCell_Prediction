# Script information ------------------------------------------------------
# title: Classification of Nk cells vs T-cells and B-cells
# author: José Alquicira Hernández
# date: 2018-10-30
# description: 
#   - scPred version (d7661592)


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary
library("Seurat")
library("scPred")


# Set output --------------------------------------------------------------


output_dir_name <- "10x_nk-cells_vs_t-cells_b-cells_" # <------ Output directory
date <- "2018-10-30_"
output <- file.path("results", paste0(date, output_dir_name))

if (!dir.exists(output)) {
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input <- file.path("results", "2018-10-12_10x_immune_tree") # <------ Input directory
filename <- "train_data.RDS" # <------ Input file
train_data <- readRDS(here(input, filename))

# Subset cells

## Extract cell indexes
i <- rownames(train_data@meta.data)[train_data@meta.data$cellType2 == "Lymphoid_cell"]

## Subset data
train_data <- SubsetData(train_data, cells.use = i, do.clean = TRUE)

## Refactorize metadata

newClasses <- ifelse(train_data@meta.data$cellType1 == "Natural_killer", "Natural_killer", "other")
newClasses <- factor(newClasses, levels = c("Natural_killer", "other"))

train_data@meta.data$cellType2 <- factor(train_data@meta.data$cellType2)
train_data@meta.data$cellType1 <- newClasses


# Process data ------------------------------------------------------------

cat("\nProcessing training data (normalization, scaling, and PCA)...\n")

train_data <- train_data %>%
  NormalizeData(
    normalization.method = "LogNormalize",
    scale.factor = 10000
  ) %>%
  FindVariableGenes(
    mean.function = ExpMean,
    dispersion.function = LogVMR,
    x.low.cutoff = 0.02,
    x.high.cutoff = 5,
    y.cutoff = 0.5
  ) %>%
  ScaleData() %>%
  RunPCA(
    pc.genes = .@var.genes,
    do.print = TRUE,
    pcs.print = 1:5,
    genes.print = 5
  )


# Get feature space -------------------------------------------------------
cat("\nGetting feature space...\n")

sc_pred <- getFeatureSpace(train_data, pVar = "cellType1")
rm(train_data)


# Train models ------------------------------------------------------------
cat("\nTraining models...\n")

library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
sc_pred <- trainModel(sc_pred)
stopCluster(cl)

cat("\nSaving scPred object...\n")

saveRDS(sc_pred, file = here(output, "scpred_train.RDS"))

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
