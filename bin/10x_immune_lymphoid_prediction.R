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


output_dir_name <- "10x_t-cells_b-cells_nk-cells" # <------ Output directory

date <- format(Sys.Date(), format = "%Y-%m-%d_")

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


filename <- "train_data.RDS" # <------ Input file
train_data <- readRDS(here(input, filename))


filename <- "pred_data.RDS" # <------ Input file
pred_data <- readRDS(here(input, filename))

# Subset cells

pred_data@meta.data$pred1 <- as.factor(predictions$predClass)

# train_data <- RunTSNE(train_data)
# TSNEPlot(train_data, group = "cellType1")


i <- rownames(train_data@meta.data)[train_data@meta.data$cellType2 %in% c("Lymphoid_cell")]
j <- rownames(pred_data@meta.data)[pred_data@meta.data$pred1 %in% c("Lymphoid_cell")]


train_data <- SubsetData(train_data, cells.use = i, do.clean = TRUE)
train_data@meta.data$cellType2 <- factor(train_data@meta.data$cellType2)

pred_data <- SubsetData(pred_data, cells.use = j, do.clean = TRUE)
pred_data@meta.data$pred1 <- factor(pred_data@meta.data$pred1)


train_data@meta.data$cellType1 <- factor(train_data@meta.data$cellType1)

gc()

# Process data ------------------------------------------------------------

train_data <- train_data %>% 
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>% 
  FindVariableGenes(mean.function = ExpMean, 
                    dispersion.function = LogVMR, 
                    x.low.cutoff = 0.02, 
                    x.high.cutoff = 5, 
                    y.cutoff = 0.5) %>% 
  ScaleData() %>% 
  RunPCA(pc.genes = .@var.genes, 
         do.print = TRUE, 
         pcs.print = 1:5, 
         genes.print = 5)

PCAPlot(train_data, group = "cellType1") -> q


# Get feature space -------------------------------------------------------

sc_pred <- getFeatureSpace(train_data, pVar = "cellType1")
rm(train_data)


# Train models ------------------------------------------------------------

library(doParallel)
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
sc_pred <- trainModel(sc_pred)
stopCluster(cl)

predictions$true <- pred_data@meta.data$cellType2


saveRDS(predictions, file = here(output, "predictions.RDS"))
saveRDS(sc_pred, file = here(output, "scpred.RDS"))


# scPred version: d7661592

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
