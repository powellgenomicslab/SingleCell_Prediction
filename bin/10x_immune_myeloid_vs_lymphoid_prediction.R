# Script information ------------------------------------------------------

# title: Classification of precursor, myeloid, and lymphoid cells
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


output_dir_name <- "10x_immune_myeloid_vs_lymphoid_prediction" # <------ Output directory

date <- format(Sys.Date(), format = "%Y-%m-%d_")
date <- "2018-10-12_" # <------ Fix date

output <- file.path("results", paste0(date, output_dir_name))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input <- file.path("results", "2018-10-12_10x_immune_tree") # <------ Input directory

# Read files

filename <- "train_data.RDS" # <------ Input file
train_data <- readRDS(here(input, filename))

filename <- "pred_data.RDS" # <------ Input file
pred_data <- readRDS(here(input, filename))



# Perform PCA -------------------------------------------------------------

train_data <- train_data %>% 
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>% 
  # FindVariableGenes( mean.function = ExpMean, 
  #                   dispersion.function = LogVMR, 
  #                   x.low.cutoff = 0.0125, 
  #                   x.high.cutoff = 3, 
  #                   y.cutoff = 0.5) %>% 
  ScaleData() %>% 
  RunPCA(pc.genes = rownames(.@data),
         do.print = TRUE, 
         pcs.print = 1:5, 
         genes.print = 5)


# plotPCA(train_data, group = "cellType1")
# plotPCA(train_data, 1, 3,group = "cellType1")
# plotPCA(train_data, 2, 3,group = "cellType1")



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


# Get performance ---------------------------------------------------------

sc_pred@predMeta <-  pred_data@meta.data
ct <- crossTab(sc_pred, true = "cellType1")


writeLines(knitr::kable(ct), 
            here(output, "cont_table.txt"))

saveRDS(predictions, file = here(output, "predictions.RDS"))
saveRDS(sc_pred, file = here(output, "scpred.RDS"))


p <- plotPredProbs(sc_pred, facet = "cellType1")

ggsave(p, filename = here(output, "prob_dist.png"), width = 9 , height = 2)


# scPred version: d7661592

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
