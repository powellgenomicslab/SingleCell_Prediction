# Script information ------------------------------------------------------

# title: Lymphocyte training
# author: José Alquicira Hernández
# date: 2018-11-09
# description:
#   scPred version : 31e6f55653d596af92eaafa0c4edab477776f58a


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary
library("Seurat")
library("scPred")
source(here("bin", "quickSeurat.R"))

# Set output --------------------------------------------------------------


output_dir_name <- "pbmc_lymphoid_prediction" # <------ Output directory

date <- "2018-11-09" # <------ Date
output <- file.path("results", paste(date, output_dir_name, sep = "_"))


if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input <- file.path("results", "2018-11-08_pbmc_pure_processed")

input %>% 
  here() %>%
  list.files(full.names = FALSE) %>%
  str_subset(regex("_train.RDS$")) %>%
  str_subset(regex("^(?!(cd14_monocytes|cd34))")) -> filenames


# Read file

readData <- function(filename) {
  readRDS(here(input, filename))
}

data <- lapply(filenames, readData)
names(data) <- str_remove(filenames, "_train.RDS")



# Merge datasets ----------------------------------------------------------

combined <- data[[1]]
for (i in 2:length(x = data)) {
  combined <- MergeSeurat(object1 = combined, object2 = data[[i]])
}

# Normalize data ----------------------------------------------------------

data <- processSeurat(combined)


# Perform PCA -------------------------------------------------------------

data <- RunPCA(object = data, 
               pc.genes = data@var.genes, 
               do.print = TRUE, 
               pcs.print = 1:5, 
               genes.print = 5)


# Get feature space -------------------------------------------------------

sc_pred <- getFeatureSpace(data, pVar = "layer2")

# Train models ------------------------------------------------------------


library(doParallel)
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
sc_pred <- trainModel(sc_pred)
stopCluster(cl)




# Perform predictions -----------------------------------------------------

# Read test data

input %>% 
  here() %>%
  list.files(full.names = FALSE) %>%
  str_subset(regex("_test.RDS$")) %>%
  str_subset(regex("^(?!(cd14_monocytes|cd34))")) -> filenames

testData <- lapply(filenames, readData)
names(testData) <- str_remove(filenames, "_test.RDS")



# Merge datasets

combined <- testData[[1]]
for (i in 2:length(x = testData)) {
  combined <- MergeSeurat(object1 = combined, object2 = testData[[i]])
}



sc_pred <- scPredict(sc_pred, combined)
sc_pred@predMeta <- combined@meta.data



# Get accuracy ------------------------------------------------------------

crossTab(sc_pred, "layer2") %>% 
  knitr::kable() %>% 
  writeLines(here(output, "crosstab_prop.txt"))

crossTab(sc_pred, "layer2", prop = FALSE) %>% 
  knitr::kable() %>% 
  writeLines(here(output, "crosstab_counts.txt"))
  


# Save scPred object ------------------------------------------------------

saveRDS(sc_pred, file = here(output, "scpred.RDS"))


# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
