# Load libraries ----------------------------------------------------------
library("tidyverse")
library("scPred")
library("here")
library("cowplot")

# Read datasets -----------------------------------------------------------

# Read training
reference <- readRDS(here("results", "2018-04-16_pancreas_reference_processing", "pancreas_processed_cpm.RDS"))
reference_metadata <- readRDS(here("results", "2018-04-16_pancreas_reference_processing", "pancreas_processed_metadata.RDS"))

# Read test
baron_cpm <- readRDS(here("results", "2018-04-15_baron_processing", "baron_processed_cpm.RDS"))
baron_metadata <- readRDS(here("results", "2018-04-15_baron_processing", "baron_processed_metadata.RDS"))


output <- file.path("results", "2018-04-16_baron_prediction")

# Eigendecompose training gene expression data ----------------------------
pancreas <- eigenDecompose(t(reference))
scPred::metadata(pancreas) <- reference_metadata

# Get informative principal components ------------------------------------
pancreas <- getInformativePCs(pancreas, pVar = "cellType", sig = 0.05)

# Train model -------------------------------------------------------------
pancreas <- trainModel(object = pancreas, seed = 66, method = "svmRadial")
saveRDS(pancreas, here(output, "object_svmRadial.RDS"))

# Classify cells in new dataset -------------------------------------------
baron_predictions <- eigenPredict(object = pancreas, newData = t(baron_cpm))
saveRDS(baron_predictions, here(output, "predictions_svmRadial.RDS"))


table(baron_predictions$class)

baron_predictions %>% 
  filter(class %in% c("alpha", "beta", "delta", "gamma")) %>% 
  group_by(class) %>% 
  summarise(n = n()) -> pred

baron_metadata %>% 
  as.data.frame() %>% 
  filter(cell_type1 %in% c("alpha", "beta", "delta", "gamma")) %>% 
  group_by(cell_type1) %>% 
  summarise(n = n()) -> true

pred$n/true$n
