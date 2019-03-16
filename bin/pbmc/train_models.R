# Script information ------------------------------------------------------

# title: Model training
# author: José Alquicira Hernández
# date:2019-03-15
# description: None


# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("here")

# Secondary
library("scPred")
source(here("bin", "quickSeurat.R"))


# Set output --------------------------------------------------------------

output_dir_name <- "trained_models" # <------ Output directory
date <- "2019-03-15" # <------ Date
output <- file.path("results", paste(date, output_dir_name, sep = "_"))
if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input
input    <- file.path("results", "2018-03-13_pbmc_assign_layers") # <------ Input directory
filename <- "pbmc.RDS" # <------ Input file


# Read file
data <- readRDS(file = here(input, filename))

# Set bootstrap seeds
set.seed(66)
seed_part <- sample(seq_len(10e4), 10)



# Train layers ------------------------------------------------------------

trainLayer <- function(seed, level){
  
  extractLevel <- paste0("level", level - 1)
  layer <- paste0("level", level)
  
  # Create train and test folds
  set.seed(seed)
  i <- createDataPartition(seq_len(nrow(data@meta.data)), times = 1, p = 0.5, list = FALSE)
  
  if(level - 1){
    if(level == 2){
      subtype <- "Lymphoid_cell"
    }else if(level == 3){
      subtype <- "T_cell"
    }
    data@meta.data[i,] %>% 
      rownames_to_column("id") %>% 
      filter(!!sym(extractLevel) == subtype) %>% 
      pull(id) -> trainCells
  }else{
    data@meta.data %>% 
      row.names() %>% 
      `[`(i) -> trainCells
  }
  
  trainData <- SubsetData(data, cells.use = trainCells, do.clean = TRUE)

  # Process data
  trainData <- processSeurat(trainData)
  
  # Perform PCA
  trainData <- RunPCA(object = trainData, 
                      pc.genes = trainData@var.genes, 
                      do.print = FALSE)
  
  # Get feature space
  trainData@meta.data[[level]] <- as.factor(trainData@meta.data[[level]])
  scpred <- getFeatureSpace(trainData, pVar = layer)
  
  # Train model
  scpred <- trainModel(scpred)
  
  scpred
  
}


res1 <- lapply(seed_part, trainLayer, 1)
res2 <- lapply(seed_part, trainLayer, 2)
res3 <- lapply(seed_part, trainLayer, 3)

res <- list(level1 = res1, level2 = res2, level3 = res3) %>% 
  lapply(function(x){names(x) <- paste0("r", seed_part); x})


# Save results ------------------------------------------------------------

saveRDS(res, file = here(output, "models.RDS"))

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))