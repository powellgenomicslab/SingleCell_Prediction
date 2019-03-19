# Script information ------------------------------------------------------

# title: Predict PBMCs using tree approach
# author: José Alquicira Hernández
# date: 2019-03-20
# description: None


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary
library("scPred")
library("Seurat")
library("BiocParallel")

# Set output --------------------------------------------------------------


output_dir_name <- "pbmc_prediction" # <------ Output directory

date <- "2019-03-20" # <------ Date

output <- file.path("results", paste(date, output_dir_name, sep = "_"))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Read prediction models
input  <- file.path("results", "2019-03-15_trained_models") # <------ Input directory
readData <- function(i) readRDS(here(input, paste0("models_layer_", i , ".RDS")))
scTree <- lapply(1:3, readData)
names(scTree) <- paste0("layer", 1:3)

# Read gene expression data
pbmc <- readRDS(here("results", "2019-03-13_pbmc_assign_layers", "pbmc.RDS"))


# Create tree prediction function -----------------------------------------

predictTree <- function(br, newData){
  
  # Layer 1 .........................................................
  
  ## Perform predictions
  predPBMC <- scPredict(object = scTree$layer1[[br]], newData = newData) %>% getPredictions()
  
  ## Create data.frame with columns as prediction layers
  preds <- split(predPBMC, predPBMC$predClass)
  cells <- predPBMC %>% 
    select(predClass) %>% 
    rename(layer1 = predClass) %>% 
    rownames_to_column("ID")
  
  ## Fill empty layers
  cells$layer2 <- NA
  cells$layer3 <- NA
  
  ## If there are cells classified as non-lymphoid cells, fill layers 
  ## 2 and 3 with predictions from layer 1
  
  noLymphoid <- names(preds) != "Lymphoid_cell"
  if(any(noLymphoid)){
    
    preds[noLymphoid] %>% 
      lapply("[", "predClass") %>% 
      reduce(rbind) %>% 
      rownames_to_column("ID") -> noLymphoid
    
    i <- match(noLymphoid$ID, cells$ID)
    cells$layer2[i] <- noLymphoid$predClass
    cells$layer3[i] <- noLymphoid$predClass
    
    
  }
  
  
  # Layer 2 .........................................................
  
  ## If there are cells classified as Lymphoid, then predictions for
  ## layer 2 are performed
  
  if(!is.null(preds$Lymphoid)){
    
    ## Extract cells predicted as Lymphoid
    i <- colnames(newData) %in% rownames(preds$Lymphoid)
    predLymphoid <- newData[,i]
    
    ## Perform predictions
    predLymphoid <- scPredict(object = scTree$layer2[[br]], newData = predLymphoid) %>% getPredictions()
    
    predLymphoid %>% 
      rownames_to_column("ID") -> predLymphoid
    
    ## Fill layer 2 with predictions
    i <- match(predLymphoid$ID, cells$ID)
    cells$layer2[i] <- predLymphoid$predClass
    
    
    ## Evaluate if there are any T cells
    preds2 <- split(predLymphoid, predLymphoid$predClass)
    noTCell <- names(preds2) != "T_cell"
    
    ## If there are cells classified as T cells, fill layers 
    ## 2 and 3 with predictions from layer 2
    
    if(any(noTCell)){
      
      preds2[noTCell] %>% 
        lapply("[", c("ID", "predClass")) %>% 
        reduce(rbind)  -> noTCell
      
      i <- match(noTCell$ID, cells$ID)
      cells$layer2[i] <- noTCell$predClass
      cells$layer3[i] <- noTCell$predClass
    }
    
    
    # Layer 3 .........................................................
    
    ## If there are cells classified as T cell, predictions for layer 3 
    ## are performed
    
    if(!is.null(preds2$T_cell)){
      
      ## Extract cells predicted as T cells
      i <- colnames(newData) %in% preds2$T_cell$ID
      predTCell <- newData[,i]
      predTCell <- scPredict(object = scTree$layer3[[br]], newData = predTCell) %>% getPredictions()
      
      predTCell %>% 
        rownames_to_column("ID") %>% 
        mutate(predClass = if_else(cytotoxic < 0.25, 
                                   "non_cytotoxic", 
                                   if_else(cytotoxic > 0.75, 
                                           "cytotoxic", 
                                           "unassigned"))) -> predTCell
      
      i <- match(predTCell$ID, cells$ID)
      cells$layer3[i] <- predTCell$predClass
      
    }
    
    
  }
  
  
  cells
  
}

# Perform predictions -----------------------------------------------------


# Set bootstrap seeds
set.seed(66)
seed_part <- sample(seq_len(10e4), 10)

# Get test folds
createTestDatasets <- function(seed){
  set.seed(seed)
  i <- createDataPartition(seq_len(nrow(pbmc@meta.data)), times = 1, p = 0.5, list = FALSE)
  pbmc@meta.data %>% 
    row.names() %>% 
    `[`(-i) -> testCells
  
  testData <- SubsetData(pbmc, cells.use = testCells, do.clean = TRUE)
  testData
}


multicoreParam <- MulticoreParam(workers = 5)
testFolds <- bplapply(seed_part, createTestDatasets, BPPARAM = multicoreParam)

testFolds %>% 
  bplapply(slot, "data", BPPARAM = multicoreParam) %>% 
  bplapply(as.matrix, BPPARAM = multicoreParam) -> testData

names(testFolds) <- paste0("r", seed_part)

applyTree <- function(br, newData){
  predictTree(br = br, newData = newData)
}


res <- bpmapply(applyTree, paste0("r", seed_part), testData, SIMPLIFY = FALSE, BPPARAM = multicoreParam)
saveRDS(res, here(output, "predictions.RDS"))


# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))