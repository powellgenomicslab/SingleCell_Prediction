# Script information ------------------------------------------------------

# title: Predict cell types for test fold of pure PBMCs
# author: José Alquicira Hernández
# date: 2018-11-13
# description:
#   scPred version : 31e6f55653d596af92eaafa0c4edab477776f58a


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary

library("scPred")
library("Seurat")

# Set output --------------------------------------------------------------


output_dir_name <- "10x_immune_tree_prediction" # <------ Output directory

date <- "2018-11-13" # <------ Date

output <- file.path("results", paste(date, output_dir_name, sep = "_"))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Read first layer model


input_dir_name <- "pbmc_myeloid_lymphoid_prediction" # <------ Output directory
date <- "2018-11-12" # <------ Date
input <- file.path("results", paste(date, input_dir_name, sep = "_"))
layer1 <- readRDS(here(input, "scPred.RDS"))



# Read second layer

input_dir_name <- "pbmc_lymphoid_prediction" # <------ Output directory

date <- "2018-11-09" # <------ Date
input <- file.path("results", paste(date, input_dir_name, sep = "_"))
layer2 <- readRDS(here(input, "scpred.RDS"))



# Read third layer model

input_dir_name <- "pbmc_t-cell_prediction" # <------ Output directory
date <- "2018-11-09" # <------ Date
input <- file.path("results", paste(date, input_dir_name, sep = "_"))
layer3 <- readRDS(here(input, "scpred_cytotoxic_vs_no-cytotoxic.RDS"))



# Read test data

input <- file.path("results", "2018-11-08_pbmc_pure_processed")

input %>% 
  here() %>%
  list.files(full.names = FALSE) %>%
  str_subset(regex("_test.RDS$")) -> filenames

readData <- function(filename) {
  readRDS(here(input, filename))
}


testData <- lapply(filenames, readData)
names(testData) <- str_remove(filenames, "_test.RDS")

# Merge datasets

combined <- testData[[1]]
for (i in 2:length(x = testData)) {
  combined <- MergeSeurat(object1 = combined, object2 = testData[[i]], add.cell.id1 = "_2", "_3")
}


# Create tree object ------------------------------------------------------

setClass("scPred-tree", representation(models = "list"), prototype(models = list()))


createPredTree <- function(...){
  models <- list(...)
  for(i in seq_len(length(models))){
    if(class(models[[i]]) != "scPred"){
      stop("Object ", i, " is not an scPred object")
    }
  }
  models
}


tree <- createPredTree(layer1 = layer1, layer2 = layer2, layer3 = layer3)



# Create prediction function ----------------------------------------------

predictTree <- function(newData){
  
  # Layer 1
  predPBMC <- scPredict(object = tree$layer1, newData = newData) %>% getPredictions()
  preds <- split(predPBMC, predPBMC$predClass)
  
  cells <- predPBMC %>% 
    select(predClass) %>% 
    rename(layer1 = predClass) %>% 
    rownames_to_column("ID")
  
  cells$layer2 <- NA
  cells$layer3 <- NA
  
  # Layer 2
  
  noLymphoid <- names(preds) != "Lymphoid"
  if(any(noLymphoid)){
    
    preds[noLymphoid] %>% 
      lapply("[", "predClass") %>% 
      reduce(rbind) %>% 
      rownames_to_column("ID") -> noLymphoid
    
    i <- match(noLymphoid$ID, cells$ID)
    cells$layer2[i] <- noLymphoid$predClass
    cells$layer3[i] <- noLymphoid$predClass
    
    
  }
  
  if(!is.null(preds$Lymphoid)){
    
    predLymphoid <- SubsetData(newData, cells.use = rownames(preds$Lymphoid), do.clean = TRUE)
    predLymphoid <- scPredict(object = tree$layer2, newData = predLymphoid) %>% getPredictions()
    
    predLymphoid %>% 
      rownames_to_column("ID") -> predLymphoid
    
    i <- match(predLymphoid$ID, cells$ID)
    cells$layer2[i] <- predLymphoid$predClass
    
    preds2 <- split(predLymphoid, predLymphoid$predClass)
    
    
    noTCell <- names(preds2) != "T_cell"
    
    if(any(noTCell)){
      
      preds2[noTCell] %>% 
        lapply("[", c("ID", "predClass")) %>% 
        reduce(rbind)  -> noTCell
      
      i <- match(noTCell$ID, cells$ID)
      cells$layer2[i] <- noTCell$predClass
      cells$layer3[i] <- noTCell$predClass
    }
    
    
    
    
    # Layer 3
    if(!is.null(preds2$T_cell)){
      
      predTCell <- SubsetData(newData, cells.use = preds2$T_cell$ID, do.clean = TRUE)
      predTCell <- scPredict(object = tree$layer3, newData = predTCell) %>% getPredictions()
      
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



# Predict cell types ------------------------------------------------------

res <- predictTree(newData = combined)
res$true <- combined@meta.data$cellType

# Relabel T cells
i <- which(combined@meta.data$layer2 == "T_cell")
j <- combined@meta.data$cellType %in% c("cytotoxic_t", "naive_cytotoxic")

k <- intersect(i,which(j))
res$true[k] <- "cytotoxic"

k <- intersect(i,which(!j))
res$true[k] <- "non_cytotoxic"

# Get accuracy ------------------------------------------------------------


crossTab <- function(tab, true, pred, fill = 0, prop = TRUE, digits = 2){
  

  
  tab %>% 
    group_by_(pred, true) %>% 
    summarise(n = n()) %>% 
    spread(key = true, value = "n", fill = fill) %>% 
    as.data.frame() %>% 
    column_to_rownames(pred) -> x
  
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

cellNames <- c("b_cells", "Myeloid_cell", "Progenitor_cell", 
               "cd56_nk", "cytotoxic", "non_cytotoxic", "unassigned")

crossTab(res, true = "true", pred = "layer3")[cellNames, ] %>% 
  knitr::kable() %>% 
  writeLines(here(output, "crosstab_prop.txt"))

crossTab(res, true = "true", pred = "layer3", prop = FALSE)[cellNames, ] %>% 
  knitr::kable() %>% 
  writeLines(here(output, "crosstab_counts.txt"))




# Save predictions --------------------------------------------------------

saveRDS(res, here(output, "predictions.RDS"))

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
