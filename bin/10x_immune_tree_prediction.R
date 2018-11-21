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


readModels <- function(input) {
  readRDS(here("results", input, "scpred.RDS"))
}

dates <- c("2018-11-12", "2018-11-09", "2018-11-09")
dirs <- c("pbmc_myeloid_lymphoid_prediction", "pbmc_lymphoid_prediction", "pbmc_t-cell_prediction")
dirs <- paste(dates, dirs, sep = "_")

tree <- lapply(dirs, readModels)
names(tree) <- paste0("layer", 1:3)


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
  combined <- MergeSeurat(object1 = combined, object2 = testData[[i]])
}



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



# Assign true layers
assignLayer1 <- function(x) switch(x,
                                   "b_cells" = "Lymphoid",
                                   "cd14_monocytes" = "Myeloid_cell",
                                   "cd34" = "Progenitor_cell",
                                   "cd56_nk" = "Lymphoid",
                                   "cytotoxic" = "Lymphoid",
                                   "non_cytotoxic" = "Lymphoid")

assignLayer2 <- function(x) switch(x,
                                   "cd14_monocytes" = "Myeloid_cell",
                                   "cd34" = "Progenitor_cell",
                                   "cytotoxic" = "T_cell",
                                   "non_cytotoxic" = "T_cell", 
                                   x)

assignLayer3 <- function(x) switch(x,
                                   "cd14_monocytes" = "Myeloid_cell",
                                   "cd34" = "Progenitor_cell", 
                                   x)


list(f1 = assignLayer1, f2 = assignLayer2, f3 = assignLayer3) %>% 
  lapply(Vectorize) %>% 
  lapply(function(x) x(res$true)) %>% 
  as.data.frame() %>% 
  set_names(paste0("true", 1:3)) %>% 
  cbind(res, .) -> res


# Assign order for contingency tables
order1 <- c("Myeloid_cell", "Lymphoid", "Progenitor_cell", "unassigned")
order2 <- c("b_cells", "T_cell", "cd56_nk", "Myeloid_cell", "Progenitor_cell", "unassigned")
order3 <- c("cytotoxic", "non_cytotoxic",  "b_cells", "cd56_nk", "Progenitor_cell", "unassigned")


order <- list(order1, order2, order3) %>% set_names(1:3)



getAccuracy <- function(i, level, prop){
  i <- i[[level]]
  true <- paste0("true", level)
  layer <- paste0("layer", level)
  crossTab(res, true = true, pred = layer, prop = prop) %>%
    .[i, i[-length(i)]] %>% 
    knitr::kable() 
}



# Save accuracy -----------------------------------------------------------


saveAccuracy <- function(level){
  level <- level %>% as.integer()
  getAccuracy(order, level, FALSE) %>% 
    writeLines(here(output, paste0("crosstab_counts_", level,".txt")))
  
  getAccuracy(order, level, TRUE) %>% 
    writeLines(here(output, paste0("crosstab_props_", level,".txt")))
  
}


lapply(names(order), saveAccuracy)


# Save predictions --------------------------------------------------------

saveRDS(res, here(output, "predictions.RDS"))

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
