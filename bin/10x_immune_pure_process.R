# Script information ------------------------------------------------------

# title: Process pure PBMC data
# author: José Alquicira Hernández
# date: 2018-11-08
# description: Processes PBMC datasets and splits the data into train and
# test for each cell type


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary
library("Seurat")

# Set output --------------------------------------------------------------


output_dir_name <- "pbmc_pure_processed" # <------ Output directory

date <- "2018-11-08" # <------ Date

output <- file.path("results", paste(date, output_dir_name, sep = "_"))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input
input <- file.path("data", "2018-10-02_10x_immune_cells") # <------ Input directory


readData <- function(id){
  # Set filename
  filename <- paste(id, "filtered_matrices_mex", sep = "_") # <------ Input file
  # Read file
  x <- Read10X(data.dir = here(input, filename, "hg19"))
  genes <- read.table(here(input, filename, "hg19", "genes.tsv")) %>% 
    reduce(paste, sep = "_")
  rownames(x) <- genes
  colnames(x) <- paste(colnames(x), id, sep = "-")
  x
}

# Get cell type ids
ids <- list.dirs(here(input), full.names = FALSE) %>% 
  str_subset("mex") %>% 
  str_subset(regex("[^hg19]$")) %>% 
  str_remove("_filtered_matrices_mex")

# Read gene expression matrices
data <- lapply(ids, readData)
names(data) <- ids



# Create Seurat objects ---------------------------------------------------
data <- mapply(data, ids, 
               FUN =  function(d, i) {
                 CreateSeuratObject(d, meta.data = data.frame(cellType = rep(i, ncol(d)), 
                                                              row.names = colnames(d)))
               })




# Add prediction layers ---------------------------------------------------

# processSeurat <- function(object){
#   object %>% 
#     NormalizeData(normalization.method = "LogNormalize", 
#                   scale.factor = 10000) %>% 
#     FindVariableGenes(mean.function = ExpMean, 
#                       dispersion.function = LogVMR, 
#                       x.low.cutoff = 0.0125, 
#                       x.high.cutoff = 3, 
#                       y.cutoff = 0.5, 
#                       do.plot = FALSE) %>% 
#     ScaleData()
# }
# 
# 
# data <- lapply(data, processSeurat)

selectLayer1 <- function(x){
  if(x == "cd14_monocytes"){
    "Myeloid_cell"
  }else if(x == "cd34"){
    "Progenitor_cell"
  }else{
    "Lymphoid"
  }
}


selectLayer2 <- function(x){
    switch(x,
           "cd4_t_helper" = "T_cell",
           "regulatory_t" = "T_cell",
           "naive_t" = "T_cell",
           "memory_t" = "T_cell",
           "cytotoxic_t" = "T_cell",
           "naive_cytotoxic" = "T_cell", x)
}




assignLayers <- function(object, name) {
  
  layer1 <- selectLayer1(name)
  object@meta.data$layer1 <- layer1
  
  if(layer1 == "Lymphoid"){
    layer2 <- selectLayer2(name)
    object@meta.data$layer2 <- layer2
  }
  
  object
  
}



data <- mapply(data, names(data), FUN = assignLayers)



splitData <- function(cellType, seed, p = 0.5){
  
  i <- rownames(cellType@meta.data)
  n <- length(i)
  set.seed(seed)
  j <- sample(seq_len(n), size = n*p)
  iTrain <- i[j]
  iTest <-  i[-j]
  
  cellTypeTrain <- SubsetData(cellType, cells.use = iTrain, do.clean = TRUE)
  cellTypeTest  <- SubsetData(cellType, cells.use = iTest, do.clean = TRUE)
  
  list(train = cellTypeTrain, test = cellTypeTest)
  
}


# Assign prediction layer labels
data <- lapply(data, splitData, seed = 66)

# Save data ---------------------------------------------------------------

saveCellType <- function(cellType, name) {
  
  saveRDS(cellType$train,
    file = here(
      output,
      paste0(name, "_train.RDS")
    )
  )
  
  saveRDS(cellType$test,
    file = here(
      output,
      paste0(name, "_test.RDS")
    )
  )
}


mapply(data, names(data), FUN = saveCellType)


# Create software test dataset --------------------------------------------

n <- lapply(data, function(x) nrow(x$train@meta.data) + nrow(x$test@meta.data)) %>% reduce(sum)
cat("\nTotal number of cells:", n, "\n")


createTest <- function(cellTypeSet, seeds = c(66, 99)){
  
  set.seed(seeds[1])
  i <- sample(seq_len(nrow(cellTypeSet$train@meta.data)), size = 100)
  cellsTrain <- rownames(cellTypeSet$train@meta.data)[i]
  
  
  set.seed(seeds[2])
  j <- sample(seq_len(nrow(cellTypeSet$test@meta.data)), size = 100)
  cellsTest <- rownames(cellTypeSet$test@meta.data)[j]
  
  train <- SubsetData(cellTypeSet$train, cells.use = i, do.clean = TRUE)
  test  <- SubsetData(cellTypeSet$test, cells.use = j, do.clean = TRUE)
  
  list(train = train, test = test)
}


subData <- lapply(data, createTest)

saveCellType <- function(cellType, name) {
  
  saveRDS(cellType$train,
          file = here(
            output,
            paste0(name, "_train_sub.RDS")
          )
  )
  
  saveRDS(cellType$test,
          file = here(
            output,
            paste0(name, "_test_sub.RDS")
          )
  )
}

mapply(subData, names(subData), FUN = saveCellType)


# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
