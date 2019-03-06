# Script information ------------------------------------------------------

# title:
# author: José Alquicira Hernández
# date:
# description:


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary
library("scmap")
library("SingleCellExperiment")


# Set output --------------------------------------------------------------


output_dir_name <- "scmap_comparison" # <------ Output directory

date <- "2019-03-04" # <------ Date

output <- file.path("results", paste(date, output_dir_name, sep = "_"))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input    <- file.path("results", "2018-05-12_pancreas_processing") # <------ Input directory

# Read training
training <- readRDS(here(input, "baron_cpm_train.RDS"))
training_metadata <- readRDS(here(input, "baron_metadata_train.RDS"))

# Read testing
test <- readRDS(here(input, "pancreas_cpm_test.RDS"))
test_metadata <- readRDS(here(input, "pancreas_metadata_test.RDS"))


# Pre-process datasets ----------------------------------------------------

sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(training)), 
                            colData = training_metadata)
logcounts(sce) <- log2(normcounts(sce) + 1)
rowData(sce)$feature_symbol <- rownames(sce)

getSCE <- function(x, meta){
  x <- SingleCellExperiment(assays = list(normcounts = as.matrix(x)), 
                            colData = meta)
  logcounts(x) <- log2(normcounts(x) + 1)
  rowData(x)$feature_symbol <- rownames(x)
  x
}

test <- mapply(test, test_metadata, FUN =  getSCE)

# Feature selection -------------------------------------------------------

sce <- selectFeatures(sce, suppress_plot = FALSE)
table(rowData(sce)$scmap_features)


# scmap cluster -----------------------------------------------------------

sce <- indexCluster(sce, cluster_col = "x.cell_type.i.")
head(metadata(sce)$scmap_cluster_index)
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))



# scmap cluster projection ------------------------------------------------

datasets <- c("muraro", "segerstolpe", "xin")

contTable <- function(x, y){
  
  islets <- c("alpha", "beta", "delta", "gamma")
  
  
  table(x, y) %>% 
    as.data.frame() %>% 
    spread(key = "y", value = "Freq", fill = 0) %>%
    as.data.frame() %>% 
    column_to_rownames("x") -> counts
  
  i <- colnames(counts) %in% islets
  orderClass <- c(islets, colnames(counts)[!i])
  counts <- counts[,orderClass]
  

  labels <- rownames(counts)
  
  counts %>% 
    mapply(function(x,d){x/d}, ., colSums(.)) %>% 
    `rownames<-`(labels) %>% 
    round(2) -> props

  list(counts = counts, props = props)
}

writeResults<- function(x, method = "cluster", dataset){
  write.table(x$counts, file = here(output, paste0(dataset, "_", method, ".txt")), 
              quote = FALSE, sep = "\t")
  cat("\n", file = here(output, paste0(dataset, "_", method, ".txt")), append = TRUE)
  write.table(x$props, file = here(output,paste0(dataset, "_", method, ".txt")), 
              quote = FALSE, sep = "\t", append = TRUE)
}

scmap_cluster <- function(dataset){
  scmapCluster <- scmapCluster(
    projection = test[[dataset]], 
    index_list = list(
      baron = metadata(sce)$scmap_cluster_index
    )
  )
  scmapCluster
}

getClusterResults <- function(dataset){
  
  scmap_cluster(dataset) %>% 
    .$combined_labs %>% 
    contTable(test[[dataset]]$cell_type1) %>% 
    writeResults(dataset = dataset)
  
}


lapply(datasets, getClusterResults)


# scmap cell --------------------------------------------------------------

scmap_cell <- function(dataset){
  set.seed(1)
  sce <- indexCell(sce)
  scmapCellRes <- scmapCell(
    projection = test[[dataset]], 
    list(
      baron = metadata(sce)$scmap_cell_index
    )
  )
  
  scmapCell2Cluster(
    scmapCellRes, 
    list(
      as.character(colData(sce)$x.cell_type.i.)
    )
  )
}


getCellResults <- function(dataset){
  scmap_cell(dataset) %>% 
    .$combined_labs %>% 
    contTable(y = test[[dataset]]$cell_type1) %>% 
    writeResults(method = "cell", dataset = dataset)
}

lapply(datasets, getCellResults)

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))