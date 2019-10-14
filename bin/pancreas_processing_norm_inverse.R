library("tidyverse")
library("here")


# Read data ---------------------------------------------------------------

## Baron
baron <- readRDS("data/2018-04-15_pancreas_baron/baron-human.rds")
baron_counts <- counts(baron)
baron_metadata <- as.data.frame(colData(baron))
rm(baron)

## Muraro
muraro <- readRDS("data/2018-04-15_pancreas_muraro/muraro.rds")
muraro_counts <- normcounts(muraro)
muraro_metadata <- as.data.frame(colData(muraro))
rm(muraro)


# Normalize gene names
muraro_counts %>% 
  rownames() %>% 
  str_split("__") %>% 
  lapply("[", 1) %>% 
  unlist() -> muraro_genenames

rownames(muraro_counts) <- muraro_genenames
rm(muraro_genenames)

## Segerstolpe

segerstolpe <- readRDS("data/2018-04-15_pancreas_segerstolpe/segerstolpe.rds")
segerstolpe_counts <- counts(segerstolpe)
segerstolpe_metadata <- as.data.frame(colData(segerstolpe))
rm(segerstolpe)


## Xin

xin <- readRDS("data/2018-04-15_pancreas_xin/xin.rds")
xin_counts <- normcounts(xin)
xin_metadata <- as.data.frame(colData(xin))
rm(xin)


output <- "results/2018-05-25_pancreas_processing_norm_inv"



datasets <- list(muraro = muraro_counts, 
                 segerstolpe = segerstolpe_counts, 
                 xin = xin_counts)


# Remove genes/cells ------------------------------------------------------

removeCells <- function(expData){
  
  filterCells <- function(filterParam){
    cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
    cellsToRemove
  }
  
  # Remove cells by library size
  libSizes <- colSums(expData)
  
  geneNames <- row.names(expData)
  
  mtID <- grepl("^MT-", geneNames, ignore.case = TRUE)
  rbID <- grepl("^RPL|^RPS", geneNames, ignore.case = TRUE)
  
  # Remove cells by mitochondrial expression
  mtPercent <- colSums(expData[mtID, ]) / libSizes
  
  # Remove cells by ribosomal expression
  rbPercent <- colSums(expData[rbID, ]) / libSizes
  
  
  # Remove cells
  lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), filterCells) %>%
    unlist() %>%
    unique() -> cellsToRemove
  

  expData <- expData[, -cellsToRemove]
  expData
}

datasets_cells <- lapply(datasets, removeCells)

lapply(datasets, ncol) 
lapply(datasets_cells, ncol) 

removeGenes <- function(expData, min = 1){
  
  genesZero <- apply(expData, 1, function(gene) all(gene == 0))
  expData <- expData[!genesZero,]
  
  
  ## Remove genes by quantile
  
  totalCells <- ncol(expData)
  percentCellsExpressed <- apply(expData, 1, function(x){sum(x > 0)}) / totalCells * 100
  meanExprs <- rowMeans(expData)
  percentCutoff <- min
  genesToRemove <- which(percentCellsExpressed < percentCutoff)
  
  if(length(genesToRemove) > 0){
    expData <- expData[-genesToRemove,]
  }
  
  expData
  
}

datasets_cells_genes <- lapply(datasets_cells, removeGenes)

lapply(datasets, nrow) 
lapply(datasets_cells_genes, nrow) 


# Calculate CPM -----------------------------------------------------------


get_cpm  <- function(expData) apply(expData, 2, function(x) (x/sum(x))*1000000)
cpm_data <- lapply(datasets_cells_genes, get_cpm)

filter_genes_cpm <- function(cpm, filter = 5){
  genesToKeep <- apply(cpm, 1, max) > filter
  cpm <- cpm[genesToKeep,]
  cpm
}


cpm_data <- lapply(cpm_data, filter_genes_cpm)
lapply(cpm_data, dim) %>% 
  as.data.frame(row.names = c("genes", "cells")) 

all_metadata <- list(muraro = muraro_metadata, 
                     segerstolpe = segerstolpe_metadata, 
                     xin = xin_metadata)

get_new_metadata <- function(cpm, metadata){
  cells <- colnames(cpm)
  cells_i <- match(cells, rownames(metadata))
  new_metadata <- metadata[cells_i, , drop = FALSE]
  new_metadata
}


all_metadata <- mapply(get_new_metadata, cpm_data, all_metadata, SIMPLIFY = FALSE)

mapply(function(dataset, metadata){
  all(colnames(dataset) == rownames(metadata))
}, cpm_data, all_metadata)


cpm_data$muraro
all_metadata$muraro


# Prepare testing dataset -------------------------------------------------

i <- baron_metadata$cell_type1 %in% c("alpha", "beta", "delta", "gamma")
baron_metadata <- data.frame(cell_type1 = baron_metadata$cell_type1[i], 
                             row.names = rownames(baron_metadata)[i])

baron_metadata$cell_type1 <- factor(baron_metadata$cell_type1, 
                                    levels = as.vector(unique(baron_metadata$cell_type1)))


i <- match(rownames(baron_metadata), colnames(baron_counts))
baron_counts <- baron_counts[,i]
all(rownames(baron_metadata) == colnames(baron_counts))


baron_counts_cells <- removeCells(baron_counts)
baron_counts_cells_genes <- removeGenes(baron_counts_cells, min = 1)

dim(baron_counts_cells_genes)

baron_cpm <- get_cpm(baron_counts_cells_genes)
baron_cpm <- filter_genes_cpm(baron_cpm)

dim(baron_cpm)

genes <- rownames(baron_cpm)
cells <- colnames(baron_cpm)


genes_i <- match(genes, rownames(baron_counts_cells_genes))
cells_i <- match(cells, colnames(baron_counts_cells_genes))
baron_counts_filter <- baron_counts_cells_genes[genes_i, cells_i]

all(rownames(baron_counts_filter) == genes) &
  all(colnames(baron_counts_filter) == cells)


cells_i <- match(cells, rownames(baron_metadata))
baron_metadata <- baron_metadata[cells_i,,drop = FALSE]

all(rownames(baron_metadata) == cells)

# Save training data ------------------------------------------------------

saveRDS(cpm_data, 
        file = here(output, "pancreas_cpm.RDS"))

saveRDS(all_metadata, 
        file = here(output, "pancreas_metadata.RDS"))

saveRDS(baron_cpm,
        file = here(output, "baron_cpm.RDS"))

saveRDS(baron_metadata,
        file = here(output, "baron_metadata.RDS"))
