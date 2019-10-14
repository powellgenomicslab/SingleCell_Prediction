# Script information ------------------------------------------------------
# title: Process pure bead-enriched PBMCs
# author: José Alquicira Hernández
# date: 2018-10-30
# description: 
#   - scPred version (d7661592)


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary
library("Seurat")
library("scPred")
library("mygene")


# Set output --------------------------------------------------------------


output_dir_name <- "pbmc_pure_processed" # <------ Output directory
date <- "2018-11-02"
output <- file.path("results", paste(date, output_dir_name, sep = "_"))

if (!dir.exists(output)) {
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input <- file.path("data", "2018-11-01_pbmc_pure") # <------ Input directory
filename <- "all_pure_select_11types.rds" # <------ Input file
pbmc <- readRDS(here(input, filename))



# Process data ------------------------------------------------------------

# Get gene ids

## All datasets contain the same gene ids 
pbmc$all_data %>% 
  lapply(function(x) x$hg19$genes) %>% 
  unlist() %>% 
  table() %>% 
  all(. == 10)

## Get gene ids
pbmc$all_data %>% 
lapply(function(ct){ct$hg19$genes}) %>% 
  unlist() %>% 
  unique() -> geneList

## Get gene symbols
symbols <- getGenes(geneid = geneList, fields = "symbol")

## Get duplicated symbols
symDup  <- symbols[duplicated(symbols$query), "query"]
i <- which(symbols$query %in% symDup)
symbolsDup <- symbols[i,]

## Bind duplicated symbols
split(symbolsDup, symbolsDup$query) %>% 
  lapply("[[", "symbol") %>% 
  lapply(function(x) paste0(x, collapse = "_")) -> symbolsFormat

## Order merged symbols according to original order
i <- match(unique(symbolsDup$query), names(symbolsFormat))
symbolsFormat <- symbolsFormat[i]

## Remove last duplicated gene id
i <- which(symbols$query %in% symDup)

rmIndex <- i[seq(2, length(i), 2)]
symbols <- symbols[-rmIndex,]

## Assing new symbols
i <- which(symbols$query %in% names(symbolsFormat))
symbols$symbol[i] <- unlist(symbolsFormat)

## Make sure all new gene ids and symbols are in the original order
all(symbols$query == geneList)




## Combine symbols with gene ids
combineSymbol <- function(query, symbol){
  if(!is.na(symbol)){
    paste(query, symbol, sep = "_")
  }else{
    query
  }
}


genes <- mapply(combineSymbol, 
                symbols$query, 
                symbols$symbol, 
                USE.NAMES = FALSE)


# Create expression matrices with corresponding barcode ids and gene names

setUpMats <- function(ct){
  mt <- ct$hg19$mat 
  mt <- t(mt)
  
  colnames(mt) <- ct$hg19$barcodes
  rownames(mt) <- genes
  mt
}

mats <- lapply(pbmc$all_data, setUpMats)


# Create Seurat objects

mats <- lapply(mats, setObject)

## assign cell type information as metadata
mats <- mapply(function(ct, md){ ct@meta.data$cellType <- md; ct }, 
               mats, 
               as.list(as.character(pbmc$all_metrics$description)), 
               SIMPLIFY = FALSE)

## Add metrics (mitochondrial and ribosomal content)
mats <- lapply(mats, addQC)

## Plot QC metrics
qcPlots <- mapply(plotQC, mats, as.list(as.character(pbmc$all_metrics$description)), SIMPLIFY = FALSE)
nGenesFilter <- list(14000, 5000, 3900, 4000, 4000, 5000, 5000, 5000, 4000, 6000)

## Filter cells
mats <- mapply(filterQC, mats, low = -Inf, high = nGenesFilter, SIMPLIFY = FALSE)


## Normalize and scale data
mats <- lapply(mats, processSeurat)


saveRDS(mats, here(output, "pbmc_pure_processed.RDS"))

lapply(mats, function(x) ncol(x@data)) %>% unlist() %>% sum()

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))



