library(tidyverse)
library(here)


# Read data ---------------------------------------------------------------

# muraro
muraro_cpm <- readRDS("results/2018-04-15_muraro_processing/muraro_processed_cpm.RDS")
muraro_metadata <- readRDS("results/2018-04-15_muraro_processing/muraro_processed_metadata.RDS")

# segerstolpe
segerstolpe_cpm <- readRDS("results/2018-04-15_segerstolpe_processing/segerstolpe_processed_cpm.RDS")
segerstolpe_metadata <- readRDS("results/2018-04-15_segerstolpe_processing/segerstolpe_processed_metadata.RDS")

# xin
xin_cpm <- readRDS("results/2018-04-15_xin_processing/xin_processed_cpm.RDS")
xin_metadata <- readRDS("results/2018-04-15_xin_processing/xin_processed_metadata.RDS")

output <- "results/2018-04-16_pancreas_reference_processing"


# Merge reference dataset -------------------------------------------------

# Merge metadata

getMetadata <- function(m){
  data.frame(cellType = as.character(m$cell_type1), row.names = rownames(m))
}

list(muraro_metadata, segerstolpe_metadata, xin_metadata) %>% 
  lapply(getMetadata) %>% 
  mapply(function(x,y){x$dataset <- y; x}, ., c("muraro", "segerstolpe", "xin"), SIMPLIFY = FALSE) %>% 
  reduce(rbind) %>% 
  as.data.frame() -> training_metadata


# Normalize gene names for muraro dataset
muraro_cpm %>% 
  rownames() %>% 
  str_split("__") %>% 
  lapply("[", 1) %>% 
  unlist() -> muraro_genenames

rownames(muraro_cpm) <- muraro_genenames
rm(muraro_genenames)

# Get shared genes
reference_genes <- intersect(rownames(muraro_cpm), 
                                      rownames(segerstolpe_cpm))

reference_genes <- intersect(reference_genes, rownames(xin_cpm))


get_reference_genes <- function(expData, genes){
  expData[match(genes, rownames(expData)), ]
}


cpm_reference <- lapply(list(muraro = muraro_cpm,  segerstolpe = segerstolpe_cpm, xin = xin_cpm), 
                        get_reference_genes, 
                        reference_genes)

# Merge datasets for training
cpm_reference %>% 
  Reduce(cbind, .) -> training

# Get reference cells
reference_cells <- intersect(muraro_metadata$cell_type1, segerstolpe_metadata$cell_type1)
reference_cells <- intersect(xin_metadata$cell_type1, reference_cells)

# Create metadata

all(rownames(training_metadata) == colnames(training))

# Filter cells
training_metadata <- training_metadata[training_metadata$cellType %in% reference_cells, , drop = FALSE]
training_metadata$cellType <- factor(training_metadata$cellType, levels = reference_cells)


# Filter cpm data
training <- training[ ,colnames(training) %in% rownames(training_metadata)]
all(rownames(training_metadata) == colnames(training))


# Save training data ------------------------------------------------------

saveRDS(training, 
        file = here(output, "pancreas_processed_cpm.RDS"))

saveRDS(training_metadata, 
        file = here(output, "pancreas_processed_metadata.RDS"))

