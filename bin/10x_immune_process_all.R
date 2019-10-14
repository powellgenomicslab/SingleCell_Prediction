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


# Do all matrices share the same genes?

data %>% 
  lapply(function(x) rownames(x)) %>% 
  reduce(cbind) %>%
  apply(1, unique) %>% 
  length() %>% 
  `==`(nrow(data$cd14_monocytes))




# Create metadata ---------------------------------------------------------

data %>% 
  lapply(ncol) %>% 
  mapply(function(n, name){rep(name, n)}, ., names(.), SIMPLIFY = FALSE) -> labels

data %>% 
  lapply(colnames) -> ids

metadata <- mapply(function(ids, labels){data.frame(labels, row.names = ids)}, 
                   ids, 
                   labels, 
                   SIMPLIFY = FALSE) %>% 
  reduce(rbind)


# Is metadata in the same order as the gene expression data?

data <- data %>% 
  reduce(cbind)

all(colnames(data) == rownames(metadata))


# QC ----------------------------------------------------------------------

# Remove genes

## Remove genes with zero counts across all cells

genesZero <- rowSums(data) == 0
data <- data[!genesZero,]


## Remove genes by quantile

totalCells <- ncol(data)
percentCellsExpressed <- apply(data, 1, function(x){sum(x > 0)}) / totalCells * 100
meanExprs <- rowMeans(data)
percentCutoff <- 1
genesToRemove <- which(percentCellsExpressed < percentCutoff)

if(length(genesToRemove) > 0){
  data <- data[-genesToRemove,]
}



# Normalize data and dimensionality reduction -----------------------------


data <- CreateSeuratObject(raw.data = data, meta.data = metadata) %>% 
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>% 
  FindVariableGenes(mean.function = ExpMean, 
                    dispersion.function = LogVMR, 
                    x.low.cutoff = 0.0125, 
                    x.high.cutoff = 3, 
                    y.cutoff = 0.5) %>% 
  ScaleData(vars.to.regress = "nUMI") %>% 
  RunPCA(pc.genes = .@var.genes, 
         do.print = TRUE, 
         pcs.print = 1:5, 
         genes.print = 5)



PCAPlot(data, group = "labels")




# Save results ------------------------------------------------------------

saveRDS(data, file = here(output, "pbmc_processed.RDS"))




# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
