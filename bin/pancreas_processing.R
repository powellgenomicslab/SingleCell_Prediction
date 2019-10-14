library(tidyverse)
library(here)


# Read data ---------------------------------------------------------------

## Baron
baron <- readRDS("data/2018-04-15_pancreas_baron/baron-human.rds")
baron_counts <- counts(baron)
baron_metadata <- colData(baron)
rm(baron)

## Muraro
muraro <- readRDS("data/2018-04-15_pancreas_muraro/muraro.rds")
muraro_counts <- normcounts(muraro)
muraro_metadata <- colData(muraro)
rm(muraro)

## Segerstolpe

segerstolpe <- readRDS("data/2018-04-15_pancreas_segerstolpe/segerstolpe.rds")
segerstolpe_counts <- counts(segerstolpe)
segerstolpe_metadata <- colData(segerstolpe)
rm(segerstolpe)


## Xin

xin <- readRDS("data/2018-04-15_pancreas_xin/xin.rds")
xin_counts <- normcounts(xin)
xin_metadata <- colData(xin)
rm(xin)


output <- "results/2018-04-15_pancreas_processing"

table(baron_metadata$cell_type1)
table(muraro_metadata$cell_type1)
table(segerstolpe_metadata$cell_type1)
table(xin_metadata$cell_type1)


# Calculate CPM -----------------------------------------------------------


get_cpm  <- function(expData) apply(expData, 2, function(x) (x/sum(x))*1000000)


cpm_data <- lapply(list(baron = baron_counts, 
                        muraro = muraro_counts, 
                        segerstolpe = segerstolpe_counts, 
                        xin = xin_counts),
                   get_cpm)

#rm(baron_counts, muraro_counts, segerstolpe_counts, xin_counts)



# Merge reference dataset -------------------------------------------------

# Normalize gene names
cpm_data$muraro %>% 
  rownames() %>% 
  str_split("__") %>% 
  lapply("[", 1) %>% 
  unlist() -> muraro_genenames

rownames(cpm_data$muraro) <- muraro_genenames
rownames(muraro_counts) <- muraro_genenames
rm(muraro_genenames)

# Get shared genes
muraro_segerstolpe_genes <- intersect(rownames(cpm_data$muraro), 
                                      rownames(cpm_data$segerstolpe))

reference_genes <- intersect(muraro_segerstolpe_genes, rownames(cpm_data$xin))


get_reference_genes <- function(expData, genes){
  expData[match(genes, rownames(expData)), ]
}


cpm_reference <- lapply(cpm_data[c("muraro", "segerstolpe", "xin")], 
                        get_reference_genes, 
                        reference_genes)

# Merge datasets for training
cpm_reference %>% 
  Reduce(cbind, .) -> training

# Get reference cells
reference_cells <- intersect(muraro_metadata$cell_type1, segerstolpe_metadata$cell_type1)
reference_cells <- intersect(xin_metadata$cell_type1, reference_cells)

# Create metadata
cellType <- c(muraro_metadata$cell_type1, segerstolpe_metadata$cell_type1, xin_metadata$cell_type1)
cell_id <- c(rownames(muraro_metadata), rownames(segerstolpe_metadata), rownames(xin_metadata))
training_metadata <- data.frame(cellType = cellType, 
                                dataset = c(rep("muraro", nrow(muraro_metadata)),
                                            rep("segerstolpe", nrow(segerstolpe_metadata)),
                                            rep("xin", nrow(xin_metadata))),
                                row.names = cell_id)

all(rownames(training_metadata) == colnames(training))

# Filter cells
training_metadata <- training_metadata[training_metadata$cellType %in% reference_cells, , drop = FALSE]
training_metadata$cellType <- factor(training_metadata$cellType, levels = reference_cells)


# Filter cpm data
training <- training[ ,colnames(training) %in% rownames(training_metadata)]
all(rownames(training_metadata) == colnames(training))



# Filter count data

cells <- colnames(training)
genes <- rownames(training)


getCounts <- function(x){
genes_i <- match(genes, rownames(x)) 
cells_i <- colnames(x) %in% cells
x <- x[genes_i, cells_i]
x
}


list(muraro = muraro_counts, segerstolpe = segerstolpe_counts, xin = xin_counts) %>%
  lapply(getCounts) %>% 
  Reduce(cbind, .) -> training_counts
  
all(colnames(training) == colnames(training_counts))
all(rownames(training) == rownames(training_counts))


# Save training data ------------------------------------------------------

saveRDS(training, 
        file = here(output, "pancreas_cpm.RDS"))

saveRDS(training_counts, 
        file = here(output, "pancreas_counts.RDS"))

saveRDS(training_metadata, 
        file = here(output, "pancreas_metadata.RDS"))
