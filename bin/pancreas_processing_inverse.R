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


output <- "results/2018-05-14_pancreas_processing"

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


baron_metadata %>% 
  as.data.frame() %>% 
  select(human, cell_type1) %>% 
  rownames_to_column("id") %>% 
  filter(cell_type1 %in% c("alpha", "beta", "delta", "gamma")) %>% 
  mutate(cell_type1 = factor(cell_type1, 
                             c("alpha", "beta", "delta", "gamma"))) %>% 
  column_to_rownames("id") -> baron_metadata


i <- match(rownames(baron_metadata), colnames(cpm_data$baron))
baron_cpm <- cpm_data$baron[,i]
all(colnames(baron_cpm) == rownames(baron_metadata))



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


# Create metadata
cellType <- c(muraro_metadata$cell_type1, segerstolpe_metadata$cell_type1, xin_metadata$cell_type1)
cell_id <- c(rownames(muraro_metadata), rownames(segerstolpe_metadata), rownames(xin_metadata))
training_metadata <- data.frame(cellType = cellType, 
                                dataset = c(rep("muraro", nrow(muraro_metadata)),
                                            rep("segerstolpe", nrow(segerstolpe_metadata)),
                                            rep("xin", nrow(xin_metadata))),
                                row.names = cell_id)

all(rownames(training_metadata) == colnames(training))



# Save training data ------------------------------------------------------

saveRDS(baron_cpm, 
        file = here(output, "baron_training_cpm.RDS"))

saveRDS(baron_metadata, 
        file = here(output, "baron_training_metadata.RDS"))

saveRDS(training, 
        file = here(output, "testing_cpm.RDS"))

saveRDS(training_metadata, 
        file = here(output, "testing_metadata.RDS"))


