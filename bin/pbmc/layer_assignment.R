# Script information ------------------------------------------------------

# title: Hiercharchy cell-type assignment of PBMCs
# author: José Alquicira Hernández
# date: 2019/03/14
# description: None


# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")
library("here")

# Secondary
library("Seurat")

# Set output --------------------------------------------------------------

output_dir_name <- "pbmc_assign_layers" # <------ Output directory
date <- "2019-03-13"  # <------ Date
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


cellLabels <- c("CD19+ B",
                "CD14+ Monocyte",
                "CD34+",
                "CD4+ T Helper2",
                "CD56+ NK",
                "CD8+ Cytotoxic T",
                "CD4+/CD45RO+ Memory",
                "CD8+/CD45RA+ Naive Cytotoxic",
                "CD4+/CD45RA+/CD25- Naive T",
                "CD4+/CD25 T Reg")

# Create Seurat objects ---------------------------------------------------
data <- mapply(data, 
               cellLabels, 
               FUN =  function(d, i) {
                 CreateSeuratObject(d, meta.data = data.frame(fname = rep(i, ncol(d)), 
                                                              row.names = colnames(d)))
               })



# Assign cell type tree levels --------------------------------------------

assign_cell_level_3 <- function(x){
  switch(x,
         "CD4+ T Helper2" = "non_cytotoxic",
         "CD4+/CD25 T Reg" = "non_cytotoxic",
         "CD4+/CD45RA+/CD25- Naive T" = "non_cytotoxic",
         "CD4+/CD45RO+ Memory" = "non_cytotoxic",
         "CD8+ Cytotoxic T" = "cytotoxic",
         "CD8+/CD45RA+ Naive Cytotoxic" = "cytotoxic", NA)
}

assign_cell_level_2 <- function(x){
  switch(x,
         "CD14+ Monocyte" = "Monocyte",
         "CD19+ B" = "B_cell",
         "CD34+" = "Blood_progenitor",
         "CD4+ T Helper2" = "T_cell",
         "CD4+/CD25 T Reg" = "T_cell",
         "CD4+/CD45RA+/CD25- Naive T" = "T_cell",
         "CD4+/CD45RO+ Memory" = "T_cell",
         "CD56+ NK" = "Natural_killer",
         "CD8+ Cytotoxic T" = "T_cell",
         "CD8+/CD45RA+ Naive Cytotoxic" = "T_cell",
         "Dendritic" = "Dendritic_cell")
}

assign_cell_level_1 <- function(x){
  switch(x,
         "Dendritic_cell" = "Myeloid_cell",
         "Monocyte" = "Myeloid_cell",
         "B_cell" = "Lymphoid_cell",
         "T_cell" = "Lymphoid_cell",
         "Natural_killer" = "Lymphoid_cell", x)
}


addCellLevel3 <- function(x){
  x@meta.data$fname %>% 
    as.character() %>% 
    sapply(assign_cell_level_3) %>% 
    as.factor() -> label
  x@meta.data$level3 <- label
  x
}

addCellLevel2 <- function(x){
  x@meta.data$fname %>% 
    as.character() %>% 
    sapply(assign_cell_level_2) %>% 
    as.factor() -> label
  x@meta.data$level2 <- label
  x
}

addCellLevel1 <- function(x){
  x@meta.data$level2 %>% 
    as.character() %>% 
    sapply(assign_cell_level_1) %>% 
    as.factor() -> label
  x@meta.data$level1 <- label
  x
}

data %>% 
  lapply(addCellLevel3) %>% 
  lapply(addCellLevel2) %>% 
  lapply(addCellLevel1) -> data


# Save data ---------------------------------------------------------------

# Merge datasets

# All data
combined <- data[[1]]
for (i in 2:length(x = data)) {
  combined <- MergeSeurat(object1 = combined, object2 = data[[i]])
}

combined

cat("Saving results...\n")

saveRDS(combined, file = here(output, "pbmc.RDS"))


# Subset data
combined@meta.data %>% 
  rownames_to_column("barcode") %>% 
  select(barcode, fname) %>% 
  split(.$fname) %>% 
  lapply(head, 100) %>% 
  reduce(rbind) %>% 
  pull("barcode") -> subCells

subCombined <- SubsetData(combined, cells.use = subCells, do.clean = TRUE)
saveRDS(subCombined, file = here(output, "pbmc_sub.RDS"))



# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
