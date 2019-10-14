# Script information ------------------------------------------------------

# title: Assign prediction labels to cells depending on hierarchy
# author: José Alquicira Hernández
# date: 2018-10-12
# description: None


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("data.table")
library("here")

# Secondary
library("Seurat")
library("caret")

# Set output --------------------------------------------------------------


output_dir_name <- "10x_immune_tree" # <------ Output directory

date <- format(Sys.Date(), format = "%Y-%m-%d_")
date <- "2018-10-12_" # <------ Fix date

output <- file.path("results", paste0(date, output_dir_name))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input    <- file.path("results", "2018-10-10_10x_immune_processed") # <------ Input directory
filename <- "pbmc_sub.RDS" # <------ Input file


# Read file

data <- readRDS(here(input, filename))



# Assign cell type tree levels --------------------------------------------


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

data@meta.data$cellType2 <- data@meta.data$cellType %>% 
  as.character() %>% 
  sapply(assign_cell_level_2) %>% 
  as.factor()


assign_cell_level_1 <- function(x){
  switch(x,
         "Dendritic_cell" = "Myeloid_cell",
         "Monocyte" = "Myeloid_cell",
         "B_cell" = "Lymphoid_cell",
         "T_cell" = "Lymphoid_cell",
         "Natural_killer" = "Lymphoid_cell", x)
}

data@meta.data$cellType1 <- data@meta.data$cellType2 %>% 
  as.character() %>% 
  sapply(assign_cell_level_1) %>% 
  as.factor()



# Layer 1. Train and test -------------------------------------------------

set.seed(66)
i <- createDataPartition(data@meta.data$cellType1, p = 0.75, list = FALSE)
j <- data@cell.names[i]
k <- data@cell.names[-i]


train_data <- SubsetData(data, cells.use = j, do.clean = TRUE)
pred_data  <- SubsetData(data, cells.use = k, do.clean = TRUE)

# Train
saveRDS(train_data, file = here(output, "train_data.RDS"))


# Prediction
saveRDS(pred_data, file = here(output, "pred_data.RDS"))


# Layer 2. Train ----------------------------------------------------------

# Myeloid
targetCells <- c("Dendritic_cell", "Monocyte")
cells <- rownames(train_data@meta.data[train_data@meta.data$cellType2 %in% 
                             targetCells,])
train_myeloid <- SubsetData(train_data, cells.use = cells, do.clean = TRUE)
train_myeloid@meta.data$cellType2 <- factor(train_myeloid@meta.data$cellType2, 
                                            levels = targetCells)
saveRDS(train_myeloid, file = here(output, "train_myeloid.RDS"))

# Lymphoid

targetCells <- c("B_cell", "T_cell", "Natural_killer")
cells <- rownames(train_data@meta.data[train_data@meta.data$cellType2 %in% 
                             targetCells,])
train_lymphoid <- SubsetData(train_data, cells.use = cells, do.clean = TRUE)
train_lymphoid@meta.data$cellType2 <- factor(train_lymphoid@meta.data$cellType2, 
                                            levels = targetCells)
saveRDS(train_lymphoid, file = here(output, "train_lymphoid.RDS"))


# Layer 3. Train ----------------------------------------------------------

targetCells <- c("CD4+ T Helper2", "CD4+/CD25 T Reg", 
                "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory",
                "CD8+ Cytotoxic T", "CD8+/CD45RA+ Naive Cytotoxic")

cells <- rownames(train_data@meta.data[train_data@meta.data$cellType %in% 
                                         targetCells,])
train_tcells<- SubsetData(train_data, cells.use = cells, do.clean = TRUE)
train_tcells@meta.data$cellType2 <- factor(train_tcells@meta.data$cellType2, 
                                             levels = targetCells)
saveRDS(train_lymphoid, file = here(output, "train_tcells.RDS"))



# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))