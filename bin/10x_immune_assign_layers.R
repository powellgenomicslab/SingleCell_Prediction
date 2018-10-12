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
filename <- "pbmc.RDS" # <------ Input file


# Read file

data <- readRDS(here(input, filename))



# Assign cell type tree levels --------------------------------------------


assign_cell_level <- function(x){
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

data@meta.data$cellType1 <- data@meta.data$cellType %>% 
  as.character() %>% 
  sapply(assign_cell_level) %>% 
  as.factor()


assign_cell_level_2 <- function(x){
  switch(x,
         "Dendritic_cell" = "Myeloid_cell",
         "Monocyte" = "Myeloid_cell",
         "B_cell" = "Lymphoid_cell",
         "T_cell" = "Lymphoid_cell",
         "Natural_killer" = "Lymphoid_cell", x)
}

data@meta.data$cellType2 <- data@meta.data$cellType1 %>% 
  as.character() %>% 
  sapply(assign_cell_level_2) %>% 
  as.factor()


# PCAPlot(data, group = "cellType2", cols.use = set_names(jcolors::jcolors(), NULL)) +
#   theme_bw() +
#   guides(color = guide_legend(title = "Cell type"))


# Create partitions -------------------------------------------------------

set.seed(66)
i <- createDataPartition(data@meta.data$cellType2, p = 0.75, list = FALSE)
j <- data@cell.names[i]
k <- data@cell.names[-i]


train_data <- SubsetData(data, cells.use = j, do.clean = TRUE)
pred_data  <- SubsetData(data, cells.use = k, do.clean = TRUE)

# Train
train_data <- train_data %>% 
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>% 
  FindVariableGenes(mean.function = ExpMean, 
                    dispersion.function = LogVMR, 
                    x.low.cutoff = 0.0125, 
                    x.high.cutoff = 3, 
                    y.cutoff = 0.5) %>% 
  ScaleData() %>% 
  RunPCA(pc.genes = .@var.genes, 
         do.print = TRUE, 
         pcs.print = 1:5, 
         genes.print = 5)

PCAPlot(train_data, group = "cellType2")
saveRDS(train_data, file = here(output, "train_data.RDS"))


# Prediction

pred_data <- pred_data %>% 
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>% 
  FindVariableGenes(mean.function = ExpMean, 
                    dispersion.function = LogVMR, 
                    x.low.cutoff = 0.0125, 
                    x.high.cutoff = 3, 
                    y.cutoff = 0.5) %>% 
  ScaleData() %>% 
  RunPCA(pc.genes = .@var.genes, 
         do.print = TRUE, 
         pcs.print = 1:5, 
         genes.print = 5)

PCAPlot(pred_data, group = "cellType2")
saveRDS(pred_data, file = here(output, "pred_data.RDS"))



# Session info ------------------------------------------------------------

options(width = 70)
devtools::session_info()
