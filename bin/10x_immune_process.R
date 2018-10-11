# Script information ------------------------------------------------------

# title: Process immune dataset from 10x
# subtitle: 68k dataset
# author: José Alquicira Hernández
# date: 2018-10-10
# description: Integrate expression data and metadata, and perform PCA


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("data.table")
library("here")

# Secondary
library("Seurat")

# Set output --------------------------------------------------------------


output_dir_name <- "10x_immune_processed" # <------ Output directory

date <- format(Sys.Date(), format = "%Y-%m-%d_")
date <- "2018-10-10_"
output <- file.path("results", paste0(date, output_dir_name))
if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

directory <- file.path("data", "2018-10-02_10x_immune_cells") # <------ Input directory


# Read file

filename <- file.path("fresh_68k_pbmc_donor_a_filtered_matrices_mex","hg19")
pbmc <- Read10X(here(directory, filename))

filename <- "68k_pbmc_barcodes_annotation.tsv"
meta.data <- fread(here(directory, filename))

all(colnames(pbmc) == meta.data$barcodes)

# Create seurat object ----------------------------------------------------

# - Keep all genes expressed in ~0.1% of the data. Keep all cells with at
# least 200 detected genes

pbmc <- CreateSeuratObject(pbmc, project = "10X_68k_PBMC")
pbmc@meta.data$cellType <- as.factor(meta.data$celltype)

# Normalize data ----------------------------------------------------------

pbmc <- NormalizeData(object = pbmc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)


# Find variable genes -----------------------------------------------------

pbmc <- FindVariableGenes(object = pbmc, 
                          mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, 
                          x.high.cutoff = 3, 
                          y.cutoff = 0.5)

cat(length(x = pbmc@var.genes), "variable genes identified\n")


# Scale data --------------------------------------------------------------

pbmc <- ScaleData(object = pbmc)



# Principal component analysis --------------------------------------------

pbmc <- RunPCA(object = pbmc, 
               pc.genes = pbmc@var.genes, 
               do.print = TRUE, 
               pcs.print = 1:5, 
               genes.print = 5)


PCElbowPlot(object = pbmc)
# According to the previous plot, 12-13 principal components are relevant

# Run Non-linear dimensional reduction (tSNE) -----------------------------

pbmc <- RunTSNE(object = pbmc, dims.use = 1:13, do.fast = TRUE)

png(filename = here(output, "pbmc_pca.png"), width = 1450, height = 800, res = 250)
PCAPlot(object = pbmc, group = "cellType")
dev.off()

png(filename = here(output, "pbmc_pca_2-3.png"), width = 1450, height = 800, res = 250)
PCAPlot(object = pbmc, group = "cellType", 2, 3)
dev.off()

png(filename = here(output, "pbmc_pca_1-3.png"), width = 1450, height = 800, res = 250)
PCAPlot(object = pbmc, group = "cellType", 1, 3)
dev.off()

# Save data ---------------------------------------------------------------

saveRDS(pbmc, file.path(output, "pbmc.RDS"))
# pbmc <- readRDS(file.path(output, "pbmc.RDS"))


# Session info ------------------------------------------------------------

options(width = 70)
devtools::session_info()