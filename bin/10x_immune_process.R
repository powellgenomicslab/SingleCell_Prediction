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
source(here("bin","quickSeurat.R"))

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

# filename <- file.path("fresh_68k_pbmc_donor_a_filtered_matrices_mex","hg19")
# pbmc <- Read10X(here(directory, filename))

filename <- file.path("pbmc68k_data.rds")
data <- readRDS(here(directory, filename))
pbmc <- data$all_data$`17820`$hg19$mat
rownames(pbmc) <- data$all_data$`17820`$hg19$barcodes
colnames(pbmc) <- paste(data$all_data$`17820`$hg19$genes, 
                        data$all_data$`17820`$hg19$gene_symbols, 
                        sep = "_")

pbmc <- t(pbmc)



filename <- "68k_pbmc_barcodes_annotation.tsv"
meta.data <- fread(here(directory, filename))

all(colnames(pbmc) == meta.data$barcodes)

# Create seurat object ----------------------------------------------------

# - Keep all genes expressed in ~0.1% of the data. Keep all cells with at
# least 200 detected genes

pbmc <- CreateSeuratObject(pbmc, project = "10X_68k_PBMC")
pbmc@meta.data$cellType <- as.factor(meta.data$celltype)



# Add QC ------------------------------------------------------------------

pbmc <- addQC(pbmc)
plotQC(pbmc, "All cell types")
pbmc <- filterQC(pbmc, attr = "percent.mito", low = -Inf, high = 0.15)
pbmc <- filterQC(pbmc, attr = "percent.ribo", low = 0.1, high = Inf)
pbmc <- filterQC(pbmc, attr = "nGene", low = -Inf, high = 2000)
plotQC(pbmc, "All cell types")



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

# According to the previous plot, 11 principal components are relevant

# Run Non-linear dimensional reduction (tSNE) -----------------------------

pbmc <- RunTSNE(object = pbmc, dims.use = 1:11, do.fast = TRUE)
TSNEPlot(pbmc, group = "cellType", colors.use = as.character(jcolors::jcolors("pal8"))) 
pbmc <- RunUMAP(pbmc)



# i <- pbmc@meta.data$cellType == "CD19+ B"
# pca <- as.data.frame(pbmc@dr$pca@cell.embeddings[,c(1,2)])
# 
# j <- pca$PC1 > -2 & pca$PC2 < 5
# subpca <- pca[i & j, ]
#  k <- rownames(subpca)
# subpbmc <- as.matrix(pbmc@raw.data[,colnames(pbmc@raw.data) %in% k])
# write.table(subpbmc, file = "~/Desktop/coco.txt", quote = FALSE, sep = "\t", row.names = TRUE)

p <- plotPCA(pbmc, group = "cellType")
ggsave(filename = here(output, "pbmc_pca.png"),  width = 7.5, height = 5, dpi = 300)

p <- plotPCA(pbmc, 2, 3, group = "cellType")
ggsave(filename = here(output, "pbmc_pca_2-3.png"),  width = 7.5, height = 5, dpi = 300)

p <- plotPCA(pbmc, 1, 3, group = "cellType")
ggsave(filename = here(output, "pbmc_pca_1-3.png"),  width = 7.5, height = 5, dpi = 300)



# Save data ---------------------------------------------------------------

saveRDS(pbmc, file.path(output, "pbmc.RDS"))

pbmc@meta.data %>% 
  rownames_to_column("id") %>% 
  group_by(cellType) %>% 
  summarise(n = n())


set.seed(66)
pbmc@meta.data %>% 
  rownames_to_column("id") %>% 
  group_by(cellType) %>% 
  sample_n(97) -> pbmc_sub


pbmc_sub <- SubsetData(pbmc, cells.use =  pbmc_sub$id, do.clean = TRUE)
saveRDS(pbmc_sub, file.path(output, "pbmc_sub.RDS"))


# Session info ------------------------------------------------------------

options(width = 70)
devtools::session_info()