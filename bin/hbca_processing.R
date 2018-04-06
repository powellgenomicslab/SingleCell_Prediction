# Import libraries --------------------------------------------------------

library("tidyverse")
library("here")
library("stringr")

# Set input/output variables ----------------------------------------------

expDataDir <- "data/2017-07-17_blood_atlas/RAW/raw_expression_matrix.txt"
metadataDir <- "data/2017-07-17_blood_atlas/RAW/metadata.txt"
output <- "results/2018-04-06_hbca_processing"

  
# Read data ---------------------------------------------------------------

expData <- read.table(here(expDataDir),
                      header = TRUE,
                      sep = "\t")

  
metadata <- read.table(here(metadataDir),
                         header = TRUE,
                         sep = "\t",
                         skip = 1)
  
# Assign rownames to gene expression matrix
rownames(expData) <- make.names(expData$Gene.ID, unique = TRUE)
expData$Gene.ID <- NULL

cat("# Raw summary\n", sep = "\n", file = here(output, "data_summary.txt"))
cat("\nNumber of genes:", nrow(expData), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nNumber of cells:", ncol(expData), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  

  
# Remove cells ------------------------------------------------------------

# Rename cells to match with provided metadata
expData %>% 
  names() %>% 
  str_split("_") %>% 
  lapply(function(x) paste0(x[seq_len(length(x) - 1)], collapse = "_")) %>% 
  unlist() -> cellNames 
  
names(expData) <- cellNames

# Match cell ids in gene expression data with metadata
expData <- expData[, names(expData) %in% metadata$TYPE]
  
dim(expData)


# Order metadata according to gene expression matrix
  
metadata <- metadata[match(names(expData), metadata$TYPE), ]
rownames(metadata) <- metadata$TYPE
metadata$TYPE <- NULL

if(!all(rownames(metadata) == names(expData))){
  stop("Cell ids do not match in metadata and expression data")
}
  

expData <- round(expData) -> expDataCounts


filterCells <- function(filterParam){
  cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
  cellsToRemove
}

# Remove cells by library size
libSizes <- colSums(expData)

geneNames <- row.names(expData)

mtID <- grepl("^MT-", geneNames, ignore.case = TRUE)
rbID <- grepl("^RPL|^RPS", geneNames, ignore.case = TRUE)

# Remove cells by mitochondrial expression
mtPercent <- colSums(expData[mtID, ]) / libSizes

# Remove cells by ribosomal expression
rbPercent <- colSums(expData[rbID, ]) / libSizes


# Remove cells
lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), filterCells) %>% 
  unlist() %>% 
  unique() -> cellsToRemove

expData <- expData[, -cellsToRemove]
cat("\n# Cell removal", "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nNumber of removed cells:", length(cellsToRemove), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nNumber of cells left:", ncol(expData), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  



# Remove genes ------------------------------------------------------------

## Remove genes with zero counts across all cells

genesZero <- apply(expData, 1, function(gene) all(gene == 0))
expData <- expData[!genesZero,]


## Remove genes by quantile

totalCells <- ncol(expData)
percentCellsExpressed <- apply(expData, 1, function(x){sum(x > 0)}) / totalCells * 100
meanExprs <- rowMeans(expData)
# (meanExprs + 1) %>% log2() %>% density() %>% plot()
percentCutoff <- 1
genesToRemove <- which(percentCellsExpressed < percentCutoff)

if(length(genesToRemove) > 0){
  expData <- expData[-genesToRemove,]
}
# rowMeans(expData) %>% log2() %>% density() %>% plot()



# Calculate CPM -----------------------------------------------------------


cpm  <- apply(expData, 2, function(x) (x/sum(x))*1000000)

genesToKeep <- apply(cpm, 1, max) > 5 
cpm <- cpm[genesToKeep,]


cat("\n# Gene removal", "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nNumber of genes left:", nrow(cpm), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  


newMetadata <- metadata[match(colnames(cpm), rownames(metadata)),]

if(!all(colnames(cpm) == rownames(newMetadata))){
  stop("Expression data and metadata cell ids do not match")
}

cat("\n# Final summary", "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  

cat("\nFinal number of cells:", ncol(cpm), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nFinal number of genes:", nrow(cpm), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  


genesMatch <- match(rownames(cpm), rownames(expDataCounts))
cellsMatch <- match(colnames(cpm), colnames(expDataCounts))

expDataCounts <- expDataCounts[genesMatch, cellsMatch]

if(!(all(rownames(expDataCounts) == rownames(cpm)) & all(colnames(expDataCounts) == colnames(cpm)))){
  stop("Normalized matrix is not equivalent to count matrix")
}


saveRDS(expDataCounts, 
        file = here(output, "hbca_processed_counts.RDS"))

saveRDS(cpm, 
        file = here(output, "hbca_processed_cpm.RDS"))

saveRDS(newMetadata, 
        file = here(output, "hbca_processed_metadata.RDS"))


options(width = 120)
devtools::session_info() %>% 
  capture.output() %>% 
  cat("\nSession Info\n", .,  sep = "\n", file = here(output, "data_summary.txt"), append = TRUE)
