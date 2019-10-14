library(tidyverse)
library(here)
library(SingleCellExperiment)

xin <- readRDS("data/2018-04-15_pancreas_xin/xin.rds")
output <- "results/2018-04-15_xin_processing"


expData <- normcounts(xin)
metadata <- colData(xin)
rm(xin)


if(!all(rownames(metadata) == colnames(expData))){
  stop("Cell ids do not match in metadata and expression data")
}  


cat("# Raw summary\n", sep = "\n", file = here(output, "data_summary.txt"))
cat("\nNumber of genes:", nrow(expData), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nNumber of cells:", ncol(expData), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  


# Remove cells ------------------------------------------------------------

expData <- round(expData)


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


newMetadata <- metadata[match(colnames(cpm), rownames(metadata)), , drop = FALSE]

if(!all(colnames(cpm) == rownames(newMetadata))){
  stop("Expression data and metadata cell ids do not match")
}



cat("\n# Final summary", "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  

cat("\nFinal number of cells:", ncol(cpm), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nFinal number of genes:", nrow(cpm), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  

saveRDS(cpm, 
        file = here(output, "xin_processed_cpm.RDS"))

saveRDS(newMetadata, 
        file = here(output, "xin_processed_metadata.RDS"))


options(width = 120)
devtools::session_info() %>% 
  capture.output() %>% 
  cat("\nSession Info\n", .,  sep = "\n", file = here(output, "data_summary.txt"), append = TRUE)

