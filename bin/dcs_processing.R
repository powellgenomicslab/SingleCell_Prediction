# Import libraries --------------------------------------------------------

library("tidyverse")
library("here")
library("stringr")

# Set input/output variables ----------------------------------------------

expDataDir <- "data/2017-09-29_dendritic_cells/GSE89232_expMatrix.txt"
output <- "results/2018-04-08_dcs_processing"


# Read data ---------------------------------------------------------------
expData <- read.table(file = expDataDir, header = TRUE, sep = "\t")
names(expData)



expData <- round(expData) -> expDataCounts

cat("# Raw summary\n", sep = "\n", file = here(output, "data_summary.txt"))
cat("\nNumber of genes:", nrow(expData), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nNumber of cells:", ncol(expData), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  


# Remove cells ------------------------------------------------------------


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
        file = here(output, "dcs_processed_counts.RDS"))

saveRDS(cpm, 
        file = here(output, "dcs_processed_cpm.RDS"))


options(width = 120)
devtools::session_info() %>% 
  capture.output() %>% 
  cat("\nSession Info\n", .,  sep = "\n", file = here(output, "data_summary.txt"), append = TRUE)
