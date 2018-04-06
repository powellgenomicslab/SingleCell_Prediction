# Import libraries --------------------------------------------------------

library("tidyverse")
library("here")
library("data.table")
library("stringr")



# Set input/output variables ----------------------------------------------

tumorFile <- "data/2017-10-95_colon/GSE81861_CRC_tumor_epithelial_cells_COUNT.csv"
normalFile <- "data/2017-10-95_colon/GSE81861_CRC_NM_epithelial_cells_COUNT.csv"
output <- "results/2018-04-04_hcc_processing"

# Read data ---------------------------------------------------------------

tumor <- fread(here(tumorFile),
               sep = ",",
               check.names = TRUE,
               data.table = FALSE)

row.names(tumor) <- tumor$V1
tumor$V1 <- NULL


normal <- fread(here(normalFile),
                sep = ",",
                check.names = TRUE,
                data.table = FALSE)


row.names(normal) <- normal$V1
normal$V1 <- NULL


if(!all(row.names(tumor) == row.names(normal))){
  stop("Gene names are not ordered in both datasets")
}


# Bind datasets -----------------------------------------------------------

expData <- round(cbind(normal, tumor)) -> expDataCounts
metadata <- data.frame(status = c(rep("normal", ncol(normal)), 
                                  rep("tumor", ncol(tumor))), 
                       row.names = colnames(expData))


# Get metadata ------------------------------------------------------------


getField <- function(i){
  metadata %>% 
    row.names() %>% 
    str_split("__") %>%
    lapply("[", i) %>% 
    unlist()
}


fields <- lapply(1:2, getField)
metadata$id <- fields[[1]]
metadata$cellType <- fields[[2]]

cat("# Raw summary\n", sep = "\n", file = here(output, "data_summary.txt"))

metadata %>% 
  group_by(status, cellType) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  knitr::kable() %>% 
  cat(sep = "\n", file = here(output, "data_summary.txt"), append = TRUE)

cat("\nNumber of genes:", nrow(expData), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nNumber of cells:", ncol(expData), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  


# Remove cells ------------------------------------------------------------

# Remove cells by library size
libSizes <- colSums(expData)

filterCells <- function(filterParam){
  cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
  cellsToRemove
}


geneNames <- row.names(expData)

mtID <- grepl("MT-", geneNames, ignore.case = TRUE)
rbID <- grepl("RPL|RPS", geneNames, ignore.case = TRUE)

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

## Calculate CPM

cpm  <- apply(expData, 2, function(x) (x/sum(x))*1000000)

genesToKeep <- apply(cpm, 1, max) > 5 
cpm <- cpm[genesToKeep,]


cat("\n# Gene removal", "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  
cat("\nNumber of genes left:", nrow(cpm), "\n", sep = " ", file = here(output, "data_summary.txt"), append = TRUE)  


keepCellType <- rownames(metadata[metadata$cellType == "stemTA", ])
cpm <- cpm[, colnames(cpm) %in% keepCellType]

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
        file = here(output, "hcc_processed_counts.RDS"))

saveRDS(cpm, 
        file = here(output, "hcc_processed_cpm.RDS"))

saveRDS(newMetadata, 
        file = here(output, "hcc_processed_metadata.RDS"))


options(width = 120)
devtools::session_info() %>% 
  capture.output() %>% 
  cat("\nSession Info\n", .,  sep = "\n", file = here(output, "data_summary.txt"), append = TRUE)