# Import libraries --------------------------------------------------------

library("tidyverse")
library("here")
library("data.table")
library("stringr")

# Read data ---------------------------------------------------------------

input <- "/Users/j.alquicira/Documents/powell_lab/projects/singlecelloz_2018/data"


expData <- fread(file.path(input, "coordinates_gene_counts_flow_cytometry.txt"), 
                      sep = "\t", 
                      data.table = FALSE)


expData %>% 
  column_to_rownames("cell.Name") -> expData


expData[,c("DC1", "DC2", "DC3")] -> diff_mapp




fc_markers <- c("FSC.H",              "CD34",               "CD16",              
                "c.Kit",              "EPCR",               "Flk2",              
                "CD150",             "CD48",               "Lin",               
                "Sca1")


fc <- expData[,fc_markers]
expData[,c("DC1", "DC2", "DC3", fc_markers)] <- NULL


# Get metadata ------------------------------------------------------------

rownames(expData) %>% 
  str_split("_") %>% 
  lapply("[", 1) %>% 
  unlist() -> cellType

metadata <- data.frame(cellType = as.factor(cellType), row.names = rownames(expData))

expData <- t(expData) -> expDataCounts


output <- "results/2018-07-05_blood_processing"



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


newMetadata <- metadata[match(colnames(cpm), rownames(metadata)), , drop = FALSE]

if(!all(colnames(cpm) == rownames(newMetadata))){
  stop("Expression data and metadata cell ids do not match")
}


genesMatch <- match(rownames(cpm), rownames(expDataCounts))
cellsMatch <- match(colnames(cpm), colnames(expDataCounts))

expDataCounts <- expDataCounts[genesMatch, cellsMatch]

if(!(all(rownames(expDataCounts) == rownames(cpm)) & all(colnames(expDataCounts) == colnames(cpm)))){
  stop("Normalized matrix is not equivalent to count matrix")
}


fc <- t(fc)
i <- match(colnames(cpm), colnames(fc))
fc <- fc[,i]

all(colnames(fc) == colnames(cpm))



saveRDS(expDataCounts, 
        file = here(output, "blood_processed_counts.RDS"))

saveRDS(cpm, 
        file = here(output, "blood_processed_cpm.RDS"))

saveRDS(fc, 
        file = here(output, "blood_processed_fc.RDS"))

saveRDS(newMetadata, 
        file = here(output, "blood_processed_metadata.RDS"))
