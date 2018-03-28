# Set up command-line arguments -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
seedPart <- args[1]
positiveClass <- args[2]


# Load libraries ----------------------------------------------------------

library("DESeq")
library("here")
library("dplyr")
library("caret")


# Read data ---------------------------------------------------------------

pbmc <- readRDS(here("data/pbmc3k_filtered_gene_bc_matrices/pbmc3k_final_list.Rda"))

pbmc$meta.data %>% 
  mutate(cellType = if_else(cell.type == positiveClass, positiveClass, "other")) %>% 
  mutate(cellType = factor(cellType, levels = c(positiveClass, "other"))) -> expMetadata 
rownames(expMetadata) <- rownames(pbmc$meta.data)


# Create results diretory -------------------------------------------------
newDir <- here(file.path("results", "2018-03-27_pbmc_degs_feature-selection", paste0("degs_", positiveClass, "_boot-seed_", seedPart)))
dir.create(newDir)


# Set up general variables ------------------------------------------------

probPart <- 0.5
phenoVar <- "cellType"


# Get expression data and metadata ----------------------------------------

# Get count matrix
cellIndexes <- colnames(pbmc$raw.data) %in% colnames(pbmc$scale.data)
expData <- pbmc$raw.data[,cellIndexes] %>% Matrix::t() %>% as.matrix()


if(!all(rownames(expData) == rownames(expMetadata))){
  stop("Expression data and metadata are not ordered by cell id")
}


set.seed(seedPart)
trainIndex <- createDataPartition(expMetadata[[phenoVar]], p = probPart,  list = FALSE, times = 1)

expTrain <- expData[trainIndex, ]
expTrainMeta <- expMetadata[trainIndex, ]

expTest  <- expData[-trainIndex, ]
expTestMeta <- expMetadata[-trainIndex, ]

dataSummary <- capture.output(cat(sprintf("Number of genes: %i\nNumber of cells: %i\n", ncol(expData), nrow(expData))))


writeLines(file.path(newDir, "expData_summary.txt"), text = dataSummary, sep = "\n")



# Get differentially expressed genes --------------------------------------

degsRes <- getDE(expData = expTrain, 
                 expMetadata = expTrainMeta, 
                 pVar = phenoVar, 
                 positiveClass = positiveClass)

## Save diagnostic plot object and eigenPred summary
saveRDS(degsRes, file = file.path(newDir, "degsRes.RDS"))