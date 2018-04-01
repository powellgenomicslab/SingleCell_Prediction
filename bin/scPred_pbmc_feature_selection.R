# Set up command-line arguments -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
seedPart <- args[1]
positiveClass <- args[2]


# Load libraries ----------------------------------------------------------

library("scPred")
library("here")


# Read data ---------------------------------------------------------------

pbmc <- readRDS(here("data/pbmc3k_filtered_gene_bc_matrices/pbmc3k_final_list.Rda"))

pbmc$meta.data %>% 
  mutate(cellType = if_else(cell.type == positiveClass, positiveClass, "other")) %>% 
  mutate(cellType = factor(cellType, levels = c(positiveClass, "other"))) -> expMetadata 
rownames(expMetadata) <- rownames(pbmc$meta.data)


# Create results diretory -------------------------------------------------
newDir <- here(file.path("results", "2018-03-27_pbmc_scPred_feature-selection", paste0("scPred_", positiveClass, "_boot-seed_", seedPart)))
dir.create(newDir)


# Set up general variables ------------------------------------------------

probPart <- 0.5
phenoVar <- "cellType"


# Get expression data and metadata ----------------------------------------
expData <- pbmc$data %>% Matrix::t() %>% as.matrix()

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


# Eigendecompse training matrix -------------------------------------------

## Remove zero-variance genes

i <- which(apply(expTrain, 2, var) == 0)

if(length(i > 0)){
  expTrain <- expTrain[,-i]
}

cat("Matrix decomposition...\n")
expTrainEigenDec <- eigenDecompose(expTrain, pseudo = TRUE)
cat("Done...\n")

# Assign metadata to eigenPred object -------------------------------------

metadata(expTrainEigenDec) <- expTrainMeta

# Get informative principal components ------------------------------------

expTrainEigenDec <- getInformativePCs(object = expTrainEigenDec, pVar = phenoVar)

## Save diagnostic plot object and eigenPred summary
saveRDS(expTrainEigenDec, file = file.path(newDir, "eigenPred_object.RDS"))