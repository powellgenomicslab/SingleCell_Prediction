# Set up command-line arguments -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
seedPart <- args[1]
positiveClass <- args[2]


# Load libraries ----------------------------------------------------------

library("scPred")
library("here")


# Set input/output variables ----------------------------------------------

expDataDir <- file.path("results", "2018-04-06_hbca_processing", "hbca_processed_cpm.RDS")
expDataMetaDir <- file.path("results", "2018-04-06_hbca_processing", "hbca_processed_metadata.RDS")


# Read data ---------------------------------------------------------------

expData <- readRDS(here(expDataDir))
expDataMeta <- readRDS(here(expDataMetaDir)) 

expDataMeta %>% 
  mutate(cellType = if_else(group == positiveClass, positiveClass, "other")) %>% 
  mutate(cellType = factor(cellType, levels = c(positiveClass, "other"))) -> tmp

rownames(tmp) <- rownames(expDataMeta)
expDataMeta <- tmp

# Create results diretory -------------------------------------------------
newDir <- here("results", "2018-04-06_hbca_scPred_feature-selection", paste0("scPred_", positiveClass, "_seed_", seedPart))
dir.create(newDir)


# Set up general variables ------------------------------------------------

probPart <- 0.5
phenoVar <- "cellType"


# Get expression data and metadata ----------------------------------------
expData <- t(expData) 

if(!all(rownames(expData) == rownames(expDataMeta))){
  stop("Expression data and metadata are not ordered by cell id")
}


set.seed(seedPart)
trainIndex <- createDataPartition(expDataMeta[[phenoVar]], p = probPart,  list = FALSE, times = 1)

expTrain <- expData[trainIndex, ]
expTrainMeta <- expDataMeta[trainIndex, ]

expTest  <- expData[-trainIndex, ]
expTestMeta <- expDataMeta[-trainIndex, ]


# Eigendecompse training matrix -------------------------------------------

## Remove zero-variance genes

i <- which(apply(expTrain, 2, var) == 0)

if(length(i > 0)){
  expTrain <- expTrain[,-i]
}

cat("Matrix decomposition...\n")
expTrainEigenDec <- eigenDecompose(expTrain)
cat("Done...\n")

# Assign metadata to eigenPred object -------------------------------------

metadata(expTrainEigenDec) <- expTrainMeta

# Get informative principal components ------------------------------------

expTrainEigenDec <- getInformativePCs(object = expTrainEigenDec, pVar = phenoVar)

## Save diagnostic plot object and eigenPred summary
saveRDS(expTrainEigenDec, file = file.path(newDir, "eigenPred_object.RDS"))
writeLines(file.path(newDir, "expData_summary.txt"), text = capture.output(expTrainEigenDec), sep = "\n")
