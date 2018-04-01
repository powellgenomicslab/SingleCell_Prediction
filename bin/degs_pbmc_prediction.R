# Set up command-line arguments -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
seedPart <- args[1]
positiveClass <- args[2]
mlMethod <- args[3]
positiveClassFormat <- gsub("\\+", "", positiveClass)


# Load libraries ----------------------------------------------------------

library("here")
library("dplyr")
library("caret")
library("pROC")
source(here("bin/degs_prediction.R"))

# Read data ---------------------------------------------------------------

dirData <- paste0("degs_", positiveClass, "_boot-seed_", seedPart)
features <- readRDS(here(file.path("results", "2018-03-27_pbmc_degs_feature-selection", dirData, "degsRes.RDS")))


# Create results diretory -------------------------------------------------
newDir <- here(file.path("results", "2018-03-27_pbmc_degs_prediction", paste0("degs_", positiveClass, "_boot-seed_", seedPart, "_", mlMethod)))
dir.create(newDir)



# Read data ---------------------------------------------------------------

pbmc <- readRDS(here("data/pbmc3k_filtered_gene_bc_matrices/pbmc3k_final_list.Rda"))

pbmc$meta.data %>% 
  mutate(cellType = if_else(cell.type == positiveClass, positiveClassFormat, "other")) %>% 
  mutate(cellType = factor(cellType, levels = c(positiveClassFormat, "other"))) -> expMetadata 
rownames(expMetadata) <- rownames(pbmc$meta.data)

# Set up general variables ------------------------------------------------

probPart <- 0.5
phenoVar <- "cellType"


# Get expression data and metadata ----------------------------------------
expData <- pbmc$data %>% Matrix::t() %>% as.matrix()
expData <- log2(expData + 1) 

if(!all(rownames(expData) == rownames(expMetadata))){
  stop("Expression data and metadata are not ordered by cell id")
}


set.seed(seedPart)
trainIndex <- createDataPartition(expMetadata[[phenoVar]], p = probPart,  list = FALSE, times = 1)

expTrain  <- expData[trainIndex, ]
expTrainMeta <- expMetadata[trainIndex, ]

expTest  <- expData[-trainIndex, ]
expTestMeta <- expMetadata[-trainIndex, ]

dataSummary <- capture.output(cat(sprintf("Number of genes: %i\nNumber of cells: %i\n", ncol(expData), nrow(expData))))


writeLines(file.path(newDir, "expData_summary.txt"), text = dataSummary, sep = "\n")


# Train model -------------------------------------------------------------

trainedModel <- trainDEGModel(expTrain, expMetadata = expTrainMeta, method = mlMethod, features = features, pVar = phenoVar, 
                              positiveClass = positiveClassFormat, seed = 66)
saveRDS(trainedModel, file = file.path(newDir, "trained_model.RDS"))

# Perform prediction in new dataset ---------------------------------------

predictions <- degPredict(features, expTest, trainedModel)
saveRDS(predictions, file = file.path(newDir, "predictions.RDS"))


rocRes <-  roc(response = expTestMeta[[phenoVar]],
               predictor = predictions[[positiveClassFormat]],
               levels = trainedModel$levels)
saveRDS(rocRes, file = file.path(newDir, "roc.RDS"))