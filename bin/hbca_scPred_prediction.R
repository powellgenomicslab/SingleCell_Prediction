# Set up command-line arguments -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
seedPart <- args[1]
positiveClass <- args[2]
mlMethod <- args[3]



# Load libraries ----------------------------------------------------------

library("scPred")
library("here")

# Read data ---------------------------------------------------------------

dirData <- paste0("scPred_", positiveClass, "_seed_", seedPart)
eigenPred <- readRDS(here("results", "2018-04-06_hbca_scPred_feature-selection", dirData, "eigenPred_object.RDS"))


# Create results diretory -------------------------------------------------
newDir <- here("results", "2018-04-06_hbca_scPred_prediction", paste0("scPred_", positiveClass, "_seed_", seedPart, "_", mlMethod))
dir.create(newDir)


# Train prediction model --------------------------------------------------

trainedModel <- trainModel(object = eigenPred, all = TRUE, method = mlMethod, number = 10, positiveClass = positiveClass, seed = 66)
saveRDS(trainedModel, file = file.path(newDir, "trained_model.RDS"))



# Read data ---------------------------------------------------------------

expDataDir <- file.path("results", "2018-04-06_hbca_processing", "hbca_processed_cpm.RDS")
expDataMetaDir <- file.path("results", "2018-04-06_hbca_processing", "hbca_processed_metadata.RDS")

expData <- readRDS(here(expDataDir))
expDataMeta <- readRDS(here(expDataMetaDir)) 


expDataMeta %>% 
  mutate(cellType = if_else(group == positiveClass, positiveClass, "other")) %>% 
  mutate(cellType = factor(cellType, levels = c(positiveClass, "other"))) -> tmp

rownames(tmp) <- rownames(expDataMeta)
expDataMeta <- tmp

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

expTest  <- expData[-trainIndex, ]
expTestMeta <- expDataMeta[-trainIndex, ]

# Project prediction dataset into training principal components -----------

expTestProj <- projectNewData(newData = expTest, referenceData = eigenPred)


# Perform prediction in new dataset ---------------------------------------

predictions <- eigenPredict(eigenPred, expTestProj, trainedModel)  
saveRDS(predictions, file = file.path(newDir, "predictions.RDS"))


# Measure prediction performance ------------------------------------------

rocRes <-  roc(response = expTestMeta[[phenoVar]],
               predictor = predictions[[positiveClass]],
               levels = trainedModel$levels)
saveRDS(rocRes, file = file.path(newDir, "roc.RDS"))
