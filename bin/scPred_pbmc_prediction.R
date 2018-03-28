# Set up command-line arguments -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
seedPart <- args[1]
positiveClass <- args[2]
mlMethod <- args[3]
positiveClassFormat <- gsub("\\+", "", positiveClass)



# Load libraries ----------------------------------------------------------

library("scPred")
library("here")

# Read data ---------------------------------------------------------------

dirData <- paste0("scPred_", positiveClass, "_boot-seed_", seedPart)
eigenPred <- readRDS(here(file.path("results", "2018-03-27_pbmc_scPred_feature-selection", dirData, "eigenPred_object.RDS")))
eigenPred@metadata %>% 
  mutate(cellType = gsub("\\+", "", cellType)) %>% 
  mutate(cellType = factor(cellType, levels = c(positiveClassFormat, "other"))) -> newMetadata
rownames(newMetadata) <- rownames(eigenPred@metadata)

metadata(eigenPred) <- newMetadata

# Create results diretory -------------------------------------------------
newDir <- here(file.path("results", "2018-03-27_pbmc_scPred_prediction", paste0("scPred_", positiveClass, "_boot-seed_", seedPart, "_", mlMethod)))
dir.create(newDir)


# Train prediction model --------------------------------------------------

trainedModel <- trainModel(object = eigenPred, top = 10, method = mlMethod, number = 10, positiveClass = positiveClassFormat, seed = 66)
saveRDS(trainedModel, file = file.path(newDir, "trained_model.RDS"))



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
expData <- pbmc$scale.data %>% Matrix::t()

if(!all(rownames(expData) == rownames(expMetadata))){
  stop("Expression data and metadata are not ordered by cell id")
}


set.seed(seedPart)
trainIndex <- createDataPartition(expMetadata[[phenoVar]], p = probPart,  list = FALSE, times = 1)

expTest  <- expData[-trainIndex, ]
expTestMeta <- expMetadata[-trainIndex, ]

dataSummary <- capture.output(cat(sprintf("Number of genes: %i\nNumber of cells: %i\n", ncol(expData), nrow(expData))))


writeLines(file.path(newDir, "expData_summary.txt"), text = dataSummary, sep = "\n")


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
