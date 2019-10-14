library("scPrediction")
library("tidyverse")
library("here")
library(DESeq)
library(BiocParallel)



# Read data ---------------------------------------------------------------

expDataDir <- "results/2017-10-14_colon_cancer_processing/coltumor_exp_data_processed_tpm.RDS"
expDataCountsDir <- "results/2017-10-14_colon_cancer_processing/coltumor_exp_data_processed_counts.RDS"
expSet <- readRDS(here(expDataDir))
expSetCounts <- readRDS(here(expDataCountsDir))



# Set prediction variables ------------------------------------------------

positiveClass <-  "tumor"
phenoVar <- "status"
negativeClass <- levels(expSet[[phenoVar]])[levels(expSet[[phenoVar]]) != positiveClass]




# Generate data partiion seeds --------------------------------------------

set.seed(100)
replicates <- sample(1:10000, 50)





# Perform predictions using DEGs ------------------------------------------

res <- predictClassDE(expData = expSet,
               expDataCounts = expSetCounts,
               phenoVar,
               positiveClass,
               seedPart = replicates[23],
               probPart = 0.5,
               log2FCFilter = 2,
               dispMethod = "per-condition",
               mlMethod = "svmPoly")


predictDE <- function(seed){
  predictClassDE(expData = expSet, expDataCounts = expSetCounts, phenoVar, positiveClass, seedPart = seed, 
                 probPart = 0.5, log2FCFilter = 2, pValFilter = 0.05, mlMethod = "svmPoly", dispMethod = "per-condition")
}

resDEG <- bplapply(replicates, predictDE, BPPARAM = MulticoreParam(workers = 2))


cat("Saving prediction results...\n")

# Save performance

outputDir <- "results/2018-03-13_hcc_DEGs_pred"
outputFileName <- paste0("performance_DEG_", positiveClass, "_vs_", negativeClass, ".RDS")
outputFile <- file.path(here(outputDir), outputFileName)
saveRDS(resDEG, file = outputFile)



# Session information


options(width = 120)
devtools::session_info()