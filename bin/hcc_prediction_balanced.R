# Set up command-line arguments -------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
seedPart <- 66
mlMethod <- "svmRadial"


# Load libraries ----------------------------------------------------------

library("scPred")
library("here")

# Read data ---------------------------------------------------------------

expData <- here("results", "2018-04-04_hcc_processing", "hcc_processed_cpm.RDS")
expDataCounts <- here("results", "2018-04-04_hcc_processing", "hcc_processed_cpm.RDS")
expDataMetadata <- here("results", "2018-04-04_hcc_processing", "hcc_processed_metadata.RDS")

expData <- readRDS(expData)
expDataCounts <- readRDS(expDataCounts)
expDataMetadata <- readRDS(expDataMetadata)


# Create results diretory -------------------------------------------------
newDir <- here("results", 
               "2018-04-05_prediction_hcc", 
               paste0("scPred_seed_", seedPart, "_", mlMethod))

dir.create(newDir)


# Set up general variables ------------------------------------------------

probPart <- 0.5
phenoVar <- "status"
positiveClass <- "tumor"


# Get expression data and metadata ----------------------------------------

predictBoot <- function(seedPart){

set.seed(seedPart)
trainIndex <- createDataPartition(expDataMetadata[[phenoVar]], p = 0.75,  list = FALSE, times = 1)

expTrain <- expData[, trainIndex]
expTrainMeta <- expDataMetadata[trainIndex,]

i <- expTrainMeta$status == "normal"

normal_i <- which(i)
tumor_i <- which(!i)

set.seed(seedPart)
newTumor_i <- sample(tumor_i, size = length(normal_i))

trainIndex <- c(newTumor_i, normal_i)
expTrain <- expData[, trainIndex]
expTrainMeta <- expDataMetadata[trainIndex,]


expTrain <- t(expTrain)

expTest  <- expData[,-trainIndex]
expTestMeta <- expDataMetadata[-trainIndex,]
expTest <- t(expTest)


# Eigendecompse training matrix -------------------------------------------

hcc <- eigenDecompose(expTrain)

# Assign metadata to eigenPred object -------------------------------------
expTrainMeta$status <- factor(expTrainMeta$status, levels = c("tumor", "normal"))
metadata(hcc) <- expTrainMeta

# Get informative principal components ------------------------------------
hcc <- getInformativePCs(object = hcc, pVar = phenoVar)

# Train prediction model --------------------------------------------------

hcc <- trainModel(object = hcc, method = mlMethod, seed = 66)
#saveRDS(hcc, file = file.path(newDir, "scPred_object.RDS"))



# Perform prediction in new dataset ---------------------------------------

predictions <- eigenPredict(hcc, expTest, threshold = 0.5)
predictions$true <- expTestMeta$status

roc <- roc(predictor = predictions$probability, 
    response = predictions$true, 
    levels = c("tumor", "normal"))

list(object = hcc, predictions = predictions, roc = roc)

}


set.seed(90)
i <- sample(1:10000, 25)
res <- lapply(i, predictBoot)


res %>% 
  lapply("[[", "predictions") %>% 
  lapply("[", c("true", "probability")) -> class_probs
  

class_probs %>% 
  lapply("[[", "true") %>% 
  lapply(function(x) ifelse(x == "tumor", 1, 0)) -> classes


class_probs %>% 
  lapply("[[", "probability") -> probs


res %>% 
  lapply("[[", "roc") %>% 
  lapply("[[", "auc") %>% 
  lapply(as.numeric) %>% 
  unlist() %>% 
  mean()



mm <- mmdata(scores = probs, labels = classes, modnames = "random", dsids = i)
mscurves <- evalmod(mm, raw_curves = TRUE)
autoplot(mscurves, show_cb = TRUE)







res %>% 
  lapply("[[", "roc") %>% 
  lapply("[", c("sensitivities", "specificities")) -> ss_metrics



ss_metrics %>% 
  lapply("[[", "sensitivities") %>% 
  as.data.frame() %>% 
  set_colnames(paste0("rep_", i)) -> sensitivities

ss_metrics %>% 
  lapply("[[", "specificities") %>% 
  as.data.frame() %>% 
  set_colnames(paste0("rep_", i)) -> specificities



sensitivities <- apply(sensitivities, 1, median)
specificities <- apply(specificities, 1, median)



resROC <- data.frame(sensitivities = sensitivities, specificities = specificities)



ggplot(resROC, aes(x = 1 - specificities, y = sensitivities)) +
  geom_line(size = 0.7) + 
  xlab("Median false polsitive rate") +
  ylab("Median true positive rate") +
  theme_bw() +
  scale_color_manual(values = c("#D84438", "#7953CC")) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11 ),  
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        strip.background = element_rect(fill = "#EEEEEE"),
        legend.position = c(0.85, 0.25))



table(predictions[,c("class", "true")])


library(precrec)

sscurves <- evalmod(scores = predictions$probability, labels =predictions$true)
autoplot(sscurves, "PRC", show_cb = TRUE)

sscurves$rocs


plot(sscurves, "PRC")
plot(sscurves)


classPred <- predictions[,c("class")]
classPred <- factor(ifelse(classPred == "Unassigned", "normal", "tumor"), levels = c("tumor", "normal"))


confusionMatrix(reference = predictions[,c("true")], classPred, positive = "tumor")


predictions %>% 
  ggplot(aes(x = probability, fill = true)) +
  geom_histogram()



# Measure prediction performance ------------------------------------------

rocRes <-  roc(response = expTestMeta[[phenoVar]],
               predictor = predictions$probability,
               levels = levels(metadata(hcc)[[phenoVar]]))
saveRDS(rocRes, file = file.path(newDir, "roc.RDS"))