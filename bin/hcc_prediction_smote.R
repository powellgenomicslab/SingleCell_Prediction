# Load libraries ----------------------------------------------------------

library("scPred")
library("here")
library("DMwR")

# Read data ---------------------------------------------------------------

expData <- here("results", "2018-04-04_hcc_processing", "hcc_processed_cpm.RDS")
expDataMetadata <- here("results", "2018-04-04_hcc_processing", "hcc_processed_metadata.RDS")

expData <- readRDS(expData)
expDataMetadata <- readRDS(expDataMetadata)


# Create results diretory -------------------------------------------------
output <- file.path("results", "2018-07-09_hcc_prediction_smote")

# Set up general variables ------------------------------------------------

phenoVar <- "status"
positiveClass <- "tumor"


# Get expression data and metadata ----------------------------------------

predictBoot <- function(seedPart){
  
  set.seed(seedPart)
  trainIndex <- createDataPartition(expDataMetadata[[phenoVar]], p = 0.75,  list = FALSE, times = 1)
  
  expTrain <- expData[, trainIndex]
  expTrainMeta <- expDataMetadata[trainIndex,]
  
  trainSplit <- as.data.frame(t(expTrain))
  trainSplit$status <- expTrainMeta$status
  

  trainSmote <- SMOTE(status ~ ., trainSplit, perc.over = 100, perc.under = 200)
  
  
  expTrain <- trainSmote[, colnames(trainSmote) != "status"]
  
  expTrainMeta <- data.frame(status = trainSmote[, colnames(trainSmote) == "status"],
                             row.names = rownames(expTrain))
  
  expTest  <- expData[,-trainIndex]
  expTestMeta <- expDataMetadata[-trainIndex,]
  expTest <- t(expTest)
  
  
  # Eigendecompse training matrix -------------------------------------------
  
  hcc <- eigenDecompose(as.matrix(expTrain))
  
  # Assign metadata to eigenPred object -------------------------------------
  expTrainMeta$status <- factor(expTrainMeta$status, levels = c("tumor", "normal"))
  metadata(hcc) <- expTrainMeta
  
  # Get informative principal components ------------------------------------
  hcc <- getInformativePCs(object = hcc, pVar = phenoVar)
  
  # Train prediction model --------------------------------------------------
  
  hcc <- trainModel(object = hcc, method = "svmRadial", seed = 66)
  #saveRDS(hcc, file = file.path(newDir, "scPred_object.RDS"))
  
  
  
  # Perform prediction in new dataset ---------------------------------------
  
  predictions <- eigenPredict(hcc, expTest, threshold = 0.9)
  predictions$true <- expTestMeta$status
  
  roc <- roc(predictor = predictions$probability, 
             response = predictions$true, 
             levels = c("tumor", "normal"))
  
  list(object = hcc, predictions = predictions, roc = roc)
  
}


set.seed(90)
i <- sample(1:10000, 50)
res <- lapply(i, predictBoot)

saveRDS(res, file = here(output, "bootstrap_replicates.RDS"))
# res <- readRDS(file = here(output, "bootstrap_replicates.RDS"))



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


library(precrec)

mm <- mmdata(scores = probs, labels = classes, modnames = "random", dsids = i)
mscurves <- evalmod(mm, raw_curves = TRUE)

mscurves %>% 
  auc() %>% 
  group_by(curvetypes) %>% 
  summarise(auc_mean = mean(aucs), sd = sd(aucs), n = n()) %>%
  group_by(curvetypes) %>% 
  summarise(auc_mean, error = qnorm(0.975) * sd / sqrt(n)) %>% 
  group_by(curvetypes) %>% 
  summarise(auc_mean,
            low = auc_mean - error, 
            high = auc_mean + error) %>% 
  as.data.frame()
  
# 0.9635517 0.9554553 0.9716481



autoplot(mscurves, show_cb = TRUE, curvetype = "ROC") +
  ggtitle("") -> roc_plot

autoplot(mscurves, show_cb = TRUE, curvetype = "PRC") +
  ggtitle("") -> prc_plot


metrics_plot <- cowplot::plot_grid(roc_plot, prc_plot)
ggsave(filename = here(output, "hcc_results.png"),
       plot = metrics_plot, 
       width = 6.5, height = 4,
       dpi = 350)


# Get metrics -------------------------------------------------------------

getMetrics <- function(x){
  conf_mat <- confusionMatrix(as.factor(ifelse(x$predictions$probability > 0.9, "tumor", "normal")), 
                              x$predictions$true, 
                              positive = "tumor")
  
  conf_mat$byClass
}

lapply(res, getMetrics) %>% 
  Reduce(rbind, .) %>% 
  as.tibble() -> res_metrics


res_metrics %>% 
  as.tibble() %>% 
  summarise_all(median)


plot(res_metrics$Sensitivity, res_metrics$Specificity)

ggplot(res_metrics, aes(x = Sensitivity)) +
  geom_histogram()


ggplot(res_metrics, aes(x = Specificity)) +
  geom_density()



res %>% 
  lapply("[[", "roc") %>% 
  lapply("[[", "auc") %>% 
  lapply(as.numeric) %>% 
  unlist() %>% 
  median()


res %>% 
  lapply("[[", "object") %>% 
  lapply(function(x){ x <- slot(x, "train"); x$tumor$results}) %>% 
  lapply(function(x){ i <- which.max(x$ROC); x[i, ]}) %>% 
  Reduce(rbind, .) %>% 
  select(-sigma, -C) %>% 
  summarise_all(mean)
  
  

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