# Load libraries ----------------------------------------------------------
library("tidyverse")
library("scPred")
library("here")
library("cowplot")

# Read datasets -----------------------------------------------------------

# Read training
reference <- readRDS(here("results", "2018-04-15_pancreas_processing", "pancreas_cpm.RDS"))
reference_metadata <- readRDS(here("results", "2018-04-15_pancreas_processing", "pancreas_metadata.RDS"))

# Read test
baron <- readRDS(here("data", "2018-04-15_pancreas_baron", "baron-human.rds"))
baron_counts <- counts(baron)
baron_cpm <- apply(baron_counts, 2, function(x) (x/sum(x))*1000000)
baron_metadata <- as.data.frame(colData(baron))
#rm(baron)

output <- file.path("results", "2018-04-15_pancreas_prediction")


# Get library information

colnames(baron_cpm) %>% 
  str_split("_") %>% 
  lapply("[", 2) %>% 
  unlist() %>% 
  str_remove(".final") %>% 
  as.factor() -> baron_metadata$lib


baron_metadata %>% 
  group_by(lib, cell_type1) %>% 
  summarise(n = n()) %>% 
  filter(cell_type1 %in% c("alpha", "beta", "delta", "gamma"))

baron_metadata %>% 
  mutate(cell_type1 = if_else(cell_type1 %in% c("alpha", "beta", "delta", "gamma"), as.character(cell_type1), "other")) %>% 
  mutate(cell_type1 = factor(cell_type1, levels = c("alpha", "beta", "delta", "gamma", "other")) ) %>% 
  group_by(cell_type1) %>% 
  summarise(n = n())

# 

baron_metadata$human <- as.factor(baron_metadata$human)

baron_metadata %>%
  rownames_to_column("id") %>%
  filter(human %in% c("1", "2", "4")) %>%
  mutate(human = factor(human, levels = c("1", "2", "4"))) %>%
  column_to_rownames("id") -> baron_metadata

baron_cpm <- baron_cpm[,match(rownames(baron_metadata), colnames(baron_cpm))]

all(colnames(baron_cpm) == rownames(baron_metadata))




# reference$dataset <- reference_metadata$dataset
# 
# 
# apply(reference, 1, function(gene){
#   batch <- reference_metadata$dataset
#   data <- data.frame(cellType, batch, pc)
#   fit <- lm(gene ~ batch, data = data)
#   fit$residuals
# }) -> reference_correction
# 
# rownames(reference_correction) <- colnames(reference)
# saveRDS(pancreas, here(output, "reference_lm_correction_noh3.RDS"))
# 

# Eigendecompose training gene expression data ----------------------------
pancreas <- eigenDecompose(t(reference_correction))
scPred::metadata(pancreas) <- reference_metadata
# saveRDS(pancreas, here(output, "object_noh3.RDS"))
# pancreas <- readRDS(here(output, "object_noh3.RDS"))


# Get informative principal components ------------------------------------
pancreas <- getInformativePCs(pancreas, pVar = "cellType", sig = 0.05)
#pancreas_batch <- getInformativePCs(pancreas, pVar = "dataset", sig = 0.05)

# pancreas_batch@features %>% 
#   lapply("[[", "PC") %>% 
#   unlist() %>% 
#   as.vector() %>% 
#   unique() -> pcs_batch

pancreas@features %>% 
  lapply("[[", "PC") %>% 
  unlist() %>% 
  as.vector() %>% 
  unique() -> pcs_cellType


#cofounded_pcs_i <- intersect(pcs_batch, pcs_cellType)


pca <- getPCA(pancreas)
#cofounded_pcs <- pca[,cofounded_pcs_i]
cofounded_pcs <- pca[,pcs_cellType]


lapply(cofounded_pcs,
       function(pc){
         batch <- pancreas@metadata[["dataset"]]
         data <- data.frame(batch, pc)
         fit <- lm(pc ~ batch, data = data)
         fit$residuals
       }) -> adj_pcs

plot(adj_pcs$PC1, adj_pcs$PC2, col = pancreas@metadata[["cellType"]])
plot(adj_pcs$PC1, adj_pcs$PC2, col = pancreas@metadata[["dataset"]])

plot(pca$PC1, pca$PC2, col = batch)


# pancreas_2 <- pancreas
pca[,pcs_cellType] <- adj_pcs

pancreas@prcomp$x <- pca
pancreas_2 <- getInformativePCs(pancreas, pVar = "cellType", sig = 0.01, correction = "bonferroni")

plotEigen(pancreas_2, group = "dataset",  pc = c(1,2))
plotEigen(pancreas_2, group = "cellType",  pc = c(1,2))


# pancreas_2 <- removeBatch(pancreas, batch = "dataset")

# Train model -------------------------------------------------------------
pancreas_2 <- trainModel(object = pancreas_2, seed = 66, method = "svmRadial")
# saveRDS(pancreas_2, here(output, "object_svmRadial_noh3.RDS"))
# pancreas <- readRDS(here(output, "object_svmRadial_noh3.RDS"))

# Classify cells in new dataset -------------------------------------------
baron_predictions <- eigenPredict(object = pancreas_2, newData = t(baron_cpm))
# saveRDS(baron_predictions, here(output, "predictions_kknn_noh3.RDS"))
# baron_predictions <- readRDS(here(output, "predictions_svmRadial.RDS"))


all(rownames(baron_predictions) == rownames(baron_metadata))




projection <- projectNewData(t(baron_cpm), pancreas_2, informative = TRUE)
projection$cellType <- ifelse(baron_metadata$cell_type1 %in% c("alpha", "beta", "delta", "gamma"),
                                as.character(baron_metadata$cell_type1), "other") %>% 
  as.factor()
projection$dataset <- "test"

plotEigen(pancreas_2, group = "cellType", geom = "density_2d", pc = c(2,3),  marginal = FALSE) +
  geom_density_2d(alpha = 0.3) +
  geom_point(data = projection, alpha = 0.3)


plotEigen(pancreas_2, group = "dataset",  pc = c(2,3),  marginal = FALSE)



projection$



res <- data.frame(pred = baron_predictions$class, true = baron_metadata$cell_type1, 
                  row.names = rownames(baron_predictions))


table(reference_metadata$cellType)

# Get number of pcs and variance explained -------------------------------

pancreas_2@features %>% 
  Reduce(rbind, .) %>% 
  pull(PC) %>% 
  unique() %>% 
  as.character() -> unique_pcs
length(unique_pcs)
sum(pancreas@expVar[names(pancreas@expVar) %in% unique_pcs])
pancreas_2@features %>% 
  lapply(nrow)

# Get accuracy results ----------------------------------------------------

res %>% 
  mutate(true = if_else(as.character(true) %in% c("alpha", "beta", "delta", "gamma"), 
                        as.character(true), "other" )) %>% 
  mutate(label = if_else(pred == as.character(true), as.character(true),
                         if_else(pred ==  "Unassigned" & true == "other", 
                                 "other", "Misclassified"))) %>% 
  group_by(true, label) %>% 
  summarise(n = n()) %>% 
  mutate(accuracy = (n / sum(n))*100) %>% 
  filter(label != "Misclassified") %>% # accuracy per group
  pull(n) %>% 
  sum() -> x


x / nrow(res) # overall accuracy




# Get number of support vectors -------------------------------------------

pancreas_2@train %>% 
  lapply("[", "finalModel") %>% 
  unlist() %>% 
  lapply("slot", "nSV")


# Intersection of support cells
pancreas_2@train %>% 
  lapply("[", "finalModel") %>% 
  unlist() %>% 
  lapply("slot", "SVindex") %>% 
  unlist() %>% 
  unique() %>% 
  length()







res$human <- as.factor(baron_metadata$human)

res$human %>% table()
res %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  group_by(true, human) %>% 
  summarise(n = n())


res %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  group_by(pred, human) %>% 
  summarise(n = n())


res$lib <- baron_metadata$lib

res %>% 
  mutate(label = if_else(pred == as.character(true), as.character(true), "Misclassified")) %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  group_by(true, label, human, lib) %>% 
  summarise(prediction = n()) %>% 
  ungroup() %>% 
  group_by(true, human, lib) %>% 
  mutate(total = sum(prediction)) %>% 
  mutate(accuracy = prediction/total) %>% 
  filter(label != "Misclassified") %>% 
  ungroup() %>% 
  select(-true) %>% 
  View()




x <- cbind(res, baron_predictions)


x %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  ggplot(aes(y = probability, x = true, fill = human)) +
  geom_boxplot(alpha = 0.3) + 
  geom_hline(yintercept = 0.70) +
  facet_wrap(~lib)

x %>% 
  filter(true %in% c("alpha", "beta", "delta", "gamma")) %>% 
  ggplot(aes(y = probability, x = true, fill = lib)) +
  geom_boxplot(alpha = 0.3) + 
  geom_hline(yintercept = 0.70) +
  facet_wrap(~human)



###### Supplementary material



pancreas_svmRadial <- readRDS(here(output, "object_svmRadial_noh3.RDS"))
pancreas_rf <- readRDS(here(output, "object_rf_noh3.RDS"))
pancreas_kknn <- readRDS(here(output, "object_kknn_noh3.RDS"))
pancreas_earth <- readRDS(here(output, "object_earth_noh3.RDS"))
pancreas_glmnet <- readRDS(here(output, "object_glmnet_noh3.RDS"))
pancreas_glm <- readRDS(here(output, "object_glm_noh3.RDS"))
pancreas_nb <- readRDS(here(output, "object_nb_noh3.RDS"))

models <- c("svmRadial", "rf", "earth", "glmnet", "glm", "kknn", "nb")
pmodels <- list(pancreas_svmRadial, pancreas_rf, pancreas_earth, pancreas_glmnet, pancreas_glm, pancreas_kknn, pancreas_nb)
names(pmodels) <- models

getTrainResults <- function(cellType){
  
lapply(pmodels, function(x) x@train[[cellType]]$results) %>%
  lapply(function(x) x[which.max(x$ROC),c("ROC", "Sens", "Spec", "ROCSD", "SensSD", "SpecSD")]) %>%
  Reduce(rbind, .) %>%
  mutate(model = models) %>%
  select(!!! rlang::syms(c("model","ROC", "Sens", "Spec", "ROCSD", "SensSD", "SpecSD"))) %>%
  arrange(-ROC)
  # knitr::kable(format = "latex", booktabs = T)
}



cellTypes <- list(alpha = "alpha", beta = "beta", delta = "delta", gamma = "gamma")
trainResults <- lapply(cellTypes, getTrainResults)

baron_predictions <- eigenPredict(object = pmodels$nb, newData = t(baron_cpm))
# saveRDS(baron_predictions, here(output, "predictions_nb_noh3.RDS"))

predictions_ <- saveRDS(baron_predictions, here(output, "predictions_nb_noh3.RDS"))


predictions_svmRadial <- readRDS(here(output, "predictions_svmRadial_noh3.RDS"))
predictions_rf <- readRDS(here(output, "predictions_rf_noh3.RDS"))
predictions_kknn <- readRDS(here(output, "predictions_kknn_noh3.RDS"))
predictions_earth <- readRDS(here(output, "predictions_earth_noh3.RDS"))
predictions_glmnet <- readRDS(here(output, "predictions_glmnet_noh3.RDS"))
predictions_glm <- readRDS(here(output, "predictions_glm_noh3.RDS"))
predictions_nb <- readRDS(here(output, "predictions_nb_noh3.RDS"))



models <- c("svmRadial", "rf", "earth", "glmnet", "glm", "kknn", "nb")
pmodels <- list(predictions_svmRadial, predictions_rf, predictions_earth, 
                predictions_glmnet, predictions_glm, predictions_kknn, predictions_nb)
names(pmodels) <- models



getAccuracy <- function(res){
  res <- data.frame(pred = res$class, true = baron_metadata$cell_type1, 
                    row.names = rownames(res))
  res %>% 
    mutate(true = if_else(as.character(true) %in% c("alpha", "beta", "delta", "gamma"), 
                          as.character(true), "other" )) %>% 
    mutate(label = if_else(pred == as.character(true), as.character(true),
                           if_else(pred ==  "Unassigned" & true == "other", 
                                   "other", "Misclassified"))) %>% 
    group_by(true, label) %>% 
    summarise(n = n()) %>% 
    mutate(accuracy = (n / sum(n))*100)
}


lapply(pmodels, getAccuracy) %>% 
  lapply(function(x) x[x$label == x$true,])

