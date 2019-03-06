# Load libraries ----------------------------------------------------------
library("tidyverse")
library("scPred")
library("here")
library("cowplot")
library("Seurat")

output <- file.path("results", "2019-03-03_pancreas_prediction_inverse")

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read datasets -----------------------------------------------------------

# Read training
test <- readRDS(here("results", "2018-05-12_pancreas_processing", 
                         "pancreas_cpm_test.RDS"))
test_metadata <- readRDS(here("results", "2018-05-12_pancreas_processing", 
                                  "pancreas_metadata_test.RDS"))

# Read testing
training <- readRDS(here("results", "2018-05-12_pancreas_processing",
                     "baron_cpm_train.RDS"))
training_metadata <- readRDS(here("results", "2018-05-12_pancreas_processing",
                              "baron_metadata_train.RDS"))


# Pre-process training dataset --------------------------------------------

training_metadata$cellType <- factor(training_metadata$x.cell_type.i., 
                                     levels = c("alpha", "beta", "delta", "gamma"))

training <- CreateSeuratObject(raw.data = training, 
                        meta.data = training_metadata) %>% 
  NormalizeData() %>% 
  FindVariableGenes(do.plot = FALSE, display.progress = FALSE) %>% 
  ScaleData() %>% 
  RunPCA()



# Get feature space -------------------------------------------------------

model <- getFeatureSpace(training, pVar = "cellType")


# Train model -------------------------------------------------------------

model <- trainModel(model)



# Predict cell types ------------------------------------------------------

createDataset <- function(x, meta){
  CreateSeuratObject(raw.data = x, 
                     meta.data = meta) %>% 
    NormalizeData()
}

test <- mapply(createDataset, test, test_metadata)

islets <- c("alpha", "beta", "delta", "gamma")


writeRes <- function(m, name){
  
  r <- crossTab(m, true = "cell_type1", prop = FALSE)
  i <- colnames(r) %in% islets
  orderClass <- c(islets, colnames(r)[!i])
  counts <- r[,orderClass]
  print(counts)
  write.table(counts, file = here(output, paste0(name, "_counts.txt")), 
              quote = FALSE, 
              sep = "\t", 
              row.names = TRUE, 
              col.names = TRUE)
  
  
  r <- crossTab(m, true = "cell_type1", prop = TRUE)
  i <- colnames(r) %in% islets
  orderClass <- c(islets, colnames(r)[!i])
  props <- r[,orderClass]
  print(props)
  write.table(props[,orderClass], file = here(output, paste0(name, "_props.txt")), 
              quote = FALSE, 
              sep = "\t", 
              row.names = TRUE, 
              col.names = TRUE)
}



alignment1 <- alignDataset(model, test$muraro)
model1 <- model
model1@projection <- alignment1
model1@predMeta <- test$muraro@meta.data
model1 <- scPredict(model1, useProj = TRUE, threshold = 0.9)


writeRes(model1, "muraro")


alignment2 <- alignDataset(model, test$segerstolpe)
model2 <- model
model2@projection <- alignment2
model2@predMeta <- test$segerstolpe@meta.data
model2 <- scPredict(model2, useProj = TRUE, threshold = 0.9)

writeRes(model2, "segerstolpe")


alignment3 <- alignDataset(model, test$xin)
model3 <- model
model3@projection <- alignment3
model3@predMeta <- test$xin@meta.data
model3 <- scPredict(model3, useProj = TRUE, threshold = 0.9)

writeRes(model3, "xin")


getPred <- function(x){ 
  x %>% getPredictions() %>% pull(predClass) -> pred
  x@predMeta$cell_type1 -> meta
  data.frame(pred, true = meta)
}

lapply(list(model1, model2, model3), getPred) %>% 
  reduce(rbind) %>% 
  mutate(true = if_else(true %in% islets, as.character(true), "other")) %>% 
  group_by(pred, true) %>%
  summarise(n = n()) %>%
  spread(key = true, value = "n", fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames("pred") -> counts

row_names <- rownames(counts)
props <- mapply(function(x,d){x/d}, counts, colSums(counts))
rownames(props) <- row_names
props %>%
  round(2) %>%
  as.data.frame() -> props

counts
props


write.table(counts, file = here(output, paste0("all", "_counts.txt")), 
            quote = FALSE, 
            sep = "\t", 
            row.names = TRUE, 
            col.names = TRUE)

write.table(props, file = here(output, paste0("all", "_props.txt")), 
            quote = FALSE, 
            sep = "\t", 
            row.names = TRUE, 
            col.names = TRUE)
