# Script information ------------------------------------------------------

# title: Get performance prediction of PBMCs using tree approach
# author: José Alquicira Hernández
# date: 2019-03-20
# description: None


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary
library("scPred")
library("Seurat")
library("BiocParallel")

# Set output --------------------------------------------------------------


output_dir_name <- "pbmc_prediction" # <------ Output directory

date <- "2019-03-20" # <------ Date

output <- file.path("results", paste(date, output_dir_name, sep = "_"))

if(!dir.exists(output)){
  dir.create(here(output))
}


# Read data ---------------------------------------------------------------

# Read gene expression data
pbmc <- readRDS(here("results", "2019-03-13_pbmc_assign_layers", "pbmc.RDS"))

# Read predictions
res <- readRDS(here(output, "predictions.RDS"))



# Create test datasets ----------------------------------------------------


# Set bootstrap seeds
set.seed(66)
seed_part <- sample(seq_len(10e4), 10)

cat("Creating test datasets.....\n")

# Get test folds
createTestDatasets <- function(seed){
  set.seed(seed)
  i <- createDataPartition(seq_len(nrow(pbmc@meta.data)), times = 1, p = 0.5, list = FALSE)
  pbmc@meta.data[-i,]
}


multicoreParam <- MulticoreParam(workers = 3)
testFolds <- bplapply(seed_part, createTestDatasets, BPPARAM = multicoreParam)

names(testFolds) <- paste0("r", seed_part)

rm(pbmc)

# Get performance ---------------------------------------------------------



crossTab <- function(tab, true, pred, fill = 0, prop = TRUE, digits = 2){
  
  tab %>% 
    group_by(!!sym(pred), !!sym(true)) %>% 
    summarise(n = n()) %>% 
    spread(key = true, value = "n", fill = fill) %>% 
    as.data.frame() %>% 
    column_to_rownames(pred) -> x
  
  if(prop){
    row_names <- rownames(x)
    x <- mapply(function(x,d){x/d}, x, colSums(x))
    rownames(x) <- row_names
    x %>% 
      #round(digits) %>% 
      as.data.frame() -> x
    
  }
  x
}


getPerformance <- function(pred, true, ...){
  predTrue <- cbind(pred, true[,paste0("level", 1:3)])
  
  r1 <- crossTab(predTrue, "level1", "layer1", ...)

  predTrue %>% 
    filter(layer1 == "Lymphoid_cell") -> predTrue
  
  
  r2 <- crossTab(predTrue, "level2", "layer2", ...)[c("B_cell", "Natural_killer", "T_cell"), c("B_cell", "Natural_killer", "T_cell")]

  predTrue %>% 
    filter(layer2 == "T_cell") -> predTrue
  
  r3 <- crossTab(predTrue, "level3", "layer3", ...)[c("cytotoxic", "non_cytotoxic"),c("cytotoxic", "non_cytotoxic")]  
  
  
  x <- list(level1 = r1, level2 = r2, level = r3)
  lapply(x, function(d){diag(as.matrix(d))})
  
}


performance <- bpmapply(getPerformance, res, testFolds, digits = 3, SIMPLIFY = FALSE, BPPARAM = multicoreParam)
performance_counts <- bpmapply(getPerformance, res, testFolds, digits = 4, prop = FALSE, SIMPLIFY = FALSE, BPPARAM = multicoreParam)


getLevel <- function(i, x){
  x %>% 
    lapply("[[", i) %>% 
    reduce(cbind)
}

lapply(1:3, getLevel, performance) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  set_names(paste0("r", seed_part)) -> performance

lapply(1:3, getLevel, performance_counts) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  set_names(paste0("r", seed_part)) -> performance_counts


performance %>% rowMeans() %>% `*`(100)
predCounts <- performance_counts %>% rowMeans()

testFolds %>% 
  lapply("[[","level1") %>% 
  lapply(table) %>% 
  reduce(rbind) %>% 
  colMeans() -> level1

testFolds %>% 
  lapply("[[","level2") %>% 
  lapply(table) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  select(B_cell, Natural_killer, T_cell) %>% 
  colMeans() -> level2
  

testFolds %>% 
  lapply("[[","level3") %>% 
  lapply(table) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  select(cytotoxic, non_cytotoxic) %>% 
  colMeans() -> level3


true <- c(level1, level2, level3)

mapply(paste, format(round(predCounts), format="d", big.mark=","), format(round(true), format="d", big.mark=","), sep = "/")


nodes <- c("Blood_progenitor", "Myeloid_cell", "B_cell", "Natural_killer", "cytotoxic", "non_cytotoxic")

paste(sum(round(predCounts)[nodes]), sum(round(true)[nodes]), sep = "/")

# Accuracy
accuracy <- {sum(round(predCounts)[nodes]) / sum(round(true)[nodes])} %>% round(2)


accuracies <- round(predCounts)[nodes] / round(true)[nodes]




performance %>% 
  apply(1, quantile, prob = c(0.025, 0.975)) %>%
  round(4) %>% 
  t() %>% 
  as.data.frame() -> CI

CI <- paste(CI$`2.5%`, CI$`97.5%`, sep = " - ")


performance %>% 
  apply(1, mean) %>% 
  data.frame(Mean = ., CI = CI) %>% 
  xtable::xtable()


performance %>% 
  rownames_to_column("cell") %>% 
  filter(cell %in% nodes) %>% 
  column_to_rownames("cell") %>% 
  apply(2, mean) %>% 
  quantile(prob = c(0.025, 0.975))

45884/47328
