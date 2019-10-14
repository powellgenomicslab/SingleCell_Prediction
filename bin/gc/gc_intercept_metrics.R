# Load libraries ----------------------------------------------------------

library("tidyverse")
library("here")
library("scPred")
library("magrittr")
library("dsLib")
library("MLmetrics")

# Set output --------------------------------------------------------------

output <- setOutput( "2019-07-06", "gc_intercept")

# Read data ---------------------------------------------------------------

## P5931

P5931_normal <- read.csv(here("data/2018-07-07_gastric_cancer/Garvan_collaboration/P5931/P5931_800B_n_epithelial.raw.data.csv"))
P5931_tumor  <- read.csv(here("data/2018-07-07_gastric_cancer/Garvan_collaboration/P5931/P5931_801B_tumor_epithelial.raw.data.csv"))

P5931_normal %<>% column_to_rownames("X")
P5931_tumor %<>% column_to_rownames("X")

shared_cells <- intersect(colnames(P5931_normal), colnames(P5931_tumor))

P5931_normal_i <- match(shared_cells, colnames(P5931_normal))
P5931_tumor_i <- match(shared_cells, colnames(P5931_tumor))

colnames(P5931_normal)[P5931_normal_i] <- paste0(shared_cells, "_normal")
colnames(P5931_tumor)[P5931_tumor_i] <- paste0(shared_cells, "_tumor")

## P6207

P6207_normal <- read.csv("data/2018-07-07_gastric_cancer/Garvan_collaboration/P6207/P6207_18522B_epithelial.raw.data.csv")
P6207_tumor  <- read.csv("data/2018-07-07_gastric_cancer/Garvan_collaboration/P6207/P6207_18523B_epithelial.raw.data.csv")

P6207_normal %<>% column_to_rownames("X")
P6207_tumor %<>% column_to_rownames("X")

# shared_cells <- intersect(colnames(P6207_normal), colnames(P6207_tumor))

getMetadata <- function(x, status){
  barcodes <- colnames(x)
  
  data.frame(barcodes, status) %>% 
    column_to_rownames("barcodes")
  
}


P5931_normal_info <- getMetadata(P5931_normal, "normal")
P5931_tumor_info  <- getMetadata(P5931_tumor, "tumor")

P6207_normal_info <- getMetadata(P6207_normal, "normal")
P6207_tumor_info  <- getMetadata(P6207_tumor, "tumor")


mergeDatasets <- function(d1, d2, by = c("r", "c")){
  
  by <- match.arg(by)
  
  if(by == "r"){
    d1_genes <- rownames(d1)
    d2_genes <- rownames(d2)
    
    shared_genes <- intersect(d1_genes, d2_genes)
    
    i <- match(shared_genes, d1_genes)
    d1 <- d1[i,]
    
    i <- match(shared_genes, d2_genes)
    d2 <- d2[i,]
    
    #all(rownames(d1) == rownames(d2))
    
    d <- t(cbind(d1, d2))
    d
  }else if(by == "c"){
    
    d1_genes <- colnames(d1)
    d2_genes <- colnames(d2)
    
    shared_genes <- intersect(d1_genes, d2_genes)
    
    i  <- match(shared_genes, d1_genes)
    d1 <- d1[,i]
    
    i  <- match(shared_genes, d2_genes)
    d2 <- d2[,i]
    
    #all(rownames(d1) == rownames(d2))
    
    d <- rbind(d1, d2)
    d
  }
}


P5931 <- mergeDatasets(d1 = P5931_normal, d2 = P5931_tumor)
P5931_metadata <- rbind(P5931_normal_info, P5931_tumor_info)
P5931_metadata$status <- factor(P5931_metadata$status, levels = c("tumor", "normal"))
all(rownames(P5931) == rownames(P5931_metadata))


P6207 <- mergeDatasets(d1 = P6207_normal, d2 = P6207_tumor)
P6207_metadata <- rbind(P6207_normal_info, P6207_tumor_info)
P6207_metadata$status <- factor(P6207_metadata$status, levels = c("tumor", "normal"))
all(rownames(P6207) == rownames(P6207_metadata))


gc <- mergeDatasets(P6207, P5931, by = "c")
gc_metadata <- rbind(P6207_metadata, P5931_metadata)

all(rownames(gc) == rownames(gc_metadata))
gc %<>% apply(1, function(x) (x/sum(x))*1000000)


# Create bootstrap prediction function -----

bootSubSamp <- function(seed, coef = 0){
  set.seed(seed)
  i <- createDataPartition(gc_metadata$status, list = FALSE)
  
  train_fold <- gc[, i]
  test_fold <- gc[, -i] 
  
  train_meta <- gc_metadata[i, , drop = FALSE]
  test_meta <- gc_metadata[-i, , drop = FALSE]
  
  scGastric <- eigenDecompose(train_fold, n = 25)
  metadata(scGastric) <- train_meta
  scGastric <- getFeatureSpace(scGastric, pVar = "status")
  scGastric <- trainModel(scGastric, seed = 66, returnData = TRUE, savePredictions = TRUE)
  
  model0 <- scGastric@train$tumor$finalModel
  model0@coef[[1]] <- rep(0, length(model0@coef[[1]]))
  scGastric@train$tumor$finalModel <- model0
  
  scGastric <- scPredict(scGastric, test_fold)
  scGastric@predMeta <-  test_meta
  
  scGastric  %>% 
    getPredictions() %>% 
    mutate(true = test_meta$status) %>% 
    mutate(predClass = ifelse(predClass == "unassigned" & tumor < 0.1, "normal", as.character(predClass))) %>% 
    mutate(prediction = ifelse(predClass == true, "correct", "incorrect")) %>% 
    group_by(true, prediction) %>% 
    summarise(n = n()) %>% 
    mutate(proportion = (n / sum(n))) -> props
  
  auroc <- AUC(scGastric@predictions$tumor, ifelse(test_meta$status == "tumor", 1, 0))
  auprc <- PRAUC(scGastric@predictions$tumor, ifelse(test_meta$status == "tumor", 1, 0))
  
  # scGastric@predictions %>% 
  #   mutate(predClass = ifelse(predClass == "unassigned", "normal", as.character(predClass))) %>% 
  #   pull(predClass) -> pred
  # 
  # debug(F1_Score)
  #   F1_Score(as.character(test_meta$status), pred, positive = "tumor")
  # 
  # 
  
  list(auroc = auroc, auprc = auprc, 
       props = props)
  
}


# Generate seeds
set.seed(66)
seeds <- sample(seq_len(10e4), 10)


res <- lapply(seeds, bootSubSamp)

saveRDS(res, here(output, "all_replicates.RDS"))

getCorrect <- function(x){
  
  xTrue <- split(x, x$true)
  
  getCorrectOnly  <- function(class){
    if(nrow(class) == 1){
      if(class$prediction == "incorrect"){
        class$proportion <- 1 - class$proportion
        class$prediction <- "correct"
        class$n <- 0
        class
      }else{
        class
      }
    }else{
      class %>% 
        filter(prediction == "correct")
    }
  }
  
  lapply(xTrue, getCorrectOnly) %>% 
    reduce(rbind) %>% 
    select(-prediction)
  
}

res %>% 
  map("props") %>% 
  lapply(getCorrect) %>% 
  mapply(function(r, rep){r$bootrep <- rep; r}, ., seeds, SIMPLIFY = FALSE) %>% 
  `names<-`(seeds) %>% 
  reduce(rbind)  -> resFormat


res %>% 
  map("auroc")

res %>% 
  map("auprc")


# Set coefficients to 1 ---------------------------------------------------


res1 <- lapply(seeds, bootSubSamp, coef = 1)

saveRDS(res1, here(output, "all_replicates_coef_1.RDS"))

res1 %>% 
  map("props") %>% 
  lapply(getCorrect) %>% 
  mapply(function(r, rep){r$bootrep <- rep; r}, ., seeds, SIMPLIFY = FALSE) %>% 
  `names<-`(seeds) %>% 
  reduce(rbind)  -> res1Format


res1 %>% 
  map("auroc") %>% 
  flatten_dbl() -> auroc

res1 %>% 
  map("auprc") %>% 
  flatten_dbl() -> auprc


mean(auroc) %>% round(3)
sd(auroc) %>% round(3)

mean(auprc) %>% round(3)
sd(auprc) %>% round(3)


res1Format %>% 
  group_by(true) %>% 
  summarize(mean = mean(proportion), sd = sd(proportion)) 

