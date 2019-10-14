# Load libraries ----------------------------------------------------------

library("tidyverse")
library("here")
library("scPred")
library("magrittr")
library("BiocParallel")
library("dsLib")
library("MLmetrics")

# Set output --------------------------------------------------------------

output <- setOutput( "2019-07-06", "gc_sample_size")


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


# Simulate sequencing depth -----------------------------------------------

bootProp <- function(seed_part, prop){
  
  set.seed(seed_part)
  i <- createDataPartition(gc_metadata$status, list = FALSE)
  
  train_fold <- gc[, i]
  test_fold <- gc[, -i]
  
  
  nCells <- ncol(train_fold)
  set.seed(66)
  j <- sample(seq_len(nCells), prop)
  train_fold <- train_fold[, j]
  
  
  train_meta <- gc_metadata[i, , drop = FALSE]
  test_meta <- gc_metadata[-i, , drop = FALSE]
  train_meta <- train_meta[j, , drop = FALSE]
  
  scGastric <- eigenDecompose(train_fold, n = 25)
  metadata(scGastric) <- train_meta
  
  
  scGastric <- getFeatureSpace(scGastric, pVar = "status")
  
  scGastric <- trainModel(scGastric, seed = 66, returnData = TRUE, savePredictions = TRUE)
  
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
  
  scGastric  %>% 
    getPredictions() %>% 
    mutate(true = test_meta$status) %>% 
    mutate(predClass = ifelse(predClass == "unassigned" & tumor < 0.1, "normal", as.character(predClass))) %>% 
    pull(predClass) %>% 
    F1_Score(test_meta$status, .) -> f1_score
  
  list(auroc = auroc, auprc = auprc, f1_score = f1_score,
       props = props)
}

set.seed(66)
seed_part <- sample(seq_len(10e4), 10)

props <- seq(100, 900, 100)

runBoot <- function(i){
  lapply(seed_part, bootProp, props[i])
}




res <- bplapply(seq_len(length(props)), runBoot, BPPARAM = MulticoreParam(workers = 3))
saveRDS(res, here(output, "all_replicates.RDS"))
# res <- readRDS(here(output, "all_replicates.RDS"))


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

names(res) <- props

res %>% 
  lapply(map, "props") %>% 
  lapply(function(x) reduce(x, rbind)) %>% 
  lapply(getCorrect) %>% 
  mapply(function(x, label){ x$prop <- label; x}, .,  props, SIMPLIFY = FALSE) %>% 
  reduce(rbind) %>% 
  mutate(prop = as.factor(prop)) %>% 
  mutate(Metric = factor(true, levels = c("tumor", "normal"), 
                         labels = c("Sensitivity", "Specificity"))) -> resFormat


res %>% 
  lapply(map, "auroc") %>% 
  lapply(function(x) reduce(x, rbind)) %>% 
  map(as.data.frame) %>% 
  mapply(function(x, label){ x$prop <- label; x}, ., props, SIMPLIFY = FALSE) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  mutate(prop = factor(prop, levels = unique(prop))) %>% 
  set_names(c("auroc", "prop")) -> auroc

res %>% 
  lapply(map, "auprc") %>% 
  lapply(function(x) reduce(x, rbind)) %>% 
  map(as.data.frame) %>% 
  mapply(function(x, label){ x$prop <- label; x}, ., props, SIMPLIFY = FALSE) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  mutate(prop = factor(prop, levels = unique(prop))) %>% 
  set_names(c("auprc", "prop")) -> auprc

res %>% 
  lapply(map, "f1_score") %>% 
  lapply(function(x) reduce(x, rbind)) %>% 
  map(as.data.frame) %>% 
  mapply(function(x, label){ x$prop <- label; x}, ., props, SIMPLIFY = FALSE) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  mutate(prop = factor(prop, levels = unique(prop))) %>% 
  set_names(c("f1_score", "prop")) -> f1_score

metrics_long <- pivot_longer(metrics, 
                             names_to = "Metric", 
                             values_to = "proportion", 
                             cols = c("auroc", "auprc", "f1_score"))

resFormat %>% 
  ungroup() %>% 
  select(prop, Metric, proportion) %>% 
  rbind(metrics_long) -> metrics_long

metrics <- cbind(auroc, auprc %>% select(-prop), f1_score %>% select(-prop))

metrics_long %>% 
  filter(Metric %in% c("Sensitivity", "Specificity", "f1_score")) %>% 
  mutate(Metric = if_else(Metric == "f1_score", "F1 Score", as.character(Metric))) %>% 
  mutate(Metric = factor(Metric, c("Sensitivity", "Specificity", "F1 Score"))) %>% 
ggplot(aes(x = prop, y = proportion, fill = Metric)) +
  geom_boxplot(alpha = 0.7) +
  xlab("Number of cells in training dataset") +
  ylab("Value") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() -> p 



metrics %>% 
  gather(key = "Metric", value = "value", c(1,3)) %>% 
ggplot(aes(x = prop, y = value, fill = Metric)) +
  geom_boxplot() +
  xlab("Number of cells in training dataset") +
  ylab("Value") +
  scale_fill_manual(labels = c("AUROC", "AUPRC"), values = c("#985277", "#FF8C61", "#F0C808")) +
  theme_bw() -> p2 

p3 <- cowplot::plot_grid(p, p2, ncol = 1)
ggsave(filename = here(output, "sample_size.png"), p3, width = 6, height = 6)

metrics_long %>% 
  group_by(Metric, prop) %>% 
  summarize(mean = mean(proportion), sd = sd(proportion)) %>% 
  mutate(mean = round(mean, 3), sd = round(sd, 3)) %>% 
  as.data.frame()
  write_csv(path = here(output, "results.csv"))

metrics %>% 
  gather(key = "Metric", value = "value", c(1,3)) %>% 
  group_by(Metric, prop) %>% 
  summarize(mean = mean(value), sd = sd(value)) %>% 
  mutate(mean = round(mean, 3), sd = round(sd, 3)) %>% 
  as.data.frame() %>% 
  write_csv(path = here(output, "results_metrics.csv"))
