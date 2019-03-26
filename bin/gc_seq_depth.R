# Load libraries ----------------------------------------------------------

library("tidyverse")
library("here")
library("scPred")
library("magrittr")
library("BiocParallel")

# Set output --------------------------------------------------------------


output_dir_name <- "gc_seq_depth" # <------ Output directory

date <- "2019-03-11" # <------ Date

output <- file.path("results", paste(date, output_dir_name, sep = "_"))

if(!dir.exists(output)){
  dir.create(here(output))
}

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


# Simulate sequencing depth -----------------------------------------------

# Create bootstrap prediction function
bootSubSamp <- function(seed, thr){
  set.seed(seed)
  i <- createDataPartition(gc_metadata$status, list = FALSE)
  
  train_fold <- gc[i, ]
  test_fold <- gc[-i, ] 
  
  train_meta <- gc_metadata[i, , drop = FALSE]
  test_meta <- gc_metadata[-i, , drop = FALSE]
  
  seqDepth <- rowSums(train_fold)
  cellMax <- max(seqDepth)
  scaleFactor <- thr / cellMax
  
  train_fold <- round(train_fold * scaleFactor)
  
  zeroDrop <- rowSums(train_fold) == 0
  
  if(any(zeroDrop)){
    train_fold <- train_fold[!zeroDrop,]
    train_meta <- train_meta[!zeroDrop, , drop = FALSE]
  }
  
  
  train_fold <- t(train_fold/rowSums(train_fold) * 1e6)
  test_fold <- t(test_fold/rowSums(test_fold) * 1e6)
  
  scGastric <- eigenDecompose(train_fold, n = 25)
  metadata(scGastric) <- train_meta
  scGastric <- getFeatureSpace(scGastric, pVar = "status")
  scGastric <- trainModel(scGastric, seed = 66, returnData = TRUE, savePredictions = TRUE)
  scGastric <- scPredict(scGastric, test_fold)
  scGastric@predMeta <-  test_meta
  
  scGastric  %>% 
    getPredictions() %>% 
    mutate(true = test_meta$status) %>% 
    mutate(predClass = ifelse(predClass == "unassigned" & tumor < 0.1, "normal", predClass)) %>% 
    mutate(prediction = ifelse(predClass == true, "correct", "incorrect")) %>% 
    group_by(true, prediction) %>% 
    summarise(n = n()) %>% 
    mutate(proportion = (n / sum(n)))
}

# Get series of sequencing depth proportions
thrs <- seq(5000, 40000, 5000)

# Generate seeds
set.seed(60)
seeds <- sample(seq_len(10e4), 10)

runBoot <- function(i){
  lapply(seeds, bootSubSamp, thrs[i])
}

multicoreParam <- MulticoreParam(workers = 3)
res <- bplapply(seq_len(length(thrs)), runBoot, BPPARAM = multicoreParam)
saveRDS(res, here(output, "all_replicates_seq-depth.RDS"))

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
  lapply(function(x) reduce(x, rbind)) %>% 
  lapply(getCorrect) %>% 
  mapply(function(x, label){ x$thr <- label; x}, ., thrs, SIMPLIFY = FALSE) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  filter(thr > 5000) %>% 
  mutate(thr = factor(thr, levels = unique(thr))) -> resFormat




ggplot(resFormat, aes(x = thr, y = proportion, fill = true)) +
  geom_boxplot() +
  xlab("Maximimum cell sequencing depth in training data") +
  ylab("Value") +
  labs(fill = "Metric") +
  scale_fill_manual(labels = c("Sensitivity", "Specificity"), values = c("red", "steelblue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p

ggsave(filename = here(output, "seq_depth.png"), p, width = 6, height = 3.5)

resFormat %>% 
  mutate(accuracy = proportion) %>% 
  group_by(thr, true) %>% 
  summarise(mean = mean(proportion), sd = sd(proportion)) %>% 
  write_csv(path = here(output, "results.csv"))






