# Script information ------------------------------------------------------

# title: Perform predictions of cancer cell status using degs
# author: José Alquicira Hernández
# date: 2019/03/05
# description:


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")
library("magrittr")

# Secondary
library("caret")
library("limma")
library("edgeR")
library("kernlab")
library("MLmetrics")

# Set output --------------------------------------------------------------


output_dir_name <- "gc_degs_prediction" # <------ Output directory

date <- "2019-03-05" # <------ Date

output <- file.path("results", paste(date, output_dir_name, sep = "_"))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# P5931

P5931_normal <- read.csv(here("data/2018-07-07_gastric_cancer/Garvan_collaboration/P5931/P5931_800B_n_epithelial.raw.data.csv"))
P5931_tumor  <- read.csv(here("data/2018-07-07_gastric_cancer/Garvan_collaboration/P5931/P5931_801B_tumor_epithelial.raw.data.csv"))

## Set row names
P5931_normal %<>% column_to_rownames("X") 
P5931_tumor %<>% column_to_rownames("X")

## Fix barcode repetition 
shared_barcodes <- intersect(colnames(P5931_normal), colnames(P5931_tumor))
P5931_normal_i <- match(shared_barcodes, colnames(P5931_normal))
P5931_tumor_i <- match(shared_barcodes, colnames(P5931_tumor))
colnames(P5931_normal)[P5931_normal_i] <- paste0(shared_barcodes, "_normal")
colnames(P5931_tumor)[P5931_tumor_i] <- paste0(shared_barcodes, "_tumor")

# P6207

P6207_normal <- read.csv("data/2018-07-07_gastric_cancer/Garvan_collaboration/P6207/P6207_18522B_epithelial.raw.data.csv")
P6207_tumor  <- read.csv("data/2018-07-07_gastric_cancer/Garvan_collaboration/P6207/P6207_18523B_epithelial.raw.data.csv")

## Set row names
P6207_normal %<>% column_to_rownames("X")
P6207_tumor %<>% column_to_rownames("X") 

# Create metadata ---------------------------------------------------------

getMetadata <- function(x, status, sample){
  barcodes <- colnames(x)
  data.frame(barcodes, status, sample) %>% 
    column_to_rownames("barcodes")
}

P5931_normal_info <- getMetadata(P5931_normal, "normal", "P5931")
P5931_tumor_info  <- getMetadata(P5931_tumor, "tumor", "P5931")

P6207_normal_info <- getMetadata(P6207_normal, "normal", "P6207")
P6207_tumor_info  <- getMetadata(P6207_tumor, "tumor", "P6207")


# Merge datasets

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

gc <- t(mergeDatasets(P6207, P5931, by = "c"))
gc_metadata <- rbind(P6207_metadata, P5931_metadata)

all(colnames(gc) == rownames(gc_metadata))


# ----


sums <- rowSums(gc)
i <- which(sums == 0)
if(length(i) > 0){
  gc2 <- gc[-i, ]
}

rpGenes <- grep(pattern = "^RP", x = rownames(gc2), value = TRUE)
propRp <- Matrix::colSums(gc2[rpGenes, ])/Matrix::colSums(gc2)
gc_metadata$propRp <- propRp 

gc_metadata %>% 
  ggplot(aes(x = propRp, fill = status)) +
  geom_histogram(color = "black") + 
  xlab("Ribosomal expression proportion") +
  theme_classic() +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(~sample)


status <- factor(gc_metadata$status, levels = c("normal", "tumor"))


# Differential expression analysis ----------------------------------------

dge <- DGEList(gc, group = status)
dge <- calcNormFactors(dge)
design <- model.matrix(~status + gc_metadata$propRp + gc_metadata$sample)
y <- new("EList")
y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
fit <- lmFit(y, design = design)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
tt <- topTable(fit, n = Inf, adjust.method = "BH")

tt %>% 
  rownames_to_column("gene") %>% 
  filter(abs(statustumor) > 2, adj.P.Val < 0.05) -> tt
degs <- tt$gene



####

# Set seed
set.seed(66)
seeds <- sample(seq_len(10e4), 10)


runPredictions <- function(seed){ 
  # Create partitions
  set.seed(seed)
  i <- createDataPartition(gc_metadata$status, list = FALSE)
  
  train_fold <- gc2[, i]
  test_fold <- gc2[, -i]
  
  train_meta <- gc_metadata[i, , drop = FALSE]
  test_meta <- gc_metadata[-i, , drop = FALSE]
  
  trCtrl <- trainControl(classProbs = TRUE,
                         method = "cv",
                         number = 10,
                         summaryFunction = twoClassSummary,
                         returnData = FALSE,
                         savePredictions = TRUE,
                         allowParallel = FALSE)
  
  
  fit <- train(t(train_fold[degs,]), 
               train_meta$status,
               method = "svmRadial",
               metric = "ROC",
               preProcess = c("center", "scale"),
               trControl = trCtrl)
  
  
  preds <- predict(fit, newdata = t(test_fold[degs,]))
  prob <- predict(fit, newdata = t(test_fold[degs,]), type = "prob")[["tumor"]]
  
  f1_score <-  F1_Score(preds, test_meta$status)
  
  
  table(x = test_meta$status, y= preds) %>% 
    as.data.frame() %>% 
    spread(key = "y", value = "Freq", fill = 0) %>%
    as.data.frame() %>% 
    column_to_rownames("x") -> counts
  
  
  labels <- rownames(counts)
  
  counts %>% 
    mapply(function(x,d){x/d}, ., colSums(.)) %>% 
    `rownames<-`(labels) -> props
  
  
  auroc <- AUC(prob, ifelse(test_meta$status == "tumor", 1, 0))
  auprc <- PRAUC(prob, ifelse(test_meta$status == "tumor", 1, 0))
  
  list(auroc = auroc, auprc = auprc, f1_score = f1_score,
       props = diag(props))
  
}

res <- lapply(seeds, runPredictions)


res %>% 
  map("props") %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  rename(sensitivity = tumor, specificity = normal) %>% 
  mutate(method = "degs") -> resDegs

res %>% 
  lapply("[", c("auroc", "auprc", "f1_score")) %>% 
  map(reduce, c) %>% 
  reduce(rbind) %>% 
  `colnames<-`(c("auroc", "auprc", "f1_score")) %>% 
  cbind(resDegs) %>% 
  `rownames<-`(1:10) -> resDegs


saveRDS(resDegs, here(output, "metrics.RDS"))



# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))