# Load libraries ----------------------------------------------------------
library("tidyverse")
library("scPred")
library("here")
library("cowplot")
library("Seurat")

output <- file.path("results", "2018-05-19_pancreas_prediction")

# Read datasets -----------------------------------------------------------

# Read training
training <- readRDS(here("results", "2018-05-12_pancreas_processing", 
                         "pancreas_cpm_train.RDS"))
training_metadata <- readRDS(here("results", "2018-05-12_pancreas_processing", 
                                  "pancreas_metadata_train.RDS"))

# Read testing
test <- readRDS(here("results", "2018-05-12_pancreas_processing",
                     "baron_cpm_test.RDS"))
test_metadata <- readRDS(here("results", "2018-05-12_pancreas_processing",
                              "baron_metadata_test.RDS"))

# Integrate training datasets ---------------------------------------------

muraro <- CreateSeuratObject(raw.data = training$muraro, 
                             meta.data = training_metadata$muraro)
muraro <- NormalizeData(muraro)
muraro <- FindVariableGenes(muraro, do.plot = F, display.progress = F)
muraro@meta.data$batch <- "muraro"
muraro@meta.data$cell_type1 <- muraro@meta.data$x.cell_type.i.

segerstolpe <- CreateSeuratObject(raw.data = training$segerstolpe,
                                  meta.data = training_metadata$segerstolpe)
segerstolpe <- NormalizeData(segerstolpe)
segerstolpe <- FindVariableGenes(segerstolpe, do.plot = F, display.progress = F)
segerstolpe@meta.data$batch <- "segerstolpe"
segerstolpe@meta.data$cell_type1 <- segerstolpe@meta.data$x.cell_type.i.


xin <- CreateSeuratObject(raw.data = training$xin,
                          meta.data = training_metadata$xin)
xin <- NormalizeData(xin)
xin <- FindVariableGenes(xin, do.plot = F, display.progress = F)
xin@meta.data$batch <- "xin"
xin@meta.data$cell_type1 <- xin@meta.data$x.cell_type.i.

baron <- CreateSeuratObject(raw.data = test,
                          meta.data = test_metadata)
baron <- NormalizeData(baron)
baron <- FindVariableGenes(baron, do.plot = F, display.progress = F)
baron@meta.data$batch <- "baron"


# Determine genes to use for CCA, must be highly variable in at least 2 datasets
ob.list <- list(muraro = muraro, segerstolpe = segerstolpe, xin = xin)

assignData <- function(x){
  x@scale.data <- x@raw.data
  x
}

ob.list %>% 
  lapply(assignData) -> ob.list

genes.use <- c()
for (i in 1:length(ob.list)) {
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}


# Run multi-set CCA
# pancreas.integrated <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 15)
pancreas.integrated <- MergeSeurat(muraro, segerstolpe)
pancreas.integrated <- MergeSeurat(pancreas.integrated, xin)
pancreas.integrated <- ScaleData(pancreas.integrated)


# CC Selection
#MetageneBicorPlot(pancreas.integrated, grouping.var = "batch", dims.eval = 1:15)

# Run rare non-overlapping filtering
#pancreas.integrated <- CalcVarExpRatio(object = pancreas.integrated, reduction.type = "pca",
#                                       grouping.var = "batch", dims.use = 1:10)
#pancreas.integrated <- SubsetData(pancreas.integrated, subset.name = "var.ratio.pca",
#                                  accept.low = 0.5)

pancreas.integrated <- RunPCA(pancreas.integrated, pc.genes = genes.use, pcs.compute = 15)

pancreas <- as.data.frame(pancreas.integrated@dr$pca@cell.embeddings)
pancreas$cellType <- as.factor(pancreas.integrated@meta.data$cell_type1)
pancreas$batch <- as.factor(pancreas.integrated@meta.data$batch)

pancreas %>% 
  ggplot(aes(x = PC1, y = PC2, col = cellType)) +
  geom_point() -> train.pca.cellType

pancreas %>% 
  ggplot(aes(x = PC1, y = PC2, col = batch)) +
  geom_point() -> train.pca.batch


# Alignment
pancreas.integrated <- AlignSubspace(pancreas.integrated,
                                     reduction.type = "pca",
                                     grouping.var = "batch",
                                     dims.align = 1:10)


pancreas.aligned <- as.data.frame(pancreas.integrated@dr$pca.aligned@cell.embeddings)
pancreas.aligned$cellType <- as.factor(pancreas.integrated@meta.data$cell_type1)
pancreas.aligned$batch <- as.factor(pancreas.integrated@meta.data$batch)

pancreas.aligned %>% 
  ggplot(aes(x = APC1, y = APC2, col = cellType)) +
  geom_point() -> train.pca.alig.cellType

pancreas.aligned %>% 
  ggplot(aes(x = APC1, y = APC2, col = batch)) +
  geom_point() -> train.pca.alig.batch


all_train_p<- cowplot::plot_grid(train.pca.cellType, train.pca.batch, 
                   train.pca.alig.cellType, train.pca.alig.batch)

loadings <- pancreas.integrated@dr$pca@gene.loadings
genes <- pancreas.integrated@scale.data[rownames(pancreas.integrated@scale.data) %in% genes.use, ]

new.center <- rowMeans(genes)
new.scale <- apply(genes, 1, sd)

i <- match(rownames(loadings), names(new.center))
new.center <- new.center[i]

i <- match(rownames(loadings), names(new.scale))
new.scale <- new.scale[i]

share_genes <- intersect(rownames(baron@raw.data), rownames(loadings))

baron_sub <- baron@data[rownames(baron@data) %in% share_genes, ]
loadings_sub <- loadings[rownames(loadings) %in% share_genes,]
new.center <- new.center[rownames(loadings) %in% share_genes]
new.scale <- new.scale[rownames(loadings) %in% share_genes]
all(rownames(baron_sub) == rownames(loadings_sub))


#projection <- scale(log1p(t(baron_sub)), center = new.center, scale = new.scale) %*% loadings_sub
projection <- scale(t(baron_sub), center = new.center, scale = new.scale) %*% loadings_sub

#projection <- t(baron_sub) %*% loadings_sub



test_projection <- rbind(pancreas.integrated@dr$pca@cell.embeddings, projection)

pancreas.integrated_2 <- pancreas.integrated

pancreas.integrated_2@dr$pca@cell.embeddings <- as.matrix(test_projection)

trainInfo <- pancreas.integrated_2@meta.data[, c("batch", "cell_type1")]

testInfo <- data.frame(batch = "baron", cell_type1 = baron@meta.data$cell_type1, row.names = rownames(baron@meta.data))

pancreas.integrated_2@meta.data <- rbind(trainInfo, testInfo)

pancreas.integrated_2@meta.data$batch <- as.factor(pancreas.integrated_2@meta.data$batch)
pancreas.integrated_2@meta.data$cell_type1 <- as.factor(pancreas.integrated_2@meta.data$cell_type1)
all(rownames(pancreas.integrated_2@meta.data)  == rownames(pancreas.integrated_2@dr$pca@cell.embeddings))



pancreas.integrated_2@cell.names <- rownames(pancreas.integrated_2@meta.data)

pancreas <- as.data.frame(pancreas.integrated_2@dr$pca@cell.embeddings)
pancreas$cellType <- as.factor(pancreas.integrated_2@meta.data$cell_type1)
pancreas$batch <- as.factor(pancreas.integrated_2@meta.data$batch)

pancreas %>% 
  ggplot(aes(x = PC1, y = PC2, col = batch)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  ggtitle("No transformation") 

pancreas %>% 
  ggplot(aes(x = PC1, y = PC2, col = cellType)) +
  geom_point() -> train_test.pca.cellType






tmp <- MergeSeurat(pancreas.integrated, baron)
tmp@meta.data <- pancreas.integrated_2@meta.data
tmp@dr$pca <- pancreas.integrated_2@dr$pca


pancreas <- as.data.frame(tmp@dr$pca@cell.embeddings)
pancreas$cellType <- as.factor(tmp@meta.data$cell_type1)
pancreas$batch <- as.factor(tmp@meta.data$batch)


pancreas %>% 
  ggplot(aes(x = PC1, y = PC2, col = cellType)) +
  geom_point() -> train_test.pca.cellType

pancreas %>% 
  ggplot(aes(x = PC1, y = PC2, col = batch)) +
  geom_point() -> train_test.pca.batch



pancreas.integrated_2 <- AlignSubspace(tmp,
                                     reduction.type = "pca",
                                     grouping.var = "batch",
                                     dims.align = 1:10)





pancreas.aligned_2 <- as.data.frame(pancreas.integrated_2@dr$pca.aligned@cell.embeddings)
pancreas.aligned_2$cellType <- as.factor(pancreas.integrated_2@meta.data$cell_type1)
pancreas.aligned_2$batch <- as.factor(pancreas.integrated_2@meta.data$batch)

pancreas.aligned_2 %>% 
  ggplot(aes(x = APC1, y = APC2, col = cellType)) +
  geom_point(alpha = 0.4) -> p

pancreas.aligned_2 %>% 
  ggplot(aes(x = APC1, y = APC2, col = batch)) +
  geom_point(alpha = 0.2) -> p2


all(pancreas.integrated_2@dr$pca@cell.embeddings == test_projection)



i <- pancreas.integrated_2@meta.data$batch == "baron"

train_eigen <- pancreas.integrated_2@dr$pca.aligned@cell.embeddings[!i, ]
test_eigen <- pancreas.integrated_2@dr$pca.aligned@cell.embeddings[i, ]


train_metadata <- pancreas.integrated_2@meta.data[!i, ]
train_metadata$cell_type1 <- factor(train_metadata$cell_type1, levels = unique(train_metadata$cell_type1))

test_metadata <- pancreas.integrated_2@meta.data[i, ]


# Get informative principal components ------------------------------------
inf_pcs <- getInformativePCs(as.data.frame(train_eigen), classes = train_metadata$cell_type1, "cell_type1", sig = 0.05)



scpred <- list(prcomp = train_eigen, 
               metadata = train_metadata,
               pVar = "cell_type1",
               features = inf_pcs)



# Train model -------------------------------------------------------------
scpred <- trainModel(object = scpred, seed = 66, method = "svmRadial")

# Classify cells in new dataset -------------------------------------------
predictions <- eigenPredict(scpred, test_eigen, threshold = 0.9)


# Get accuracy results ----------------------------------------------------

getAccuracy <- function(predictions, metadata){
  res <- data.frame(pred = predictions$class, true = test_metadata$cell_type1, 
                    row.names = rownames(predictions))
  
  res %>% 
    mutate(true = if_else(as.character(true) %in% c("alpha", "beta", "delta", "gamma"), 
                          as.character(true), "other" )) %>% 
    mutate(label = if_else(pred == as.character(true), as.character(true),
                           if_else(pred ==  "Unassigned" & true == "other", 
                                   "other", "Misclassified"))) %>% 
    group_by(true, label) %>% 
    summarise(n = n()) %>% 
    mutate(accuracy = (n / sum(n))*100) %>% 
    filter(label != "Misclassified") 
}


getAccuracy(predictions, test_metadata)

getAccuracy(segerstolpe_predictions, segerstolpe_metadata)
getAccuracy(xin_predictions, xin_metadata)



getAccuracy <- function(predictions, metadata){
  res <- data.frame(pred = predictions$class, true = metadata$cell_type1, 
                    row.names = rownames(predictions))
  
  res %>% 
    mutate(label = if_else(pred == as.character(true), as.character(true),
                           if_else(pred ==  "Unassigned" & true == "other", 
                                   "other", "Misclassified")))
  
}



x <- getAccuracy(muraro_predictions, muraro_metadata)

table(x[as.character(x$true) != as.character(x$pred),"pred"])

table(x$label)



res <- data.frame(pred = muraro_predictions$class, true = muraro_metadata$cell_type1, 
                  row.names = rownames(muraro_predictions))


res %>% 
  mutate(true = if_else(as.character(true) %in% c("alpha", "beta", "delta", "gamma"), 
                        as.character(true), "other")) %>% 
  mutate(label = if_else(pred == as.character(true), as.character(true),
                         if_else(pred ==  "Unassigned" & true == "other", 
                                 "other", "Misclassified"))) -> x
group_by(true, label) %>% 
  summarise(n = n()) %>% 
  mutate(accuracy = (n / sum(n))*100)


table(x[x$true == "other","pred"])









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



res %>% 
  mutate(true = if_else(as.character(true) %in% c("alpha", "beta", "delta", "gamma"), 
                        as.character(true), "other" )) %>% 
  group_by(true, pred) %>% 
  summarise(n = n()) %>% 
  mutate(accuracy = (n / sum(n))*100)



# Get number of pcs and variance explained -------------------------------

pancreas@features %>% 
  Reduce(rbind, .) %>% 
  pull(PC) %>% 
  unique() %>% 
  as.character() -> unique_pcs
length(unique_pcs)
sum(pancreas@expVar[names(pancreas@expVar) %in% unique_pcs])
pancreas@features %>% 
  lapply(nrow)

# Get number of support vectors -------------------------------------------

pancreas@train %>% 
  lapply("[", "finalModel") %>% 
  unlist() %>% 
  lapply("slot", "nSV")


# Intersection of support cells
pancreas@train %>% 
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



#########

# Get gamma cells from reference
i <- reference_metadata$cellType == "gamma"
ref_gamma_metadata <- reference_metadata[i,,drop = FALSE]
ref_gamma <- reference[,i]
all(colnames(ref_gamma) == rownames(ref_gamma_metadata))

# Get gamma cells from test
i <- baron_metadata$cell_type1 == "gamma"
test_gamma_metadata <- as.data.frame(baron_metadata[i,])
test_gamma <- baron_cpm[,i]
all(colnames(test_gamma) == rownames(test_gamma_metadata))

ref_gamma <- ref_gamma[apply(ref_gamma, 1, var) != 0, ]
ref_gamma_cor <- cor(t(ref_gamma))
ref_gamma_cor <- abs(ref_gamma_cor)


test_gamma <- test_gamma[apply(test_gamma, 1, var) != 0, ]
test_gamma_cor <- cor(t(test_gamma))
test_gamma_cor <- abs(test_gamma_cor)


x <- apply(ref_gamma_cor, 1, 
           function(x){ 
             unlist(lapply(x, function(y){if(y >= 0.8){ 1}else{0}}))
           }
)

y <- apply(test_gamma_cor, 1, 
           function(x){ 
             unlist(lapply(x, function(y){if(y >= 0.8){ 1}else{0}}))
           }
)



<- x %*% diag(0,  nrow(x), ncol(x))

