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


# Merge gene expression data
reference <- merge(training$muraro, training$segerstolpe, by = 0)
reference <- column_to_rownames(reference, "Row.names")
reference <- merge(reference, training$xin, by = 0)
reference <- column_to_rownames(reference, "Row.names")
reference_test <- merge(reference, test, by = 0)
reference_test <- column_to_rownames(reference_test, "Row.names")


# Merge metadata

test_metadata_ss <- data.frame(x.cell_type.i. = test_metadata$cell_type1, 
                            row.names = rownames(test_metadata))

reference_test_metadata <- training_metadata
reference_test_metadata$baron <- test_metadata_ss
#reference_test_metadata <- reference_test_metadata[c("xin", "muraro", "segerstolpe", "baron")]

CreateSeuratObject(raw.data = training$)


# Batch effect correction -------------------------------------------------
muraro_ref <- reference_test[,colnames(reference_test) %in% row.names(reference_test_metadata$muraro)]
segerstolpe_ref <- reference_test[,colnames(reference_test) %in% row.names(reference_test_metadata$segerstolpe)]
xin_ref <- reference_test[,colnames(reference_test) %in% row.names(reference_test_metadata$xin)]
test_ref <- reference_test[,colnames(reference_test) %in% row.names(reference_test_metadata$baron)]

all(rownames(muraro_ref) == rownames(test_ref))

filterGenes <- function(x, quantile = 0.9){
  x_var <- apply(x, 1, var)
  cutoff <- quantile(x_var, quantile)
  i <- which(x_var > cutoff)
  rownames(x)[i]
}


vgenes <- lapply(list(muraro_ref, segerstolpe_ref, xin_ref, test_ref), filterGenes)
vgenes %>% 
 unlist() %>% 
  unique() -> vgenes

vgenes_i <- which(rownames(test_ref) %in% vgenes)


mnn_results <- mnnCorrect(as.matrix(log(xin_ref + 1)), 
                          as.matrix(log(muraro_ref + 1)), 
                          as.matrix(log(segerstolpe_ref + 1)), 
                          as.matrix(log(test_ref + 1)), 
                          subset.row = vgenes_i,
                          sigma = 1,
                          cos.norm.out = FALSE)

saveRDS(mnn_results, file = here(output, "mnn_results.RDS"))

training_correction <-mnn_results$corrected[1:3] %>% Reduce(cbind, .)
test_correction <- mnn_results$corrected[[4]]



mapply(function(dataset, batch){
  dataset$batch <- batch
  dataset
},reference_test_metadata, as.list(names(reference_test_metadata)), SIMPLIFY = FALSE) %>% 
  Reduce(rbind, .) %>% 
  set_colnames(c("cell_type1", "batch")) -> reference_test_metadata

all(rownames(reference_test_metadata) == colnames(reference_test))


# Build training
training_metadata <- reference_test_metadata[reference_test_metadata$batch != "baron", ]
all(rownames(training_metadata) == colnames(training_correction))


# Build test
all(rownames(test_metadata) == colnames(test_correction))



# Eigendecompose training gene expression data ----------------------------
pancreas <- eigenDecompose(t(training_correction), pseudo = FALSE)

training_metadata$batch <- as.factor(training_metadata$batch)
cell_types <- as.vector(unique(training_metadata$cell_type1))
training_metadata$cell_type1 <- factor(training_metadata$cell_type1, levels = cell_types)

scPred::metadata(pancreas) <- training_metadata


library(ggExtra)
plotEigen(pancreas, group = "batch", pc = c(1,3))
plotEigen(pancreas, group = "batch")


# Get informative principal components ------------------------------------
pancreas <- getInformativePCs(pancreas, pVar = "cell_type1", sig = 0.05)


# Train model -------------------------------------------------------------
pancreas <- trainModel(object = pancreas, seed = 66, method = "svmRadial")
# saveRDS(pancreas, here(output, "object_svmRadial.RDS"))
# pancreas <- readRDS(here(output, "object_svmRadial.RDS"))

# Classify cells in new dataset -------------------------------------------

#test_predictions <- eigenPredict(object = pancreas, newData = t(adj_test))

projection <- projectNewData(t(test_correction), pancreas, informative = TRUE)
projection$cell_type1 <- ifelse(test_metadata$cell_type1 %in% c("alpha", "beta", "delta", "gamma"),
                                as.character(test_metadata$cell_type1), "other") %>% 
  as.factor()
projection$batch <- "test"

plotEigen(pancreas, group = "cell_type1", geom = "density_2d", pc = c(1,2),  marginal = FALSE) +
  geom_density_2d(alpha = 0.3) +
  geom_point(data = projection, alpha = 0.3)


plotEigen(pancreas, group = "cell_type1",  marginal = TRUE)


predictions <- eigenPredict(pancreas, t(test_correction))


# Get accuracy results ----------------------------------------------------

getAccuracy <- function(predictions, metadata){
res <- data.frame(pred = predictions$class, true = metadata$cell_type1, 
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

