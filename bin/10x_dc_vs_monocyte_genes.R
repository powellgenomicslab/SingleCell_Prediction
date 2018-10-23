# Script information ------------------------------------------------------

# title:
# author: José Alquicira Hernández
# date:
# description:


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("data.table")
library("here")

# Secondary
library("caret")


# Set output --------------------------------------------------------------


output_dir_name <- "dc_vs_monocyte_gene_pred" # <------ Output directory

date <- format(Sys.Date(), format = "%Y-%m-%d_")
date <- "2018-10-22_" # <------ Fix date

output <- file.path("results", paste0(date, output_dir_name))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input    <- file.path("data", "2018-10-22_human_markers") # <------ Input directory
filename <- "Human_cell_markers.txt" # <------ Input file


# Read file

data <- fread(input = here(input, filename))
data[cellName3 %like% "Dendritic|dendritic" & 
       cellType == "Normal cell" & 
       tissueType2 %like% "blood|Blood",][["cellMarker"]] %>% 
  str_split(", ") %>% 
  reduce(c) %>% 
  table() %>% 
  sort(decreasing = TRUE) %>%
  #.[. > 1] %>% 
  names() -> dc


data[cellName3 %like% "Monocyte|monocyte" & 
       !(cellName3 %like% "derived|progenitor") &
       cellType == "Normal cell" & 
       tissueType2 %like% "blood|Blood",][["cellMarker"]] %>% 
  str_split(", ") %>% 
  reduce(c) %>% 
  table() %>% 
  sort(decreasing = TRUE) %>% 
  names -> monocytes


#####

# Read data ---------------------------------------------------------------

# Input

input <- file.path("results", "2018-10-12_10x_immune_myeloid_vs_lymphoid_prediction") # <------ Input directory

# Read files

filename <- "predictions.RDS" # <------ Input file
predictions <- readRDS(here(input, filename))


input <- file.path("results", "2018-10-12_10x_immune_tree") # <------ Input directory


filename <- "train_data.RDS" # <------ Input file
train_data <- readRDS(here(input, filename))


filename <- "pred_data.RDS" # <------ Input file
pred_data <- readRDS(here(input, filename))

# Subset cells

pred_data@meta.data$pred1 <- as.factor(predictions$predClass)

# train_data <- RunTSNE(train_data)
# TSNEPlot(train_data, group = "cellType1")


i <- rownames(train_data@meta.data)[train_data@meta.data$cellType2 %in% c("Lymphoid_cell", "Myeloid_cell")]
j <- rownames(pred_data@meta.data)[pred_data@meta.data$pred1 %in% c("Myeloid_cell")]


train_data <- SubsetData(train_data, cells.use = i, do.clean = TRUE)
train_data@meta.data$cellType2 <- factor(train_data@meta.data$cellType2)

pred_data <- SubsetData(pred_data, cells.use = j, do.clean = TRUE)
pred_data@meta.data$pred1 <- factor(pred_data@meta.data$pred1)



i <- rownames(train_data@meta.data)[train_data@meta.data$cellType1 %in% c("Dendritic_cell", "Monocyte")]
train_data <- SubsetData(train_data, cells.use = i, do.clean = TRUE)
train_data@meta.data$cellType1 <- factor(train_data@meta.data$cellType1)



# Process data ------------------------------------------------------------

# Train
## Subset marker genes and normalize data

genes.use <- which(row.names(train_data@raw.data) %in% unique(c(dc, monocytes)))
train_genes <- train_data@raw.data[genes.use, ]
train_genes <- CreateSeuratObject(train_genes)
train_genes <- AddMetaData(object = train_genes, metadata = train_data@meta.data[, c(4,6,7)]) # Add the idents to the meta.data slot

train_genes <- train_genes %>% 
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>% 
  ScaleData()


# Test
## Subset marker genes and normalize data

genes.use <- which(row.names(pred_data@raw.data) %in% unique(c(dc, monocytes)))
test_genes <- pred_data@raw.data[genes.use, ]
test_genes <- CreateSeuratObject(test_genes)
test_genes <- AddMetaData(object = test_genes, metadata = pred_data@meta.data[, c(4,6,7)]) # Add the idents to the meta.data slot

test_genes <- test_genes %>% 
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>% 
  ScaleData()

saveRDS(train_genes, file = here(output, "train_dataset.RDS"))
saveRDS(test_genes, file = here(output, "test_dataset.RDS"))


# Set training and test datasets

train_set <- t(as.matrix(train_genes@data))
test_set <- t(as.matrix(test_genes@data))


# Remove near-zero variance features

# nearZeroVar(train_set, saveMetrics = TRUE)
# nvz <- nearZeroVar(train_set)
# train_set <- train_set[, -nvz]

# Check correlated predictors

## > No correlated predictors
cor(train_set) 

# Check linear combinations

## > No linear combinations
findLinearCombos(train_set)

data.frame(scale(train_set), cellType = train_genes@meta.data$cellType1) %>% 
  ggplot(aes(CD14, FCN1, color = cellType)) +
  geom_point() +
  scale_color_manual(values = set_names(jcolors::jcolors(), NULL)) +
  theme_bw()



# pp <- preProcess(train_set, method = c("YeoJohnson", "nzv"))



fitControl <- trainControl(method = "cv",
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           number = 10)

train_set <- data.frame(train_set, cellType1 = train_data@meta.data$cellType1)

svmFit <- train(cellType1 ~ ., data = train_set, 
                trControl = fitControl,
                method = "svmRadial", 
                metric = "ROC",
                preProcess = c("center", "scale", "YeoJohnson", "nzv")
                )

new_names <- make.names(colnames(test_set))
colnames(test_set) <- new_names


preds <- predict(svmFit, newdata = test_set, type = "prob")

preds <- preds %>% 
  mutate(predClass = if_else(Dendritic_cell > 0.9, "DC", 
                             if_else(Monocyte > 0.9, "Monocyte", "Unassigned")))

preds$true <- factor(test_genes@meta.data$cellType1)


######

tidy_table <- function(x, left, top, fill = 0, prop = FALSE, digits = 2){
  x %>% 
    group_by_(left, top) %>% 
    summarise(n = n()) %>% 
    spread(key = top, value = "n", fill = fill) %>% 
    as.data.frame() %>% 
    column_to_rownames(left) -> x
  
  if(prop){
    row_names <- rownames(x)
    x <- mapply(function(x,d){x/d}, x, colSums(x))
    rownames(x) <- row_names
    x %>% 
      round(digits) %>% 
      as.data.frame() -> x
    
  }
  
  x
  
}

svmFit

tidy_table(preds, "predClass", "true", prop = TRUE) %>% 
  select(-B_cell, -Blood_progenitor, -T_cell, -Natural_killer)


saveRDS(preds, here(output, "predictions.RDS"))
saveRDS(svmFit, here(output, "pred_model.RDS"))


# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))