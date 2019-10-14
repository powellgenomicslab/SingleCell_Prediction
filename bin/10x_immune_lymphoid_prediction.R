# Script information ------------------------------------------------------
# title: Classification of dendritic cells and monocytes
# author: José Alquicira Hernández
# date: 2018-10-12
# description: None


# Import libraries --------------------------------------------------------

# Primary

library("tidyverse")
library("here")

# Secondary
library("Seurat")
library("scPred")



# Set output --------------------------------------------------------------


output_dir_name <- "10x_t-cells_b-cells_nk-cells" # <------ Output directory
date <- "2018-10-13_"
output <- file.path("results", paste0(date, output_dir_name))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input <- file.path("results", "2018-10-12_10x_immune_myeloid_vs_lymphoid_prediction") # <------ Input directory

# Read files

filename <- "predictions.RDS" # <------ Input file
predictions <- readRDS(here(input, filename))


input <- file.path("results", "2018-10-12_10x_immune_tree") # <------ Input directory


filename <- "train_lymphoid.RDS" # <------ Input file
train_data <- readRDS(here(input, filename))




# Subset cells



# Process data ------------------------------------------------------------

# cat("\nProcessing training data (normalization, scaling, and PCA)...\n")
# 
train_data <- train_data %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000) %>%
  # FindVariableGenes(mean.function = ExpMean,
  #                   dispersion.function = LogVMR,
  #                   x.low.cutoff = 0.02,
  #                   x.high.cutoff = 5,
  #                   y.cutoff = 0.5) %>%
  ScaleData() %>%
  RunPCA(pc.genes = rownames(.@data),
         do.print = TRUE,
         pcs.print = 1:5,
         genes.print = 5)

PCAPlot(train_data, group = "cellType2")

# Get feature space -------------------------------------------------------

# cat("\nGetting feature space...\n")
sc_pred <- getFeatureSpace(train_data, pVar = "cellType2")
# rm(train_data)
# 
# cat("\nSaving scPred object...\n")
# saveRDS(sc_pred, file = here(output, "scpred.RDS"))
# 

# Train models ------------------------------------------------------------
#sc_pred <- readRDS(file = here(output, "scpred.RDS"))

# library(doParallel)
# cl <- makePSOCKcluster(5)
# registerDoParallel(cl)
# sc_pred <- trainModel(sc_pred)
# stopCluster(cl)

#saveRDS(sc_pred, file = here(output, "scpred_train.RDS"))
sc_pred <- readRDS(file = here(output, "scpred_train.RDS"))


pca <- as.data.frame(getPCA(sc_pred))

gene <- "CD8"
exp <- sc_pred@trainData[rownames(sc_pred@trainData) == gene,]
pca[[make.names(gene)]] <- exp
pca$class <- sc_pred@metadata$cellType1

pca %>% 
ggplot() +
  aes_string(x = "PC2", y = "PC3", color = make.names(gene)) +
  geom_point(size = 0.5) +
  geom_rug() +
  scale_color_continuous(low = "steelblue", high = "red") +
  facet_wrap(~class) +
  theme_bw()


x <- readRDS("/Users/jose/Downloads/all_pure_select_11types.RDS")


plotExp(sc_pred, gene = "CD19")
plotLoadings(sc_pred, pc = 2)



plotEigen(sc_pred, group = "cellType1", geom = "points") %>% 
  plotly::ggplotly()

plotPCVar(sc_pred)


plotTrainProbs(sc_pred, nrow = 1)

trainPreds <- getTrainPred(sc_pred)

trainPreds$T_cell %>% 
  rownames_to_column("id") %>% 
  filter(obs != pred, obs == "other") %>% 
  filter(T_cell > 0.9) -> tCellMiss



pca <- getPCA(sc_pred) %>% 
  as.data.frame() %>% 
  rownames_to_column("id") %>% 
  mutate(class = sc_pred@metadata$cellType1) %>% 
  mutate(tCellMiss = as.factor(id %in% tCellMiss$id))


pca %>% 
  ggplot() +
  aes(x = PC11, fill = tCellMiss) +
  geom_histogram(color = "black", alpha = 0.5) +
  geom_rug() + 
  scale_fill_jcolors() +
  theme_bw()




ggplot(pca, aes(PC10, PC16, color = class)) +
  geom_point(size = 0.1) +
  geom_density_2d() +
  scale_color_jcolors() +
  facet_wrap(~tCellMiss) +
  theme_bw()


plotEigen(sc_pred, group = "cellType1", 10, 16)



sc_pred <- scPredict(sc_pred, pred_data)

predictions <- getPredictions(sc_pred)
predictions$true <- factor(pred_data@meta.data$cellType1)
#saveRDS(predictions, file = here(output, "predictions.RDS"))
predictions <- readRDS(file = here(output, "predictions.RDS"))


# Update probabilities




############

tidy_table <- function(x, left, top, fill = 0, prop = FALSE, digits = 2) {
  x %>%
    group_by_(left, top) %>%
    summarise(n = n()) %>%
    spread(key = top, value = "n", fill = fill) %>%
    as.data.frame() %>%
    column_to_rownames(left) -> x

  if (prop) {
    row_names <- rownames(x)
    x <- mapply(function(x, d) {
      x / d
    }, x, colSums(x))
    rownames(x) <- row_names
    x %>%
      round(digits) %>%
      as.data.frame() -> x
  }

  x
}

############
 

i <- apply(predictions[1:3], 1, which.max)


max.props <- vector(length = nrow(predictions))

for(j in seq_len(nrow(predictions))){
  max.props[[j]] <- predictions[j, 1:3][i[j]]
}

max.props <- unlist(max.props)


predictions$max.prop <- max.props
predictions$max.prop.class <- names(max.props)


# Plot ROC curves
library(plotROC)

plot_roc <- function(x) {
  x$pred$obs <- factor(x$pred$obs, levels = as.character(x$levels))
  x$pred$obs <- ifelse(x$pred$obs == x$levels[1], 1, 0)

  ggplot(
    x$pred,
    aes_string(m = x$levels[1], d = "obs")
  ) +
    geom_roc(hjust = -0.4, vjust = 1.5) +
    coord_equal() +
    style_roc() +
    ggtitle(x$levels[1]) +
    geom_abline(slope = -1, intercept = 1, color = "red")
}



# Get performance results


## Set thresholds for each cell type
thr1 <- 0.98
thr2 <- 0.92
thr3 <- 0.75

# thr1 <- thr2 <- thr3 <- 0.9


## Assign table


predictions <- predictions %>%
  rownames_to_column("id") %>%
  mutate(class = if_else(max.prop > thr1 & max.prop.class == "B_cell", max.prop.class,
    ifelse(max.prop > thr2 & max.prop.class == "Natural_killer", max.prop.class,
      ifelse(max.prop > thr3 & max.prop.class == "T_cell", max.prop.class, "Unassigned")
    )
  ))

predClasses <- split(predictions, predictions$predClass)
x <- predClasses$T_cell
  

ggplot(x, aes(x = T_cell, fill = true)) +
  geom_histogram(color = "gray2") +
  scale_fill_jcolors("pal3") +
  theme_bw() +
  xlab("Probability") +
  ylab("Number of cells")






i <- rownames(pred_data@meta.data)[pred_data@meta.data$cellType1 %in% c("B_cell")]
TCells <- SubsetData(pred_data, cells.use = i, do.clean = TRUE)

predictions %>% 
  filter(true != class) %>% 
  filter(true == "B_cell", class == "T_cell") -> bCellMiss


TCells@meta.data$bCellMiss <- rownames(TCells@meta.data) %in% bCellMiss$id


VlnPlot(TCells, features.plot = c("CCR7"), group.by = "bCellMiss")



tidy_table(predictions, "class", "true", prop = TRUE) %>% 
  select(-Blood_progenitor, -Dendritic_cell)






####################------------------------------------------------------


getProb(sc_pred) %>% 
  lapply("[[", 1) %>% 
  lapply(mean) %>% 
  unlist() -> x


getPred(sc_pred) %>%
  lapply(function(x){
    x %>% 
    filter(obs != "other") %>% 
      pull(other) %>% 
      mean()
  }) %>% 
  lapply(function(x) 1 - x) %>% 
  unlist() -> mean_props


mean_props * mean(mean_props)



a %>% 
  filter(obs != "other") %>% 
  
  select(obs, B_cell) 
#  gather(key = "class", value = "prob", 1:2) %>% 
  ggplot(aes(x = other)) +
  geom_density()


a %>% 
  #  gather(key = "class", value = "prob", 1:2) %>% 
  ggplot(aes(x = other)) +
  geom_density()


x[1] / sum(x)
x[2] / sum(x)
x[3] / sum(x)




props_positive <- getProb(sc_pred) %>% 
  lapply("[", 1) %>% 
  reduce(cbind)

props_negative <- getProb(sc_pred) %>% 
  lapply("[", 2) %>% 
  reduce(cbind)

props_positive <- (props_positive)
props_negative <- props_negative / rowSums(props_negative)


train_models <- sc_pred@train

preds <- getPred(sc_pred)

preds$B_cell[,c(5,6)] <- cbind(props_positive$B_cell, props_negative$other)
preds$Natural_killer[,c(5,6)] <- cbind(props_positive$Natural_killer, props_negative$other)
preds$T_cell[,c(5,6)] <- cbind(props_positive$T_cell, props_negative$other)

train_models$B_cell$pred <- preds$B_cell
train_models$T_cell$pred <- preds$T_cell
train_models$Natural_killer$pred <- preds$Natural_killer



lapply(train_models, plot_roc) %>% 
  plot_grid(plotlist = ., nrow = 1)




####################------------------------------------------------------




## Get metrics









# scPred version: d7661592

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
