# Load libraries ----------------------------------------------------------
library("tidyverse")
library("scPred")
library("here")
library("cowplot")


# Read datasets -----------------------------------------------------------
hbca_cpm <- readRDS(here("results", "2018-04-06_hbca_processing", "hbca_processed_cpm.RDS"))
hbca_counts <- readRDS(here("results", "2018-04-06_hbca_processing", "hbca_processed_counts.RDS"))
hbca_metadata <- readRDS(here("results", "2018-04-06_hbca_processing", "hbca_processed_metadata.RDS"))
dcs_cpm <- readRDS(here("results", "2018-04-08_dcs_processing", "dcs_processed_cpm.RDS"))
output <- file.path("results", "2018-04-08_dcs_prediction")

# Label cells in training dataset -----------------------------------------
hbca_metadata %>% 
  rownames_to_column('cell') %>% 
  mutate(cellType = if_else(grepl("DC", group), "DC", "Mono")) %>% 
  mutate(cellType = as.factor(cellType)) %>% 
  column_to_rownames("cell") -> hbca_metadata

colnames(dcs_cpm) %>% 
  str_split("_") %>% 
  lapply("[", c(1,2)) %>% 
  reduce(rbind) %>% 
  set_rownames(colnames(dcs_cpm)) %>% 
  set_colnames(c("origin", "subtype")) %>% 
  as.data.frame() -> dcs_metadata

# Eigendecompose training gene expression data ----------------------------
hbca <- eigenDecompose(t(hbca_cpm))
metadata(hbca) <- hbca_metadata

# Get informative principal components ------------------------------------
hbca <- getInformativePCs(hbca, pVar = "cellType", sig = 0.05)

# Train model -------------------------------------------------------------
hbca <- trainModel(object = hbca, seed = 66, method = "nb")

# Train model -------------------------------------------------------------
# models <- c("svmRadial", "rf", "earth", "glmnet", "glm", "kknn", "nb")
# pmodels <- lapply(models, function(x) trainModel(object = hbca, seed = 66, method = x))
# names(pmodels) <- models
# 
# options(knitr.table.format = "latex")
# lapply(pmodels, function(x) x@train$DC$results) %>% 
#   lapply(function(x) x[which.max(x$ROC),c("ROC", "Sens", "Spec", "ROCSD", "SensSD", "SpecSD")]) %>% 
#   reduce(rbind) %>% 
#   mutate(model = models) %>%
#   select(!!! rlang::syms(c("model","ROC", "Sens", "Spec", "ROCSD", "SensSD", "SpecSD"))) %>% 
#   arrange(-ROC) %>% 
#   knitr::kable(format = "latex", booktabs = T)

# Classify cells in new dataset -------------------------------------------
dcs_predictions <- eigenPredict(object = hbca, newData = t(dcs_cpm))

all(rownames(dcs_predictions) == rownames(dcs_metadata))

dcs_predictions$origin <- dcs_metadata$origin
dcs_predictions %>% 
  group_by(origin, class) %>% 
  summarise(n = n()) %>% 
  knitr::kable()

dcs_predictions %>% 
  group_by(class) %>% 
  summarise(n = n())


acc <- table(dcs_predictions$class)
acc[1]/ sum(acc)


dcs_predictions %>% 
  mutate(label = if_else(class == "DC", "DC", "Other")) %>% 
  group_by(origin, label) %>% 
  summarise(n = n()) %>% 
  group_by(origin) %>% 
  mutate(classN = sum(n)) %>% 
  mutate(ratio = n/classN) %>% 
  filter(label == "DC")

dcs_predictions %>% 
  mutate(label = if_else(class == "DC", "DC", "Other")) %>% 
  group_by(label) %>% 
  summarise(n = n()) %>% 
  mutate(ratio = n/sum(n))

# 797/ (797 + 24 +  60)
# 
# round(456 / (456 + 5), 2)
# round(341 / (341 + 19 + 60), 2)

round(451 / (451 + 10), 2)
round(346 / (346 + 14 + 60), 2)
# Plot data ---------------------------------------------------------------
dcsProj <- projectNewData(newData = t(dcs_cpm), referenceData = hbca)


dcs_predictions  %>% 
  rownames_to_column("id") %>% 
  mutate(predictions = if_else(class == "Mono", 
                               "Misclassified",
                               if_else(class == "DC",
                                       "DC test",
                                       as.character(class)))) %>% 
  mutate(predictions = as.factor(predictions))%>% 
  column_to_rownames("id") -> dcs_predictions




plotEigen(hbca, group = "cellType", pc = c(1,2), geom = "points") -> p

test <- dcsProj[,c("PC1", "PC2")]
test$cellType <- "DC test"
test$step <- "test"
test$predictions <- dcs_predictions$predictions
p$data$step <- "train"
p$data$predictions <- p$data$cellType
p$data <- rbind(p$data, test)
p$layers[[1]]$aes_params$alpha <- 0.5

p + 
  scale_colour_manual(name = "Cell Type", 
                      labels = c("DCs-Train (Villani)", "Monocytes-Train (Villani)", "DCs-Test (Breton)"), 
                      values = c("#fa8775", "#ffd700", "#cd34b5")) +
  geom_density_2d() +
  xlab(paste0("PC1 ", "(", hbca@expVar[1], "% exp. var.)")) +
  ylab(paste0("PC2 ", "(", hbca@expVar[2], "% exp. var.)")) +
  xlim(c(-55, 80)) +
  ylim(c(-65, 35)) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15 ),  
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)) -> p


# ggsave(filename = "dcs_pca_projection.png", plot = p, width = 6.3, height = 4, dpi = 350, path = here(output))



removeCells <- which(p$data$predictions == "Misclassified" | p$data$predictions == "Unassigned")

densitySubset <- p$data[-removeCells,]
densitySubset$predictions <- factor(densitySubset$predictions, levels = levels(densitySubset$predictions)[1:3])


ggplot(p$data, aes(x = PC1, y = PC2, color = predictions)) +
  geom_point(alpha = 0.5) +
  geom_density_2d(aes(x = PC1, y = PC2, color = predictions), data = densitySubset) +
  scale_colour_manual(name = "Cell Type", 
                      labels = c("DCs-Train (Villani)", "Monocytes-Train (Villani)", "DCs-Test (Breton)", "Misclassified cells", "Unassigned"), 
                      values = c("#fa8775", "#ffd700", "#cd34b5", "blue", "black")) +
  xlab(paste0("PC1 ", "(", hbca@expVar[1], "% exp. var.)")) +
  ylab(paste0("PC2 ", "(", hbca@expVar[2], "% exp. var.)")) +
  xlim(c(-55, 80)) +
  ylim(c(-65, 35)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15 ),  
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11))-> p2



incorrect <- p$data[removeCells,]
incorrect$origin <- as.factor(ifelse(grepl("^cb_.*", rownames(incorrect), perl = TRUE), "cb", "blood"))

incorrect %>% 
  group_by(predictions, origin) %>% 
  summarise(n = n())


pm <- plot_grid(p, p2, hjust = TRUE, vjust = TRUE, nrow = 1)
ggsave("dcs_prediction.png", plot = pm, width = 13, height = 3.8, dpi = 350, path = here(output))

# scPred version: 60836d183836fe1132694952e230402e6963907f



hbca_metadata %>% 
  rownames_to_column('cell') %>% 
  filter(cellType == "DC") %>% 
  mutate(cellType = factor(group, levels = sort(unique(group)))) %>% 
  column_to_rownames("cell") -> hbca_dc_metadata



# Filter data -------------------------------------------------------------
hbca_dc_cpm <- hbca_cpm[ ,colnames(hbca_cpm) %in% rownames(hbca_dc_metadata)]
dcs_predictions_subset <- dcs_predictions[dcs_predictions$class == "DC",]

dcs_subset_cpm <- dcs_cpm[,colnames(dcs_cpm) %in% rownames(dcs_predictions_subset)]

i <- rownames(dcs_metadata) %in% rownames(dcs_subset_cpm)
dcs_subset_metadata <- dcs_metadata[i,]

all(rownames(dcs_subset_metadata) == rownames(dcs_subset_cpm))

# Eigendecompose training gene expression data ----------------------------
hbca_dc <- eigenDecompose(t(hbca_dc_cpm))
metadata(hbca_dc) <- hbca_dc_metadata

# Get informative principal components ------------------------------------
hbca_dc <- getInformativePCs(hbca_dc, pVar = "cellType", sig = 0.05)

# Train model -------------------------------------------------------------
hbca_dc <- trainModel(object = hbca_dc, seed = 66)

# Classify cells in new dataset -------------------------------------------
predictions_sub_dcs <- eigenPredict(object = hbca_dc, newData = t(dcs_subset_cpm))


