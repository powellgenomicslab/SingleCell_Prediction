# Script information ------------------------------------------------------

# title: Classification of precursor, myeloid, and lymphoid cells
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


output_dir_name <- "10x_immune_myeloid_vs_lymphoid_prediction" # <------ Output directory

date <- format(Sys.Date(), format = "%Y-%m-%d_")
date <- "2018-10-12_" # <------ Fix date

output <- file.path("results", paste0(date, output_dir_name))

if(!dir.exists(output)){
  dir.create(here(output))
}

# Read data ---------------------------------------------------------------

# Input

input <- file.path("results", "2018-10-12_10x_immune_tree") # <------ Input directory

# Read files

filename <- "train_data.RDS" # <------ Input file
train_data <- readRDS(here(input, filename))

filename <- "pred_data.RDS" # <------ Input file
pred_data <- readRDS(here(input, filename))


# Get feature space -------------------------------------------------------

sc_pred <- getFeatureSpace(train_data, pVar = "cellType2")
rm(train_data)


# Train models ------------------------------------------------------------

library(doParallel)
cl <- makePSOCKcluster(2)
registerDoParallel(cl)
sc_pred <- trainModel(sc_pred)
stopCluster(cl)


# Make predictions --------------------------------------------------------
sc_pred <- scPredict(sc_pred, pred_data)
predictions <- getPredictions(sc_pred)


# Get performance ---------------------------------------------------------

sc_pred@predMeta <-  pred_data@meta.data

sc_pred %>% 
  getPredictions() %>% 
  mutate(true = pred_data@meta.data$cellType2) %>% 
  mutate(predClass = ifelse(predClass == "unassigned", "normal", predClass)) %>% 
  mutate(prediction = ifelse(predClass == true, "correct", "incorrect")) %>% 
  group_by(true, prediction) %>% 
  summarise(n = n()) %>% 
  mutate(accuracy = (n / sum(n))*100) -> performance

performance %>%   
  filter(prediction == "incorrect") %>% 
  select(-prediction) %>% 
  write.table(file = here(output, "accuracy.txt"), quote = FALSE, row.names = FALSE)

performance %>%   
  filter(prediction == "correct") %>% 
  select(-prediction) %>% 
  set_colnames(c("true", "n", "misclassification")) %>% 
  write.table(file = here(output, "misclassification.txt"), quote = FALSE, row.names = FALSE)

saveRDS(predictions, file = here(output, "predictions.RDS"))
saveRDS(sc_pred, file = here(output, "scpred.RDS"))

predictions$true <- pred_data@meta.data$cellType2

predictions %>% 
  gather(key = "class", value = "prob", 1:3) %>% 
ggplot(aes(x = prob, fill = class)) +
  geom_histogram(color = "black") +
  facet_wrap(~true, scales = "free") +
  scale_fill_manual(values = set_names(jcolors::jcolors(), NULL)) +
  theme_bw() +
  xlab("Probability") +
  ylab("Number of cells") -> p

ggsave(p, filename = here(output, "prob_dist.png"), width = 9 , height = 2)


# scPred version: d7661592

# Session info ------------------------------------------------------------

options(width = 70)
capture.output(devtools::session_info(), file = here(output, "session_info.txt"))
