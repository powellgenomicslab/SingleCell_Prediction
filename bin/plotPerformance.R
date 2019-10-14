library(here)
library(tidyverse)

replicates <- list.dirs(path = here("results/2018-04-05_prediction_hcc/"), full.names = FALSE)

models <- c("svmLinear", "svmPoly", "svmRadial", "glmnet", "rf")
names(models) <- models


# Get scPred results ------------------------------------------------------


## Test

getPredROC <- function(model){
  i <- grep(paste0("^scPred_seed_.+_", model, "$"), replicates, perl = TRUE)
  replicates[i] %>% 
    lapply(function(x) readRDS(here("results","2018-04-05_prediction_hcc", x, "roc.RDS"))) %>%
    lapply(function(x) as.numeric(x$auc)) %>% 
    unlist()
}

lapply(models, getPredROC) %>% 
  as.data.frame() %>% 
  gather(key = "model", value = "AUROC") %>% 
  mutate(step = "test", method = "scPred") -> rocTest


## Train

getTrainROC <- function(model){
  i <- grep(paste0("^scPred_seed_.+_", model, "$"), replicates, perl = TRUE)
  replicates[i] %>% 
    lapply(function(x) readRDS(here(file.path("results","2018-04-05_prediction_hcc", x, "trained_model.RDS")))) %>%
    lapply(function(x) max(x$results$ROC)) %>% 
    unlist()
}



lapply(models, getTrainROC) %>% 
  as.data.frame() %>% 
  gather(key = "model", value = "AUROC") %>% 
  mutate(step = "train", method = "scPred") -> rocTrain



# Get DEGs results --------------------------------------------------------

# Test

getPredROC <- function(model){
  i <- grep(paste0("^degs_seed_.+_", model, "$"), replicates, perl = TRUE)
  replicates[i] %>% 
    lapply(function(x) readRDS(here(file.path("results","2018-04-05_prediction_hcc", x, "roc.RDS")))) %>%
    lapply(function(x) as.numeric(x$auc)) %>% 
    unlist()
}

lapply(models, getPredROC) %>% 
  as.data.frame() %>% 
  gather(key = "model", value = "AUROC") %>% 
  mutate(step = "test", method = "DEGs") -> rocDEGTest


## Train

getTrainROC <- function(model){
  i <- grep(paste0("^degs_seed_.+_", model, "$"), replicates, perl = TRUE)
  replicates[i] %>% 
    lapply(function(x) readRDS(here(file.path("results","2018-04-05_prediction_hcc", x, "trained_model.RDS")))) %>%
    lapply(function(x) max(x$results$ROC)) %>% 
    unlist()
}


lapply(models, getTrainROC) %>% 
  as.data.frame() %>% 
  gather(key = "model", value = "AUROC") %>% 
  mutate(step = "train", method = "DEGs") -> rocDEGTrain




###

rocRes <- rbind(rocTrain, rocTest, rocDEGTrain, rocDEGTest)
rocRes$step <- factor(rocRes$step, levels = c("train", "test"))



rocRes %>% 
  filter(step == "test") %>% 
  group_by(model, method) %>% 
  summarise(roc = median(AUROC))



ggplot(rocRes, aes(x = model, y = AUROC)) +
  geom_boxplot(aes(fill = step)) +
  facet_wrap(~method) +
  scale_fill_brewer(palette = "Set1")



ggplot(rocRes, aes(x = model, y = AUROC)) +
  geom_boxplot(aes(fill = method)) +
  facet_wrap(~step) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw()



rocRes %>% 
  filter(model == "rf", step == "test") %>% 
  ggplot(aes(x = AUROC, fill = method)) +
  geom_density(alpha = 0.5)


rocRes %>% 
  filter(model == "rf", step == "test") %>% 
  split(.$method) -> methods

wilcox.test(methods$DEGs$AUROC, methods$scPred$AUROC, alternative = "greater")



rocRes %>% 
  filter(model == "rf", step == "test") %>% 
  group_by(method) %>% 
  summarise(median = median(AUROC)) %>% 
  pull(median) %>% 
  diff() %>% 
  abs()



############


# Plot ROC curves ---------------------------------------------------------




getPredROC <- function(model){
  i <- grep(paste0("^degs_log2_seed_.+_", model, "$"), replicates, perl = TRUE)
  replicates[i] %>% 
    lapply(function(x){ 
      x <- readRDS(here(file.path("results","2018-04-05_prediction_hcc", x, "roc.RDS")))
      x[c("specificities", "sensitivities")] %>% as.data.frame() %>% mutate(model = model)
    })
  }

lapply(models, getPredROC) %>% 
  lapply(function(x) Reduce(rbind, x)) %>% 
  Reduce(rbind, .) %>% 
  mutate(method = "DEGs") -> rocDEG

getPredROC <- function(model){
  i <- grep(paste0("^scPred_seed_.+_", model, "$"), replicates, perl = TRUE)
  replicates[i] %>% 
    lapply(function(x){ 
      x <- readRDS(here(file.path("results","2018-04-05_prediction_hcc", x, "roc.RDS")))
      x[c("specificities", "sensitivities")] %>% as.data.frame() %>% mutate(model = model)
    })
}

lapply(models, getPredROC) %>% 
  lapply(function(x) Reduce(rbind, x)) %>% 
  Reduce(rbind, .) %>% 
  mutate(method = "scPred") -> rocscPred

resROC <- rbind(rocDEG, rocscPred)
resROC$Method<- factor(resROC$method, levels = c("scPred", "DEGs"))


p <- ggplot(resROC, aes(x = 1 - specificities, y = sensitivities, col = Method)) +
  stat_summary(fun.y = "median", geom = "line", size = 0.7) +
  facet_wrap(~model) +
  xlab("Median false polsitive rate") +
  ylab("Median true positive rate") +
  theme_bw() +
  scale_color_manual(values = c("#D84438", "#7953CC")) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11 ),  
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        strip.background = element_rect(fill = "#EEEEEE"),
        legend.position = c(0.85, 0.25))


ggsave("hcc_results.png", plot = p, width = 9, height = 5, dpi = 350, 
       path = here("results", "2018-04-05_prediction_hcc"))


resROC %>% 
  filter(model == "svmRadial") %>% 
ggplot(aes(x = 1 - specificities, y = sensitivities, col = Method)) +
  stat_summary(fun.y = "median", geom = "line", size = 0.7) +
  xlab("Median false polsitive rate") +
  ylab("Median true positive rate") +
  theme_bw() +
  scale_color_manual(values = c("#D84438", "#7953CC")) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11 ),  
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        strip.text.x = element_text(size = 11),
        strip.background = element_rect(fill = "#EEEEEE")) -> p2

p2  + annotate("text", x = 0.75, y = 0.75, label = "AUC = 0.97", col = "#D84438") +
  annotate("text", x = 0.75, y = 0.68, label = "AUC = 0.77", col = "#7953CC") -> p2



ggsave("hcc_results_svmRadial.png", plot = p2, width = 5.5, height = 4, dpi = 350, 
       path = here("results", "2018-04-05_prediction_hcc"))


rocTest %>% group_by(model) %>% summarise(m = median(AUROC))
rocDEGTest %>% group_by(model) %>% summarise(m = median(AUROC))
# scPred version: 60836d183836fe1132694952e230402e6963907f