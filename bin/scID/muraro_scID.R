library(tidyverse)
library(dsLib)
library(scID)


output <- setOutput("2019-07-04", "scID_muraro")

trainData <- readRDS("data/2018-04-15_pancreas_baron/baron-human.rds")
testData <- readRDS("data/2018-04-15_pancreas_muraro/muraro.rds")



rownames(testData) <- rownames(testData) %>% str_remove("__.*$")

sourceCellTypes <-  colData(trainData)[,"cell_type1"]
sourceCellTypes <- as.factor(ifelse(sourceCellTypes %in% c("alpha", "beta", "delta", "gamma"), 
                                    as.character(sourceCellTypes), "other"))


commonGenes <-  intersect(rownames(trainData), rownames(testData))


trainData <- trainData [commonGenes,]
testData <- testData[commonGenes,]


trainData@colData$cell_type1 <- sourceCellTypes
meta.data <- as.data.frame(trainData@colData)
meta.data$id <- rownames(meta.data)

true <- testData@colData$cell_type1 %>% as.character()


# Exctract islets of langerhans

i <- sourceCellTypes %in% c("alpha", "beta", "delta", "gamma")
expTrain <- counts(trainData)[,i]
meta.data <- meta.data[i,]

cl <- meta.data[,"cell_type1"]
names(cl) <- rownames(meta.data)

scID_output <- scid_multiclass(target_gem = normcounts(testData), 
                               reference_gem = expTrain, 
                               reference_clusters = cl, 
                               logFC = 0.6, 
                               only_pos = FALSE,  
                               estimate_weights_from_target = FALSE, 
                               normalize_reference = TRUE)


pred <- scID_output$labels



table(pred, true) %>% 
  as.data.frame() %>% 
  set_names(c("pred", "true", "Freq")) %>% 
  spread(key = "true", value = "Freq") %>% 
  column_to_rownames("pred") -> props

mapply(function(x,y) x/y, props, colSums(props)) %>% 
  `rownames<-`(rownames(props)) %>% 
  round(2) %>% 
  .[, c("alpha", "beta", "delta", "gamma", "acinar", "ductal", "endothelial", "epsilon", "mesenchymal")
] %>% write.table(file = here(output, "muraro_prop.txt"), sep = "\t", row.names = FALSE)  


props %>% 
  .[, c("alpha", "beta", "delta", "gamma", "acinar", "ductal", "endothelial", "epsilon", "mesenchymal")] %>% 
  write.table(file = here(output, "muraro_counts.txt"), sep = "\t", row.names = FALSE) 
