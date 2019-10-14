library(tidyverse)
library(dsLib)
library(singleCellNet)

output <- setOutput("2019-07-04", "scn_segerstolpe")


trainData <- readRDS("data/2018-04-15_pancreas_baron/baron-human.rds")
testData <- readRDS("data/2018-04-15_pancreas_segerstolpe/segerstolpe.rds")



sourceCellTypes <-  colData(trainData)[,"cell_type1"]
sourceCellTypes <- as.factor(ifelse(sourceCellTypes %in% c("alpha", "beta", "delta", "gamma"), 
                                    as.character(sourceCellTypes), "other"))


commonGenes <-  intersect(rownames(trainData), rownames(testData))


trainData <- trainData [commonGenes,]
testData <- testData[commonGenes,]


trainData@colData$cell_type1 <- sourceCellTypes
meta.data <- as.data.frame(trainData@colData)
meta.data$id <- rownames(meta.data)

# targetCellTypes <- testData@colData$cell_type1 %>% as.character()
# true <- as.factor(ifelse(targetCellTypes %in% c("alpha", "beta", "delta", "gamma"), 
#                          as.character(targetCellTypes), "other"))
# 


true <- testData@colData$cell_type1 %>% as.character()


# Exctract islets of langerhans

i <- sourceCellTypes %in% c("alpha", "beta", "delta", "gamma")
expTrain <- logcounts(trainData)[,i]
meta.data <- meta.data[i,]


class_info <- scn_train(stTrain = meta.data, expTrain = expTrain, dLevel = "cell_type1", 
                        colName_samp = "id")


classRes_val_all = scn_predict(cnProc = class_info[["cnProc"]], expDat = logcounts(testData))
classRes_val_all <- classRes_val_all[,colnames(testData)]

pred <- apply(classRes_val_all, 2, function(x){
  if(any(x > 0.9)){ 
    which.max(x)
  }else{
    6
  }
}) %>% 
  c(rownames(classRes_val_all), "other")[.]

length(pred)
length(true)

level <- c("alpha", "beta", "delta", "gamma", "acinar", "co-expression", "ductal", "endothelial", "epsilon", "mast", "MHC class II", "PSC")

table(pred, true) %>% 
  as.data.frame() %>% 
  spread(key = "true", value = "Freq") %>% 
  column_to_rownames("pred") -> props

mapply(function(x,y) x/y, props, colSums(props)) %>% 
  `rownames<-`(rownames(props)) %>% 
  round(2) %>% 
  .[,level] %>% 
  write.table(file = here(output, "segerstolpe_prop.txt"), sep = "\t", row.names = FALSE)  


props %>% 
  .[, level] %>% 
  write.table(file = here(output, "segerstolpe_counts.txt"), sep = "\t", row.names = FALSE)









