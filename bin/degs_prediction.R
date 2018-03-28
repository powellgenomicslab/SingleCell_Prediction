# expData <- expTrainCounts
# expMetadata <- expTrainMeta
# pVar <- "status"
# positiveClass <- "tumor"
# pValFilter = 0.05
# log2FCFilter = 2
# pseudocount = TRUE
# fitType = 'local'
# method = "per-condition"

getDE <- function(expData, expMetadata, pVar, positiveClass, pValFilter = 0.05, log2FCFilter = 2, pseudocount = TRUE,
                  fitType = 'local', method = "per-condition"){
  

  if(!all(rownames(expData) == rownames(expMetadata))){
    stop("Cell ids in expression matrix and metadata do not match")
  }
  
  
    if(!pVar %in% names(expMetadata)){
    stop("phenoVar must be included in 'expData' variable metadata")
  }  
  
  
  classes <- levels(expMetadata[[pVar]])
  
  if(!positiveClass %in% classes){
    stop("'positiveClass' is not included in phenoVar metadata variable")
  }
  
  negativeClass <- classes[!classes %in% positiveClass]
  
  
  
  
  # Add pseudocount of 1, otherwise sizeFactors of zero will make the estimateDispersions() fail
  # deseq v1
  if(pseudocount){
    exp <- expData + 1
  }else{
    exp <- expData
  }
  
  
  exp %>% 
    t() %>% 
    newCountDataSet(conditions = expTrainMeta[[phenoVar]]) %>% 
    estimateSizeFactors() %>% 
    estimateDispersions(fitType = fitType, method = method) %>% 
    nbinomTest(condA = positiveClass, condB = negativeClass) %>% 
    mutate(adjFoldChange = (baseMeanB - 1) /(baseMeanA - 1)) %>% 
    mutate(adjlog2FoldChange = log2(adjFoldChange))-> deg
  
  # Filter genes by p-value and log2 fold change
  deg %>% 
    filter(padj < pValFilter, abs(adjlog2FoldChange) >  log2FCFilter) %>% 
    arrange(-padj, abs(adjlog2FoldChange)) ->  sigDEG
  
  sigDEG
}

trainDEGModel <- function(expData,
                       expMetadata,
                       features,
                       pVar,
                       positiveClass = NULL,
                       method = "svmPoly",
                       resampleMethod = "cv",
                       seed = NULL,
                       number = 10,
                       returnData = FALSE,
                       savePredictions = FALSE){

  
  if(nrow(expMetadata) == 0){
    stop("No metadata is available")
  }
  
  if(nrow(features) == 0){
    stop("No features are available")
  }
  
  features <- expData[,colnames(expData) %in% features$id]

  # Get response variable
  response <- expMetadata[[pVar]]
  
  if(!is.null(positiveClass)){
    # Get and refactor response variable according to positive class
    # According to twoClassSummary() documentation
    ## "If assumes that the first level of the factor variables corresponds to a relevant result 
    ## but the lev argument can be used to change this."
    
    if(!any(levels(response) == positiveClass)){
      stop(paste0("positiveClass '", positiveClass, "' is not included in prediction variable '", object@pVar, "'"))
    }
    orderedLevels <- c(levels(response)[levels(response) == positiveClass], levels(response)[levels(response) != positiveClass])
    response <- factor(response, levels = orderedLevels)
  }
  
  # Train model
  
  if(!is.null(positiveClass)){
    if(!is.null(seed)) set.seed(seed)
    trCtrl <- trainControl(classProbs = TRUE,
                           method = resampleMethod,
                           number = number,
                           summaryFunction = twoClassSummary,
                           returnData = returnData,
                           savePredictions = savePredictions,
                           allowParallel = FALSE)
    
    fit <- train(x = as.matrix(features), 
                 y = response, 
                 method = method,
                 metric = "ROC",
                 trControl = trCtrl)
  }else{
    if(!is.null(seed)) set.seed(seed)
    trCtrl <- trainControl(method = resampleMethod,
                           number = number,
                           savePredictions = savePredictions,
                           allowParallel = FALSE)
    
    fit <- train(x = as.matrix(features), 
                 y = response, 
                 method = method,
                 trControl = trCtrl,
                 returnData = returnData)
    
  }
  
  
  fit
  
  
  
}
degPredict <- function(features, predData, trainedModel, classes = NULL){
  
  
  if(!is(trainedModel, "train")){
    stop("'trainedModel' object must be of class 'train'")
  }
  
  
  if(length(features) == 0){
    stop("No informative features are present")
  }
  

  predData <- predData[,colnames(predData) %in% features$id]
  
  prediction <- predict(trainedModel, newdata = as.matrix(predData), type = "prob")
  
  if(!is.null(classes) & is.numeric(classes)){
    prediction <- factor(ifelse(prediction[,1] > classes, names(prediction)[1], names(prediction)[2]),
                         levels = names(prediction))
    prediction <- data.frame(class = prediction)
  }
  
  row.names(prediction) <- row.names(predData)
  
  prediction
}
