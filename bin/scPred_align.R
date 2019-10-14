getInformativePCs <- function(object, classes, pVar, correction = "fdr", sig = 0.1){
  
  
  # Validations -------------------------------------------------------------
  
  if(!any(correction %in% stats::p.adjust.methods)){
    stop("Invalid multiple testing correction method. See ?p.adjust function")
  }
  
  if(!is.factor(classes)){
    stop("Prediction variable must be a factor object")
  }else if(!all(levels(classes) %in% unique(classes))){
    stop("Not all levels are included in prediction variable")
  }else if(length(levels(classes)) == 1){
    stop("No training is possible with only one classification class. Check prediction variable")
  }
  
  
  # Filter principal components by variance ---------------------------------
  
  pca <- object
  
  if(length(levels(classes)) == 2){
    res <- .getPcByClass(levels(classes)[1], object, classes, pca, correction, sig)
    res <- list(res)
    names(res) <- levels(classes)[1]
  }else{
    res <- lapply(levels(classes), .getPcByClass, object, classes, pca, correction, sig)
    names(res) <- levels(classes)
  }
  

  message("DONE!")
  res
  
}

.getPcByClass <- function(positiveClass, object, classes, pca, correction, sig){
  
  i <- classes != positiveClass
  newClasses <- as.character(classes)
  newClasses[i] <- "other"
  newClasses <- factor(newClasses, levels = c(positiveClass, "other"))
  
  
  lapply(pca, function(pc) wilcox.test(pc[newClasses == positiveClass], pc[newClasses == "other"])) %>% 
    lapply('[[', "p.value") %>% 
    as.data.frame() %>% 
    gather(key = "PC", value = "pValue") %>% 
    mutate(pValueAdj = p.adjust(pValue, method = correction, n = nrow(.))) %>% 
    arrange(pValueAdj) %>% 
    filter(pValueAdj < sig) -> sigPCs
  
  sigPCs
}




trainModel <- function(object,
                       method = "svmPoly",
                       resampleMethod = "cv",
                       seed = NULL,
                       number = 10,
                       returnData = FALSE,
                       savePredictions = FALSE){
  
  # Validate class
  
  if(nrow(object$metadata) == 0){
    stop("No metadata has been assigned to object")
  }
  
  if(length(object$features) == 0){
    stop("No features have been determined. Use 'getInformativePCs' function")
  }
  
  
  classes <- object$metadata[[object$pVar]]
  
  if(length(levels(classes)) == 2){
    modelsRes <-  .trainModelByClass(levels(classes)[1],
                                     classes,
                                     object,
                                     method,
                                     resampleMethod,
                                     seed = seed,
                                     number,
                                     returnData,
                                     savePredictions)
    modelsRes <- list(modelsRes)
    names(modelsRes) <- levels(classes)[1]
    
    
  }else{
    modelsRes <- lapply(levels(classes), .trainModelByClass,
                        classes,
                        object,
                        method,
                        resampleMethod,
                        seed,
                        number,
                        returnData,
                        savePredictions)
    names(modelsRes) <- levels(classes)
  }
  
  object$train <- modelsRes
  object
}

.trainModelByClass <- function(positiveClass,
                               classes,
                               object,
                               method,
                               resampleMethod,
                               seed,
                               number,
                               returnData,
                               savePredictions){
  
  if(nrow(object$features[[positiveClass]]) == 0){
    message("No informative principal components were identified for class: ", positiveClass)
  }
  
  
  features <- object$prcomp[, as.character(object$features[[positiveClass]]$PC)]
  
  
  # Get and refactor response variable according to positive class
  # According to twoClassSummary() documentation
  ## "If assumes that the first level of the factor variables corresponds to a relevant result 
  ## but the lev argument can be used to change this."
  
  i <- classes != positiveClass
  response <- as.character(classes)
  response[i] <- "other"
  response <- factor(response, levels = c(positiveClass, "other"))
  
  
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
  fit
}



# Prediction --------------------------------------------------------------


eigenPredict <- function(object, newData, threshold = 0.7, informative = TRUE){
  
  
  if(!(is(newData, "matrix") | is(newData, "data.frame"))){
    stop("'predData' object must be a dataframe or a matrix")
  }
  
  if(length(object$features) == 0){
    stop("No informative principal components have been obtained yet.\nSee getInformativePCs() function")
  }
  
  projection <- newData
  
  
  classes <- names(object$features)
  
  res <- lapply(classes, .predictClass, object, projection)
  names(res) <- levels(classes)
  res <- as.data.frame(res)
  row.names(res) <- rownames(projection)
  
  if(length(classes) == 1){
    classes <- levels(object$metadata[[object$pVar]])
    res[[classes[2]]] <- 1 - res[[classes[1]]]
    
    res$class <- ifelse(res[,1] > threshold, classes[1], 
                        ifelse(res[,2] > threshold, classes[2], "Unassigned")) %>% 
      as.factor()
    
    res[,2] <- NULL 
    names(res) <- c("probability", "class")
    return(res)
  }
  
  i <- apply(res, 1, which.max)
  
  prob <- c()
  for(j in seq_len(nrow(res))){
    prob[j] <- res[j,i[j]]
  }
  
  predictions <- data.frame(probability = prob, prePrediction = names(res)[i])
  rownames(predictions) <- rownames(projection)
  
  predictions %>% 
    rownames_to_column("id") %>% 
    mutate(class = ifelse(probability > threshold, as.character(prePrediction), "Unassigned")) %>%
    select(-prePrediction) %>% 
    column_to_rownames("id") -> finalPrediction
  
  return(finalPrediction)
  
  
}


.predictClass <- function(positiveClass, object, projection){
  # Get features for positive class
  featureList <- object$features[[positiveClass]]
  features <- projection[,as.character(featureList$PC)]
  
  # Perform presictions
  prediction <- predict(object$train[[positiveClass]], 
                        newdata = features, 
                        type = "prob")
  
  prediction[,1, drop = FALSE]
  
}


