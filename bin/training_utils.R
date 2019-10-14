getTrainPred <- function(object){
  
  if(!is(object, "scPred")){
    stop("'object' must be of class 'scPred'")
  }
  
  if(!length(object@train)){
    stop("No training model found!")
  }
  
  
  preds <- lapply(object@train, "[[", "pred")
  ids <- lapply(object@train, "[[", "trainingData") %>% 
    lapply(rownames)
  
  mapply(function(x, y){ rownames(x) <- y; x}, preds, ids, SIMPLIFY = FALSE)
  
}

plotTrainProbs <- function(object, ...){
  
  preds <- getTrainPred(object)
  namesClasses <- names(preds)
  
  plotProbs <- function(data, class){
    
    ggplot(data, aes_string(class, fill = "obs")) +
      geom_histogram(color = "black", size = 0.1) +
      xlab(paste0("P(", class,")")) +
      ylab("Number of cells") +
      labs(fill = "Prediction") +
      scale_fill_manual(values = c("#1874CD", "#EE2C2C")) +
      theme_bw()
    
  }
  
  plots <- mapply(plotProbs, preds, namesClasses, SIMPLIFY = FALSE) -> plot_res
  n <- length(plots)
  
  if(n > 1){
    legend <- get_legend(plots[[1]] +
                           scale_fill_manual(labels = c("+ class", "- class"),
                                              values = c("#1874CD", "#EE2C2C"))
    )

    plots <- lapply(plots, function(p) p + theme(legend.position = "none"))
    plots[2:n] <- lapply(plots[2:n], function(p) p + ylab(""))
    plots[2:n] <- lapply(plots[2:n], function(p) p + ylab(""))
    plots <- plot_grid(plotlist = plots, ...)
    plots <- plot_grid(plots, legend, rel_widths = c(n, 1/n), rel_heights = c(n, 1/n), ...)


  }
  
  plots
  
}



object@features %>% 
  lapply("[[", 1) %>% 
  unlist() %>% 
  as.character() %>% 
  unique()


plotPCVar <- function(object){
  
  if(!is(object, "scPred")){
    stop("'object' must be of class 'scPred'")
  }
  
  expVar <- object@expVar %>% 
    as.data.frame() %>% 
    rownames_to_column()
  
  names(expVar) <- c("PC", "expVar")
  expVar$PC <- factor(expVar$PC, levels = as.character(expVar$PC), labels = seq_len(nrow(expVar)))
  
  
  p <- expVar %>% 
    ggplot() +
    aes(x = PC, y = expVar, group = 1) +
    geom_point() +
    geom_line(size = 0.1) +
    xlab("Principal component") +
    ylab("Explained variance") +
    theme_bw()
  
  p
}








