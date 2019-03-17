require(Seurat)
require(pryr)
require(Matrix)
require(gridGraphics)


setObject <- function(ct) {
  CreateSeuratObject(raw.data = ct, 
                     min.cells = round(ncol(ct)*0.01),
                     min.genes = 200)
}

addQC <- function(data){
  mito.genes <- grep(pattern = "_MT-", x = rownames(x = data@data), value = TRUE)
  percent.mito <- Matrix::colSums(data@data[mito.genes, ])/Matrix::colSums(data@data)
  
  ribo.genes <- grep(pattern = "_Rps|_Rpl", x = rownames(x = data@data), value = TRUE, ignore.case = TRUE)
  percent.ribo <- Matrix::colSums(data@data[ribo.genes, ])/Matrix::colSums(data@data)
  
  data <- AddMetaData(object = data, metadata = percent.mito, col.name = "percent.mito")
  data <- AddMetaData(object = data, metadata = percent.ribo, col.name = "percent.ribo")
  data
}

plotQC <- function(object, cellType){
  
  p1 <- object@meta.data %>% 
    ggplot(aes_string("nUMI", "percent.mito")) +
    geom_point() +
    ggtitle(cellType) +
    theme_bw()
  
  p2 <- object@meta.data %>% 
    ggplot(aes_string("nUMI", "percent.ribo")) +
    geom_point() +
    ggtitle(cellType) +
    theme_bw() 
  
  p3 <- object@meta.data %>% 
    ggplot(aes_string("nUMI", "nGene")) +
    geom_point() +
    ggtitle(cellType) +
    theme_bw()
  
  cowplot::plot_grid(p1, p2, p3, nrow = 1)
  
}


filterQC <- function(object, attr =  "nGene", low, high){
  FilterCells(object = object,
              subset.names = attr, 
              low.thresholds = low, 
              high.thresholds = high)
}

processSeurat <- function(object){
  object %>% 
  NormalizeData(normalization.method = "LogNormalize", 
                scale.factor = 10000) %>% 
    FindVariableGenes(mean.function = ExpMean, 
                      dispersion.function = LogVMR, 
                      x.low.cutoff = 0.0125, 
                      x.high.cutoff = 3, 
                      y.cutoff = 0.5, 
                      do.plot = FALSE) %>% 
    ScaleData()
}


plotPCA <- function(object, dim1 = 1, dim2 = 2, group = NULL, size = 0.5, digits = 2){
  
  if(!"pca" %in% names(object@dr)){
    stop("No PCA computed!")
  }
  
  pca <- as.data.frame(object@dr$pca@cell.embeddings)
  expVar1 <- round(object@dr$pca@sdev[dim1] / sum(object@dr$pca@sdev) * 100, digits)
  expVar2 <- round(object@dr$pca@sdev[dim2] / sum(object@dr$pca@sdev) * 100, digits)
  
  if(!is.null(group)){
    
    if(!group %in% names(object@meta.data)){
      stop("Grouping variable absent in metadata")
    }
    
    pca[[group]] <- object@meta.data[[group]]
    
  }
  
  
  p <- ggplot(pca) +
    aes_string(x = paste0("PC", dim1), y = paste0("PC", dim2)) 
  
  if(!is.null(group)){  
    p <- p + geom_point(aes_string(color = group), size = size)
    pal <- getPalette(length(unique(pca[[group]])))
    p <- p + scale_color_manual(values = pal)
    
  }else{
    p <- p + geom_point(size = size)
  }
  
  p + 
    xlab(paste0("PC", dim1, " (exp.var ", expVar1, "%)")) + 
    ylab(paste0("PC", dim2, " (exp.var ", expVar2, "%)")) +
    theme_bw() + guides(colour = guide_legend(override.aes = list(size = 2)))
  
  
  
}  

getPalette <- function(n){
  
  if(n < 6){
    c("#29BF12", "#00A5CF", "#DE1A1A", "#574AE2", "#FFBF00")
  }else if(n < 9){
    c("#558aa6", "#B1740F", "#D5006A", "#08585A", "#FFFD98", "#9449d2", "#BBBE64", "#D7263D")
  }else if(n < 13){
    c("#943CB4", "#194D44", "yellow", "#5B6DC8", "#3CA437", "#6B244C", "#6ACDC5", "#DE1A1A", "#BBB53E", "#2A297A", "#995533", "#D590DA")
  }else{
    stop("Too many classes")
  }
  
  
}



