library(shiny)
library(shinythemes)
library(tidyverse)
library(ggjoy)
library(knitr)
library(gridExtra)
library(grDevices)
library(RColorBrewer)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinytheme("sandstone"), {
  
  # Application title
  titlePanel("Human Blood Cell Atlas")
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      uiOutput(outputId = 'x_selector'),
      uiOutput(outputId = 'y_selector'),
      uiOutput(outputId = 'fill'),
      uiOutput(outputId = 'facet'),
      tags$hr(),
      uiOutput(outputId = 'eigenvector'),
      uiOutput(outputId = 'de_cutoff')
    ),
    mainPanel(
      HTML('<center>'),
      
      plotOutput(outputId = "plot", width = 800, height = 700),
      
      HTML('</center>')
    )
    
  )
})

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  performance <- reactive({
    
    dc.cells <- paste0("DC", 1:6) 
    mono.cells <- paste0("Mono", 1:4)
    all.cells <- c(dc.cells, mono.cells)
    
    
    # For all clusters...
    # Read DE performance 
    de.perf <- lapply(all.cells, function(cluster){
      readRDS(paste0(getwd(), "/../../results/2017-08-24_hbca_DE_eigenvectors/performance_cluster", cluster, ".RDS"))
    })
    names(de.perf) <- all.cells
    
    # Read PCA performance
    pca.perf <- lapply(all.cells, function(cluster){
      readRDS(paste0(getwd(), "/../../results/2017-08-24_hbca_pca_eigenvectors/performance_cluster", cluster, ".RDS"))
    })
    names(pca.perf) <- all.cells
    
    ## Rename eigenvector column
    
    pca.perf <- pca.perf %>%
      lapply(function(cluster){names(cluster) <- c(names(cluster)[-ncol(cluster)], "n.ev"); cluster})
    
    
    # Read MDS performance
    mds.perf <- lapply(all.cells, function(cluster){
      readRDS(paste0(getwd(), "/../../results/2017-08-24_hbca_mds_eigenvectors/performance_cluster", cluster, ".RDS"))
    })
    names(mds.perf) <- all.cells
    
    # Read DM performance
    dm.perf <- lapply(all.cells, function(cluster){
      readRDS(paste0(getwd(), "/../../results/2017-08-24_hbca_dm_eigenvectors/performance_cluster", cluster, ".RDS"))
    })
    names(dm.perf) <- all.cells
    
    # Aggregate dimension reduction results
    dim.red.perf <- list(pca = pca.perf, mds = mds.perf, dm = dm.perf)
    
    mapply(function(method, name.method){
      
      mapply(function(cluster, name.cluster){
        cluster$cluster <- name.cluster
        cluster
      }, method, names(method), SIMPLIFY = FALSE) %>% reduce(rbind) -> agg.clusters
      
      agg.clusters$method <- name.method
      agg.clusters
      
    }, dim.red.perf, names(dim.red.perf), SIMPLIFY = FALSE) %>% reduce(rbind) -> dim.red.perf
    
    names(dim.red.perf) <- c("ML_model", "Accuracy", "Sensitivity", 
                             "Specificity", "Kappa",  "Eigenvectors", "Cell_type", "Reduction_method")
      dim.red.perf$n.ev.num <- as.numeric(dim.red.perf$Eigenvectors) + 1
      dim.red.perf
  })
  
  
  DE.performance <- reactive({
    dc.cells <- paste0("DC", 1:6) 
    mono.cells <- paste0("Mono", 1:4)
    all.cells <- c(dc.cells, mono.cells)
    
    de.perf <- lapply(all.cells, function(cluster){
      readRDS(paste0(getwd(), "/../../results/2017-08-24_hbca_DE_eigenvectors/performance_cluster", cluster, ".RDS"))
      #readRDS(paste0("results/2017-08-24_hbca_DE_eigenvectors/performance_cluster", cluster, ".RDS"))
    })
    names(de.perf) <- all.cells   
      mapply(function(cluster, name.cluster){
        cluster$Cell_type <- name.cluster
        cluster
      }, de.perf, names(de.perf), SIMPLIFY = FALSE) %>% 
      reduce(rbind) -> de.performance
      names(de.performance) <- c("ML_model", "Accuracy", "Sensitivity", 
      "Specificity", "Kappa", "Cell_type")
      de.performance
  })
  
  
  ## First variable
  observe({
    if(!is.null(performance())){
      output$x_selector <- renderUI({
        choices <-list("Reduction method" = "Reduction_method",
                       "ML model" = "ML_model",
                       "Eigenvectors" = "Eigenvectors",
                       "Cell type" = "Cell_type"
                       )
        selectInput(inputId = 'x', 
                    label = 'Choose x variable',
                    choices = choices)
      })
    }
    
  })
  
  observe({
    if(!is.null(performance())){
      output$eigenvector <- renderUI({
        sliderInput(inputId = 'eigenvector', 
                    label = 'Pick number of eigenvectors', 
                    min = min(as.numeric(performance()[,"Eigenvectors"])) + 1,
                    max = max(as.numeric(performance()[,"Eigenvectors"]) + 1),
                    value = max(as.numeric(performance()[,"Eigenvectors"]) + 1), 
                    animate = TRUE)
      })
    }
    
  })
  
  
  
  ## Second variable
  observe({
    if(!is.null(performance())){
      output$y_selector <- renderUI({
        choices <- c("Accuracy", "Sensitivity", "Specificity", "Kappa")
        selectInput(inputId = 'y', 
                    label = 'Choose performance metric',
                    choices = choices)
      })
    }
    
  })
  
  
  ## Color variable
  observe({
    if(!is.null(performance())){
      output$fill <- renderUI({
        choices <- list('None' = 'None', 
                        "ML model" = "ML_model", 
                        "Cell type" = "Cell_type", 
                        "Reduction method" = "Reduction_method")
        selectInput(inputId = 'fill', 
                    label = 'Choose color/fill variable',
                    choices = choices)
      })
    }
    
  })
  
  ## Facet variable
  observe({
    if(!is.null(performance())){
      output$facet <- renderUI({
        choices <- list('None' = 'None', 
                        "Cell type" = "Cell_type",
                        "ML model" = "ML_model", 
                        "Reduction method" = "Reduction_method", 
                        "Eigenvectors" = "Eigenvectors")
        selectInput(inputId = 'facet', 
                    label = 'Choose faceting variable',
                    choices = choices)
      })
    }
    
  })
  
  observe({
    if(!is.null(performance())){
      output$de_cutoff <- renderUI({
        choices <- as.character(unique(DE.performance()[,"ML_model"]))
        radioButtons(inputId = 'de_cutoff', 
                    label = 'Draw median cutoff according to ML model',
                    choices = choices)
      })
    }
    
  })
  
  
  
  
  facet <- reactive({
    if(input$facet == "None"){
      NA
    }
    else{ 
      input$facet
    }
  })
  
  fill.color <- reactive({
    if(input$fill == "None"){
      NA
    }
    else{ 
      input$fill
    }
  })
  
  
  
  ## Render plot
  output$plot <- renderPlot({
    if(!is.null(performance()) & !is.null(DE.performance()) & !is.null(input$y) & !is.null(input$y)){
        performance() %>% 
          filter(n.ev.num <= input$eigenvector) -> performance.data


      DE.performance() %>% 
        filter(ML_model == input$de_cutoff) %>% 
        gather(key = "Metric", value = "value", 2:5) %>% 
        filter(Metric == input$y) -> DE
        
        # de.performance %>%
        # filter(ML_model == "glm") %>%
        # gather(key = "Metric", value = "value", 2:5) %>% 
        # filter(Metric == "Accuracy") -> DE
      
      
      if(is.na(fill.color())){
        p <- ggplot(performance.data, aes_string(input$x, input$y)) +
          geom_boxplot()
      }else{
        p <- ggplot(performance.data, aes_string(input$x, input$y)) + 
          geom_boxplot(aes_string(fill = fill.color()))
      }
      
      if(!is.na(facet())){
        p <- p + facet_wrap(as.formula(paste("~", facet())))
      }
      p <- p + theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, colour = "black"),
              axis.text.y = element_text(angle = 90, hjust = 1, size = 12, colour = "black"),
              axis.title.x = element_text(size = 20),
              axis.title.y = element_text(size = 20),
              plot.title = element_text(size = 22, face = "bold"),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 13)) +
        scale_fill_brewer(palette = "Set1") +
        geom_hline(data = DE, 
                  aes(yintercept = median(value)), alpha = 0.5)
      p
    }
  })
  
}
# Run the application 
shinyApp(ui = ui, server = server)

