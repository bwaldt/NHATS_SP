#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(haven)
library(PCAmixdata)
library(klaR)
library(cluster)
library(car)
library(nnet)
library(rmngb)


load(file = "spfile.Rda")
newSp <- as.data.frame(unclass(spfile))
newSp <- zap_formats(newSp)
newSp[sapply(newSp, is.atomic)] <- lapply(newSp[sapply(newSp, is.atomic)],as.factor)


newSp <-  newSp %>% 
  dplyr::select(-spid) %>% # remove id 
  dplyr::select(-r1dresidr, -r1breakoffst, -r1breakoffqt) %>%  # remove logistical type questions
  dplyr::select(-hc1agehrtsrg, -hc1agebcksrg, -hc1agecatsur, -hc1agehipsur, -hc1ageknesur, -hc1hosovrnht, -hc1dementage) %>% # remove time since last surgey
  dplyr::filter(hc1health != -9 & hc1health != -8 & hc1health != -1 ) %>% # people whose overall health we dont now
  dplyr::select(-ht1mthslived, -ht1yrslived, -ht1spacename, -fl1facility) %>% # 
  dplyr::select(-hh1yrendmarr,-hh1proxlivsp) %>%
  dplyr::select(-cs1dreconcil) 


sampleSp <- newSp %>% filter(is1resptype != 2)


sample_subset <- sampleSp[,grepl( "hc1|ss1|pc1|cg1|mo1|dm1|sc1|ds1|mc1meds|pa1|wb1" , names( sampleSp ) )]

sample_subset <- sample_subset %>% dplyr::select(-mc1medsrem, -mc1medsdif, -mc1medsyrgo, -mc1medsslf)


spAnal <- sampleSp

convertedDistance <- readRDS(file = "sample_subset.rds")

colNames <- colnames(spAnal)
checkboxGroupList <- setNames(as.list(colNames[1:length(colNames)]), colNames[1:length(colNames)])
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("NHATS - Sample Survey Participants Only (n=7020)"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("clusters",
                     "Number of Clusters:",
                     min = 1,
                     max = 10,
                     value = 7),
         sidebarPanel(
           checkboxGroupInput(inputId = "displayVariables",
                              label = h4("Select Variables"),
                              choices = checkboxGroupList,
                              selected = 2
           ),
           verbatimTextOutput("varLength"),
           fluidRow("Wait for about 10s for plotting when you select more than 2 attributes."),
           width = 2
         )
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        h4("Size of Each Cluster in Percantage"),
        plotOutput("distPlot"),
        h4("P-Value on Plots shows the probabilty that ALL the clusters are the same with respect to variable"),
        uiOutput("plots"),
        actionButton("Plot", "Draw Tables"),
        h4("Table shows the percentage breakdown of each Cluster by given Variable. Same info as graphs above"),
        dataTableOutput("table"),
        h4("Pairwise Chi-Square Test, Only shows cluster Pairs with p-value over 0.05"),
        h5("Interpretation: P-Value is the Probability that Cluster A and Cluster B are the same with respect to varibale in row."),
        dataTableOutput("tablePValues")
      )
      
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  changeCluster<-reactive({return(pam(convertedDistance, k=input$clusters))})
  
  output$distPlot <- renderPlot({
    # generate bins based on input$bins from ui.R
    
    spAnal$clusterDim <- changeCluster()$cluster
    
    spAnal %>% 
      group_by(clusterDim) %>%
      summarise(
        count = n()
      ) %>%
      mutate(freq = count / sum(count)) %>%
      ggplot(aes(x=as.factor(clusterDim), y=freq)) + geom_bar(position = 'dodge' ,stat='identity') + labs(title = "Cluster Size") + geom_text(aes(label=round(freq,2)), position=position_dodge(width=0.9), vjust=-0.25)
  })
  
  output$value <- renderPrint({ changeCluster()$objective })
  # output$selectedVars <- renderPrint({ colNames[input$displayVariables] })
  output$varLength <- renderPrint({ length(input$displayVariables) })
  
  v <- reactiveValues(doPlot = FALSE)
  # Detect whether the button is clicked
  observeEvent(input$Plot, {
    v$doPlot <- input$Plot
  })
  
  makePlots <- reactive({
    # If the button is not clicked, do not draw
    if (v$doPlot == FALSE) {return()}
    
    # Isolate the plot code, maintain the old plot until the button is clicked again
    isolate({
      var_len <- length(input$displayVariables)
      varNames <- input$displayVariables
      
      if (var_len == 1) {
        spAnal$clusterDim <- changeCluster()$cluster
        varName <-  varNames[1]
        grp_cols <- c("clusterDim",varName)
        
        # Convert character vector to list of symbols
        dots <- lapply(grp_cols, as.symbol)
        
        # Perform frequency counts
        plotData <- spAnal %>%
          group_by_(.dots=dots) %>%
          summarise(
            count = n()
          ) %>%
          mutate(freq = count / sum(count))
        
        p_value <- chisq.test(table(spAnal$clusterDim, spAnal[[varName]]))[3]
        
          plot <- ggplot(plotData) + 
            geom_bar(aes(x=as.factor(clusterDim), y=freq, fill = as.factor(plotData[[varName]]) ),position = 'dodge' ,stat='identity') + 
            labs(title = paste0("Cluster By ",varName, " P-Value: ", p_value))
          return(plot)
          }
      
      if (var_len > 1) {
        
        
        plot_list = list()
        for (i in 1:var_len){
          varName <-  varNames[i]
          spAnal$clusterDim <- clusterPam$cluster
        
          grp_cols <- c("clusterDim",varName)
        
          # Convert character vector to list of symbols
          dots <- lapply(grp_cols, as.symbol)
        
          # Perform frequency counts
          plotData <- spAnal %>%
            group_by_(.dots=dots) %>%
            summarise(
              count = n()
            ) %>%
            mutate(freq = count / sum(count))
          
          p_value <- chisq.test(table(spAnal$clusterDim, spAnal[[varName]]))[3]
          
          p <- ggplot(plotData) + 
            geom_bar(aes(x=as.factor(clusterDim), y=freq, fill = as.factor(plotData[[varName]]) ),position = 'dodge' ,stat='identity') + 
            labs(title = paste0("Cluster By ",varName, " P-Value: ", p_value))

        
        }
        return(multiplot(plot_list))
      }
      
    })
    
  })
  
  output$table <- renderDataTable(makeTable())
  makeTable <- reactive({
    # If the button is not clicked, do not draw
    if (v$doPlot == FALSE) {return()}
    
    # Isolate the plot code, maintain the old plot until the button is clicked again
    isolate({
      var_len <- length(input$displayVariables)
      varNames <- input$displayVariables
      if (var_len == 1) {
        spAnal$clusterDim <- changeCluster()$cluster
        varName <- varNames[1]
        table1 <- table(spAnal[[varName]], spAnal$clusterDim)
        table1 <- prop.table(table1, margin=2)*100
        table <- round(table1,2)
        rownames(table) <- paste0(varName," ",levels(spAnal[[varName]]))
        table <- cbind(names = rownames(table), table)
        return(table)
        }
      
      if (var_len > 1) {
        table <- NULL
        for (i in 1:var_len){
            spAnal$clusterDim <- changeCluster()$cluster
            varName <- varNames[i]
            if (nlevels(spAnal[[varName]]) < 12){
            
              table1 <- table(spAnal[[varName]], spAnal$clusterDim)
              table1 <- prop.table(table1, margin=2)*100
              table1 <- round(table1,2)
              rownames(table1) <- paste0(varName," ",levels(spAnal[[varName]]))
              table1 <- cbind(names = rownames(table1), table1)
              
              table <- rbind(table,table1)
            } else {
              df <-  data.frame(cluster = spAnal$clusterDim, var = spAnal[[varName]])
              means <-  df %>%
                group_by(cluster) %>%
                summarise(
                  mean = mean(as.numeric(var))
                ) 
              
              newTab <- as.table(t(means$mean))
              newTab <- round(newTab,2)
              rownames(newTab) <- varName
              newTab <- cbind(names = rownames(newTab), newTab)
              table <- rbind(table, newTab) 
            }
            
        }
        return(table)
      }

    })
    
  })

  output$tablePValues <- renderDataTable(makeTablePvalue())
  makeTablePvalue <- reactive({
    # If the button is not clicked, do not draw
    if (v$doPlot == FALSE) {return()}
    
    # Isolate the table code, maintain the old table until the button is clicked again
    isolate({
      varNames <- input$displayVariables
      returnedTable <- NULL
      for (i in 1:length(varNames)){
        spAnal$clusterDim <- changeCluster()$cluster
        varName <- varNames[i] 
        x <- spAnal$clusterDim
        y <-  spAnal[[varName]]
        table <- table(spAnal$clusterDim , spAnal[[varName]] )
        i <- (colSums(table) != 0)
        p_values <- pairwise.chisq.test(table[,i])
        
        x <- which(p_values$p.value >.05 ,arr.ind = TRUE)
        
        finalTable <- cbind(x , p_values$p.value[which(p_values$p.value >.05)])
        
        
        colnames(finalTable) <- c("Cluster A","Cluster B","P-Value")
        
        if (nrow(finalTable) > 0) {
            for (j in 1:nrow(finalTable)){
              rownames(finalTable)[j] <- varName
            }
        }
        finalTable <- cbind(names = rownames(finalTable), finalTable)
        
        returnedTable <- rbind(returnedTable, finalTable)
      }
      return(returnedTable)
      
      
      
    })
    
  })
  
    
  
  plotInput <- reactive({
    varNames <- input$displayVariables
    n_plot <- length(varNames)
    return (list("n_plot"=n_plot, "varNames"=varNames))
  })
  
  # Insert the right number of plot output objects into the web page
  output$plots <- renderUI({
    
    plot_output_list <- lapply(1:plotInput()$n_plot, function(i) {
      plotname <- paste("plot", i, sep="")
      plotOutput(plotname)
    })
    
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, plot_output_list)
  })
  
  
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  observe({
    varNames <- plotInput()$varNames
    for (i in 1:plotInput()$n_plot) {
      # Need local so that each item gets its own number. Without it, the value
      # of i in the renderPlot() will be the same across all instances, because
      # of when the expression is evaluated.
      local({
        my_i <- i
        plotname <- paste("plot", my_i, sep="")
        
        output[[plotname]] <- renderPlot({
                  varName <-  varNames[my_i]
                  spAnal$clusterDim <- changeCluster()$cluster
                  
                  grp_cols <- c("clusterDim",varName)
                  
                  # Convert character vector to list of symbols
                  dots <- lapply(grp_cols, as.symbol)
                  
                  if (nlevels(spAnal[[varName]]) < 12) {
                  
                  # Perform frequency counts
                  plotData <- spAnal %>%
                    group_by_(.dots=dots) %>%
                    summarise(
                      count = n()
                    ) %>%
                    mutate(freq = count / sum(count))
                  
                  M <- table(spAnal$clusterDim, spAnal[[varName]])
                  M <- M[,which(!apply(M,2,FUN = function(x){all(x == 0)}))]
                  
                  p_value <- chisq.test(M)[3]
                  
                  ggplot(plotData) + 
                    geom_bar(aes(x=as.factor(clusterDim), y=freq, fill = as.factor(plotData[[varName]]) ),position = 'dodge' ,stat='identity') + 
                    labs(title = paste0("Cluster By ",varName, " P-Value: ", p_value))
                  
                  } else{
                    spAnal[[varName]] <- as.numeric(as.character(spAnal[[varName]]))
                    
                    df      <-    data.frame(cluster = spAnal$clusterDim, var = spAnal[[varName]])
                
                    
                    p <- Anova(multinom(as.factor(cluster) ~ as.numeric(var), data = df))
                    p_value <- p$`Pr(>Chisq)`[1]
                    
                    df %>%
                    group_by(cluster) %>%
                      summarise(
                        mean = mean(var)
                      )  %>%
                      ggplot() + 
                      geom_bar(aes(x=as.factor(cluster), y=mean),position = 'dodge' ,stat='identity') + 
                      labs(title = paste0("Cluster By ",varName, " P-Value: ", p_value))
                  
                    
                    }
                  
                  
                  
                    })
              })
          }
      })
}


# Run the application 
shinyApp(ui = ui, server = server)

