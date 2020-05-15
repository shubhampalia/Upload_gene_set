library(shiny)


ui<- navbarPage(
  title = "Heatmap",
  
  tabPanel(
    title = "Upload Gene Set",
    
    sidebarLayout(
      sidebarPanel(
        fileInput("file1", "Choose CSV/ .txt File",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        # Horizontal line ----
        tags$hr(),
        
        # Input: Select number of rows to display ----
        radioButtons("disp", "Display",
                     choices = c(Head = "head",
                                 All = "all"),
                     selected = "head")
      ),
      
      mainPanel(
        tabsetPanel(type = "tab",
                    
                    tabPanel("User defined heatmap",
                             plotOutput("userhm")),
                    
                    
                    
                    tabPanel("Data-set",
                             DT::dataTableOutput("genesample"))
        )
      )
    )
  )
)
                    
      
server<- function(input, output){
  output$userhm <- renderPlot({
    req(input$file1)
    
    f2 <- read.csv(input$file1$datapath)
    
    mm31= match(f2[,1], ee1$Name)
    mm3 = match( paste0("gene:",f2[,1]), ee1$Name)
    
    for ( i in c(1:length(mm3))){
      if(is.na(mm3[i]==TRUE)){
        mm3[i]=mm31[i]
      }
    }
    
    
    colids3 = grep("RPKM",colnames(ee1))
    userdf = as.matrix( ee1[mm3,colids3])
    coln3 = colnames(ee1)[colids3]
    
    rownames(userdf) = f2[,1]
    colnames(userdf) = sapply( coln3, function(x) strsplit(x,"_sequence")[[1]][[1]])
    
    goodid3 = setdiff( rownames(userdf), rownames(userdf)[is.na(userdf[,1])])
    userdf = userdf[which(rownames(userdf) %in% goodid3),]
    dict3 = unlist(tt[,2])
    names(dict3) = sapply( unlist(tt[,1]), function(x) strsplit(x,"_sequence")[[1]][[1]] )
    colnames(userdf) = dict3[colnames(userdf)]
    pheatmap::pheatmap(log10(userdf+0.1),fontsize_row = 2, cluster_cols = FALSE, cluster_rows = T,
                       annotation_names_col = TRUE, annotation_names_row = TRUE,
                       labels_col = c("endosperm 1 ", "endosperm 2", "endosperm 3",
                                      "pericarp 1", "pericaprp 2", "pericarp 3",
                                      "embryo 1", "embryo 2", "embryo 3",
                                      "aleurone 1", "aleurone 2", "aleurone 3"), 
                       fontsize_col = 14, angle_col = 90, treeheight_row = 0, treeheight_col = 0,
                       main = "HEATMAP AS PER THE GENE FILE UPLOADED BY THE USER")
    
    output$genesample<- DT::renderDataTable({
      req(userdf) 
      
      #df <- read.csv(userdf)
      
      if(input$disp == "head") {
        return(head(userdf))
      }
      else {
        return(userdf)
      }
    })
    
  })
  
  
  
}

shinyApp(ui= ui, server = server)