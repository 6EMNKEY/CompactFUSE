#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyr)
library(dplyr)
library(stringr)
#INSTAALLL -v-
library(shinydashboard)
library(ggplot2)
#INSTALL -v-
library(DT)
library(RColorBrewer)
library(ggtext)

# Define UI for application that draws a histogram
ui <- dashboardPage(
    # Application title
    dashboardHeader(title = "CompactFuse"),

    # Sidebar with a slider input for number of bins 

        #h3/4/5()for headers
      dashboardSidebar(
        sidebarMenu(
          menuItem("Known", tabName = "Known"),
          menuItem("Novel", tabName = "Novel"),
          menuItem("Individual", tabName = "Individual"),
          menuItem("Transposones(dFam)", tabName= "TransposonesTab"),
          menuItem("Breakpoint Analisis", tabName= "breakpoint")
        )
        #actionButton(label="toggle_input", "Help")
      ),
      dashboardBody(
        
        tabItems(
          
          tabItem(tabName = "Known", 
                wellPanel(style = "background: #e9edef",
                  fluidRow(
                  column(6,selectInput("DiseaseT", "Select target disease:",append("All",MorphDict$V2),selected = "All", multiple = TRUE)),
                  column(6,selectInput("TopologyT", "Select target tissue:",append('All',TopDict$V2),selected = "All",multiple = TRUE)),
                  column(6,sliderInput("SpanningR", "Number of spanning reads:", min =0,max = max(dataM$spanning.reads),value =c(0,max(dataM$spanning.reads)))),
                  column(6,sliderInput("SpanningP", "Number of spanning pairs:",min =0,max = max(dataM$spanning.pairs),value =c(0,max(dataM$spanning.pairs)))),
                  column(6,sliderInput("Rank", "Rank from 0 to 4", min = 0 , max = 4, value = 0)),
                  column(6,selectInput("Oncogene", "Oncogeneicity: ", c("ALL", "ONE", "BOTH"))),
                  column(6,submitButton("Search"))
                  )
              ),
          DT::dataTableOutput(("filtered")),
          downloadButton("downloadF", "Download Data")
                  
                  
 
                  ),
          tabItem(tabName = "Novel", 
                    #INPUTSES --> Sliders, Buttons etc and position =left/right
                wellPanel(style = "background: #e9edef",
                  fluidRow(
                  column(6,sliderInput("SpanningRN", "Number of spanning reads:", min =0,max = max(dataM$spanning.reads),value =c(0,max(dataM$spanning.reads)))),
                  column(6,sliderInput("SpanningPN", "Number of spanning pairs:",min =0,max = max(dataM$spanning.pairs),value =c(0,max(dataM$spanning.pairs)))),
                      # column(6,sliderInput("Rank", "Rank from 0 to 4", min = 0 , max = 4, value = 0)),
                  column(6,selectInput("OncogeneN", "Oncogeneicity: ", c("ALL", "ONE", "BOTH"))),
                  column(6,submitButton("Search"))
                  )
                  )
                      
                    # Use conditionalPanel to hide/show input page
,
          DT::dataTableOutput(("filteredNovels")),
          downloadButton("downloadFN", "Download Data"),
                  ),
          tabItem(tabName = "Individual", 
                wellPanel(style = "background: #e9edef",
                  fluidRow(
                  column(6,selectInput("FusionIn", "Select an individal Fusion:", c(TOTAL$Fusion),selected = TOTAL$Fusion[1] )),
                  column(6,selectInput("Modification","Modifier(select a column):", c(append("none",colnames(vcf_txt))), selected="none")),
                  column(6,selectInput("TransChrom","Ending chromosome for translocations:", c(append("none",colnames(vcf_txt))), selected="none")),
                  column(6,radioButtons("ModType", "Type of modification:", choiceNames = c("End Position","Confidence interval (n,n) format", "Length", "Sequence in nucleotides (length function applied)"), choiceValues = c(1,2,3,4))),
                  column(6,sliderInput("depth", "Depth of the search in bp:", min = 2000, max = 40000000, value = 20000, step = 1)),
                  column(6,submitButton("Search"))
                  )
          ),
          DT::dataTableOutput(("selected")),
                  tags$head(
                    tags$style(
                      HTML(".scrollable-box {
             max-height: 600px;
             overflow-y: auto;
             overflow-x: hidden;
           }
           .scrollable-box pre {
             white-space: pre-wrap;
             word-wrap: break-word;
           }")
                    )
                  ),
                  box(
                    title = "VCF information",
                    width = NULL,
                    status = "primary",
                    solidHeader = TRUE,
                    collapsible = TRUE,
                    class = "scrollable-box",
                    height = "650px",
                    verbatimTextOutput("long_text")
                  )
                  
                  ),
        tabItem(tabName = "TransposonesTab",
                wellPanel(style = "background: #e9edef",
                  fluidRow(
                  column(6,selectInput("FusionT", "Select an individal Fusion:", c(TOTAL$Fusion),selected = TOTAL$Fusion[1] )),
                  column(6,submitButton("Search"))
                  )
        ),
        plotOutput("Transposones1"),
        plotOutput("Transposones2")
        ),
        tabItem(tabName = "breakpoint",
                wellPanel(style = "background: #e9edef",
                  fluidRow(
                  column(6,selectInput("FusionB", "Select an individal Fusion:", c(TOTAL$Fusion),selected = TOTAL$Fusion[1] )),
                  # column(6,radioButtons("GeneB", "Gene selection:", choiceNames = c("Gene 1","Gene 2"), choiceValues = c(1,2))),
                  column(6,submitButton("Search"))
                  )
        ),
        tags$head(
          tags$style(
            HTML(".scrollable-boxb {
             max-height: 950px;
             overflow-y: auto;
             overflow-x: auto;
           }
           .scrollable-box pre {
             white-space: pre-wrap;
             word-wrap: break-word;
           }")
          )
        ),
        box(
          title = "ARRIBA || Breakpoint Reads in Alignment GENE 1:",
          width = NULL,
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          class = "scrollable-boxb",
          height = "500px",
          htmlOutput("reads_text1")
        ),tags$head(
          tags$style(
            HTML(".scrollable-boxbc {
             max-height: 500px;
             overflow-y: auto;
             overflow-x: auto;
           }
           .scrollable-box pre {
             white-space: pre-wrap;
             word-wrap: break-word;
           }")
          )
        ),
        box(
          title = "ARRIBA ||Breakpoint Reads in Alignment GENE 2:",
          width = NULL,
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          class = "scrollable-boxbc",
          height = "450px",
          htmlOutput("reads_text2")
        ),
        tags$head(
          tags$style(
            HTML(".scrollable-boxb {
             max-height: 450px;
             overflow-y: auto;
             overflow-x: auto;
           }
           .scrollable-box pre {
             white-space: pre-wrap;
             word-wrap: break-word;
           }")
          )
        ),
        box(
          title = "STARFUSION || Breakpoint Reads in Alignment GENE 1:",
          width = NULL,
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          class = "scrollable-boxb",
          height = "1000px",
          htmlOutput("reads_text3")
        ),tags$head(
          tags$style(
            HTML(".scrollable-boxbc {
             max-height: 950px;
             overflow-y: auto;
             overflow-x: auto;
           }
           .scrollable-box pre {
             white-space: pre-wrap;
             word-wrap: break-word;
           }")
          )
        ),
        box(
          title = "STARFUSION ||Breakpoint Reads in Alignment GENE 2:",
          width = NULL,
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          class = "scrollable-boxbc",
          height = "1000px",
          htmlOutput("reads_text4")
        )
        
        )

        )
      )

)

















server <- function(input, output) {

  observeEvent(input$toggle_input, {
    if (input$toggle_input %% 2 == 1) {
      # Show input row
      shinyjs::show("input_row")
    } else {
      # Hide input row
      shinyjs::hide("input_row")
    }
  })
    output$filtered <- renderDT({
      rank_ass <- c()
      rank_p <- c()
      xlist <- check_mitelman(dataM$Fusion,input$DiseaseT, input$TopologyT)
      for(i in 1:length(xlist)){
        rank_ass[i] <- sum(xlist[[i]])
        rank_p[i] <- paste(xlist[i], collapse = ":")
      }
      dataM$ranking <- rank_ass
      dataM$rank_esp <- rank_p
      dataM$onco <- check_onco(dataM$Gene1,dataM$Gene2)
      dataM$Pseudo1 <- pseudogenes$P1
      dataM$Pseudo2 <- pseudogenes$P2
      dataM$Repetitive1 <- Repetitiv$R1
      dataM$Repetitive2 <- Repetitiv$R2
      
      dataM <- dataM[dataM$ranking > 0,]
      print(length(dataM$Fusion))
      dataM <- dataM[dataM$ranking >= input$Rank[1],]
      print(length(dataM$Fusion))
      if (input$Oncogene == "BOTH"){
        dataM <- dataM[dataM$onco ==2,]
      }
      else if (input$Oncogene == "ONE"){
        dataM <- dataM[dataM$onco == 1,]
      }
      dataM <- filter(dataM, spanning.pairs >= input$SpanningP[1])
      dataM <- filter(dataM,spanning.pairs <= input$SpanningP[2])
      dataM <- filter(dataM,spanning.reads >= input$SpanningR[1])
      dataM <- filter(dataM, spanning.reads <= input$SpanningR[2])
      dataM["Quality"] <- dataM$onco
      dataM$Quality <- ifelse(dataM$Pseudo1 == "Gene" & dataM$Pseudo2 == "Gene", dataM$Quality + 1, dataM$Quality)
      dataM$Quality <- ifelse(dataM$Repetitive1 == "Normal region" & dataM$Repetitive2 == "Normal region", dataM$Quality + 1, dataM$Quality)
      dataM$Quality <- ifelse(dataM$Proven1 > 0 | dataM$Proven2 > 0,  dataM$Quality + 1, dataM$Quality)
      dataM$Quality <- ifelse(dataM$Callers > 1 , dataM$Quality + 1, dataM$Quality)
      dataM['Priority'] = dataM['Quality'] + dataM['ranking']
      write.csv(as.matrix(dataM), "KNOWN_CALLS.csv", row.names = FALSE)
      desired_order <- c("Fusion", "base1", "chrom1", "base2", "chrom2", "Callers", "ranking", "onco", "Proven1", "Proven2", "Pseudo1", "Pseudo2", "Repetitive1", "Repetitive2", "Quality", "Priority","Gene1", "Gene2", "strand1", "strand2", "spanning.pairs", "spanning.reads", "GeneId1", "GeneId2", "caller")
      dataM <- dataM[,desired_order]
      datatable(dataM,options = list(scrollX = TRUE,
                                     scrollXInner = "100%")) %>% formatStyle(
        'ranking',
        backgroundColor = styleEqual(c(0, 1, 2, 3 , 4), c('red', 'orange', 'yellow', 'lightgreen', 'green'))) %>% formatStyle(
          'onco',
          backgroundColor = styleEqual(c(0, 1, 2), c('white','cyan' ,'blue'))) %>%     formatStyle(
            'Pseudo1',
            backgroundColor = styleEqual(c("Gene","Pseudogene"), c('lightgrey','darkgrey'))) %>%     formatStyle(
              'Pseudo2',
              backgroundColor = styleEqual(c("Gene","Pseudogene"), c('lightgrey','darkgrey'))) %>%     formatStyle(
                'Repetitive1',
                backgroundColor = styleEqual(c("Normal region","Repeated region"), c('lightgrey','darkgrey'))) %>%     formatStyle(
                  'Repetitive2',
                  backgroundColor = styleEqual(c("Normal region","Repeated region"), c('lightgrey','darkgrey'))) 

      })
    
    
    output$selected <- renderDT({
      rank_ass <- c()
      rank_p <- c()
      xlist <- check_mitelman(dataM$Fusion,"All", "All")
      for(i in 1:length(xlist)){
        rank_ass[i] <- sum(xlist[[i]])
        rank_p[i] <- xlist[i]
      }
      dataM$ranking <- rank_ass
      dataM$rank_esp <- rank_p
      dataM$onco <- check_onco(dataM$Gene1,dataM$Gene2)
      dataM$Pseudo1 <- pseudogenes$P1
      dataM$Pseudo2 <- pseudogenes$P2
      dataM$Repetitive1 <- Repetitiv$R1
      dataM$Repetitive2 <- Repetitiv$R2
      #print(input$Fusion)
      dataM <- dataM[dataM$Fusion %in% input$FusionIn,]
      dataM["Quality"] <- dataM$onco
      dataM$Quality <- ifelse(dataM$Pseudo1 == "Gene" & dataM$Pseudo2 == "Gene", dataM$Quality + 1, dataM$Quality)
      dataM$Quality <- ifelse(dataM$Repetitive1 == "Normal region" & dataM$Repetitive2 == "Normal region", dataM$Quality + 1, dataM$Quality)
      dataM$Quality <- ifelse(dataM$Proven1 > 0 | dataM$Proven2 > 0,  dataM$Quality + 1, dataM$Quality)
      dataM['Priority'] = dataM['Quality'] + dataM['ranking']
      desired_order <- c("Fusion", "base1", "chrom1", "base2", "chrom2", "Callers", "ranking", "onco", "Proven1", "Proven2", "Pseudo1", "Pseudo2", "Repetitive1", "Repetitive2", "Quality", "Priority","Gene1", "Gene2", "strand1", "strand2", "spanning.pairs", "spanning.reads", "GeneId1", "GeneId2", "caller")
      dataM <- dataM[,desired_order]
      datatable(dataM,options = list(scrollX = TRUE,
                                     scrollXInner = "100%")) %>% formatStyle(
        'ranking',
        backgroundColor = styleEqual(c(0, 1, 2, 3 , 4), c('red', 'orange', 'yellow', 'lightgreen', 'green'))) %>% formatStyle(
          'onco',
          backgroundColor = styleEqual(c(0, 1, 2), c('white','cyan' ,'blue'))) %>%     formatStyle(
            'Pseudo1',
            backgroundColor = styleEqual(c("Gene","Pseudogene"), c('lightgrey','darkgrey'))) %>%     formatStyle(
              'Pseudo2',
              backgroundColor = styleEqual(c("Gene","Pseudogene"), c('lightgrey','darkgrey'))) %>%     formatStyle(
                'Repetitive1',
                backgroundColor = styleEqual(c("Normal region","Repeated region"), c('lightgrey','darkgrey'))) %>%     formatStyle(
                  'Repetitive2',
                  backgroundColor = styleEqual(c("Normal region","Repeated region"), c('lightgrey','darkgrey'))) 
    })
    
    output$long_text <- renderText({
      dataM <- dataM[dataM$Fusion %in% input$FusionIn,]
      ch1 <- dataM$chrom1
      bk1 <- dataM$base1
      ch2 <- dataM$chrom2
      bk2 <- dataM$base2
      
      results1 <- vcf_search(vcf_txt,bk1, ch1, input$ModType, input$Modification, input$TransChrom, input$depth)
      results2 <- vcf_search(vcf_txt,bk2, ch2, input$ModType, input$Modification,input$TransChrom, input$depth)
      #
      text = "Gene 1:\n"
      c = 1
      for (it in results1$Information){
        text <-paste(text,paste0("Result #", c ), "\n")
        text <- paste(text,paste0("\n", it),"\n")
        c = c +1
      }
      text <- paste(text,"Gene 2: \n" , "\n")
      c = 1
      for (it in results2$Information){
        text <-paste(text,paste0("Result #", c ), "\n")
        text <- paste(text,paste0("\n", it),"\n")
        c = c+1
      }
      text
    })
    output$Transposones1 <- renderPlot({
      plot_transposons(alignment_df,input$FusionT, "G1")
    })
    output$Transposones2 <- renderPlot({
      plot_transposons(alignment_df,input$FusionT, "G2")
    })
    
    output$filteredNovels <- renderDT({
      rank_ass <- c()
      rank_p <- c()
      xlist <- check_mitelman(dataM$Fusion,"All", "All")
      for(i in 1:length(xlist)){
        rank_ass[i] <- sum(xlist[[i]])
        rank_p[i] <- paste(xlist[[i]], colapse = ";")
      }
      dataM$ranking <- rank_ass
      dataM$onco <- check_onco(dataM$Gene1,dataM$Gene2)
      
      dataM$Pseudo1 <- pseudogenes$P1
      dataM$Pseudo2 <- pseudogenes$P2
      dataM$Repetitive1 <- Repetitiv$R1
      dataM$Repetitive2 <- Repetitiv$R2
      
      dataM <- dataM[dataM$ranking == 0,]
      if (input$OncogeneN == "BOTH"){
        dataM <- dataM[dataM$onco ==2,]
      }
      else if (input$OncogeneN == "ONE"){
        dataM <- dataM[dataM$onco == 1,]
      }
      dataM <- filter(dataM, spanning.pairs >= input$SpanningPN[1])
      dataM <- filter(dataM,spanning.pairs <= input$SpanningPN[2])
      dataM <- filter(dataM,spanning.reads >= input$SpanningRN[1])
      dataM <- filter(dataM, spanning.reads <= input$SpanningRN[2])
      dataM["Quality"] <- dataM$onco
      dataM$Quality <- ifelse(dataM$Pseudo1 == "Gene" & dataM$Pseudo2 == "Gene", dataM$Quality + 1, dataM$Quality)
      dataM$Quality <- ifelse(dataM$Repetitive1 == "Normal region" & dataM$Repetitive2 == "Normal region", dataM$Quality + 1, dataM$Quality)
      dataM$Quality <- ifelse(dataM$Proven1 > 0 | dataM$Proven2 > 0,  dataM$Quality + 1, dataM$Quality)
      dataM['Priority'] = dataM['Quality'] + dataM['ranking']
      desired_order <- c("Fusion", "base1", "chrom1", "base2", "chrom2", "Callers", "ranking", "onco", "Proven1", "Proven2", "Pseudo1", "Pseudo2", "Repetitive1", "Repetitive2", "Quality", "Priority","Gene1", "Gene2", "strand1", "strand2", "spanning.pairs", "spanning.reads", "GeneId1", "GeneId2", "caller")
      dataM <- dataM[,desired_order]
      write.csv(as.matrix(dataM), "NOVEL_CALLS.csv", row.names = FALSE)
      datatable(dataM,options = list(scrollX = TRUE,
                                     scrollXInner = "100%")) %>% formatStyle(
                                       'ranking',
                                       backgroundColor = styleEqual(c(0, 1, 2, 3 , 4), c('red', 'orange', 'yellow', 'lightgreen', 'green'))) %>% formatStyle(
                                         'onco',
                                         backgroundColor = styleEqual(c(0, 1, 2), c('white','cyan' ,'blue'))) %>%     formatStyle(
                                           'Pseudo1',
                                           backgroundColor = styleEqual(c("Gene","Pseudogene"), c('lightgrey','darkgrey'))) %>%     formatStyle(
                                             'Pseudo2',
                                             backgroundColor = styleEqual(c("Gene","Pseudogene"), c('lightgrey','darkgrey'))) %>%     formatStyle(
                                               'Repetitive1',
                                               backgroundColor = styleEqual(c("Normal region","Repeated region"), c('lightgrey','darkgrey'))) %>%     formatStyle(
                                                 'Repetitive2',
                                                 backgroundColor = styleEqual(c("Normal region","Repeated region"), c('lightgrey','darkgrey'))) 


    })
    
    output$reads_text1 <- renderText({
      dataM <- dataM[dataM$Fusion %in% input$FusionB,]
      setwd("/home/sapatri/")
      print(paste("description:",which(dataM$Fusion == input$FusionB,arr.ind=TRUE)))
      text1 = paste(readLines(paste0("1.2practiques/BREAK/Fusions0/", paste0(which(TOTAL$Fusion == input$FusionB,arr.ind=TRUE)-1, "_1.csv"))),collapse="\n")

      
      text1 <- gsub(" ", "<span style='font-family: monospace;'>&nbsp;</span>", text1)
      # Add line breaks after each newline
      text1 <- gsub("\n", "<br>", text1)
      
      HTML( gsub("I", "<span style='font-family: monospace;'>I</span>",
        gsub("E", "<span style='font-family: monospace;'>E</span>",
        gsub("Fusion","<strong>Fusion</strong>",
              gsub("A", "<span style='color: red;font-family: monospace;'>A</span>", 
                      gsub("C", "<span style='color: green;font-family: monospace;'>C</span>", 
                           gsub("G", "<span style='color: blue;font-family: monospace;'>G</span>", 
                                gsub("T", "<span style='color: orange;font-family: monospace;'>T</span>", text1))))))))

    })
    output$reads_text2 <- renderText({
      dataM <- dataM[dataM$Fusion %in% input$FusionB,]
      setwd("/home/sapatri/")
      print(paste("description:",which(dataM$Fusion == input$FusionB,arr.ind=TRUE)))
      text2 = paste(readLines(paste0("1.2practiques/BREAK/Fusions0/", paste0(which(TOTAL$Fusion == input$FusionB,arr.ind=TRUE)-1, "_2.csv"))),collapse="\n")
      text2 <- gsub(" ", "<span style='font-family: monospace;'>&nbsp;</span>", text2)
      
      # Add line breaks after each newline
      text2 <- gsub("\n", "<br>", text2)
      
      HTML( gsub("I", "<span style='font-family: monospace;'>I</span>",
                 gsub("E", "<span style='font-family: monospace;'>E</span>",
                      gsub("Fusion","<strong>Fusion</strong>",
                           gsub("A", "<span style='color: red;font-family: monospace;'>A</span>", 
                                gsub("C", "<span style='color: green;font-family: monospace;'>C</span>", 
                                     gsub("G", "<span style='color: blue;font-family: monospace;'>G</span>", 
                                          gsub("T", "<span style='color: orange;font-family: monospace;'>T</span>", text2))))))))
    })
    
    
    output$reads_text3 <- renderText({
      dataM <- dataM[dataM$Fusion %in% input$FusionB,]
      setwd("/home/sapatri/")
      print(paste("description:",which(dataM$Fusion == input$FusionB,arr.ind=TRUE)))
      text3 = paste(readLines(paste0("1.2practiques/BREAK/Fusions2/", paste0(which(TOTAL$Fusion == input$FusionB,arr.ind=TRUE)-1, "_1.csv"))),collapse="\n")
      
      
      text3 <- gsub(" ", "<span style='font-family: monospace;'>&nbsp;</span>", text3)
      # Add line breaks after each newline
      text3 <- gsub("\n", "<br>", text3)
      
      HTML( gsub("I", "<span style='font-family: monospace;'>I</span>",
                 gsub("E", "<span style='font-family: monospace;'>E</span>",
                      gsub("Fusion","<strong>Fusion</strong>",
                           gsub("A", "<span style='color: red;font-family: monospace;'>A</span>", 
                                gsub("C", "<span style='color: green;font-family: monospace;'>C</span>", 
                                     gsub("G", "<span style='color: blue;font-family: monospace;'>G</span>", 
                                          gsub("T", "<span style='color: orange;font-family: monospace;'>T</span>", text3))))))))
      
    })
    output$reads_text4 <- renderText({
      dataM <- dataM[dataM$Fusion %in% input$FusionB,]
      setwd("/home/sapatri/")
      print(paste("description:",which(dataM$Fusion == input$FusionB,arr.ind=TRUE)))
      text4 = paste(readLines(paste0("1.2practiques/BREAK/Fusions2/", paste0(which(TOTAL$Fusion == input$FusionB,arr.ind=TRUE)-1, "_2.csv"))),collapse="\n")
      text4 <- gsub(" ", "<span style='font-family: monospace;'>&nbsp;</span>", text4)
      
      # Add line breaks after each newline
      text4 <- gsub("\n", "<br>", text4)
      
      HTML( gsub("I", "<span style='font-family: monospace;'>I</span>",
                 gsub("E", "<span style='font-family: monospace;'>E</span>",
                      gsub("Fusion","<strong>Fusion</strong>",
                           gsub("A", "<span style='color: red;font-family: monospace;'>A</span>", 
                                gsub("C", "<span style='color: green;font-family: monospace;'>C</span>", 
                                     gsub("G", "<span style='color: blue;font-family: monospace;'>G</span>", 
                                          gsub("T", "<span style='color: orange;font-family: monospace;'>T</span>", text4))))))))
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
