#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

##########################
######## Packages ########
##########################

library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)

##########################
######## Functions #######
##########################

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


##########################
######## Files ########
##########################

accessibility_peaks_data <- loadRData("data/accessibility_peaks_data.rda")
UMI_mean_filtred_data <- loadRData("data/UMI_mean_filtred_data.rda")

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel(title = "Alert MARS-ATAC on chromosomes !",
    windowTitle = "Un dauphin"),
  
  h3("App description"),
  
  wellPanel(
    p("This is my firt ever made shiny.app.
           Please feel free to use it and report any bugs at romualdparmentier@orange.fr",
      style = 'font-size:18px')
  ),
  
  h3("Notice"),
  
  wellPanel(
    p("- First select your chromosome of choice, your time points and your donor and then click on the update button",
      style = 'font-size:18px'),
    p("- You can actively choose you window of interest once your input parameters have been set",
      style = 'font-size:18px'),
  ),
  
  h3("Input parameters"),
  
  wellPanel(
    fluidRow(
      column(2,
        selectInput(
          inputId = "donor",
          label = "Donors",
          choices = list("Donor 1" = 1,
            "Donor 2" = 2),
          selected = 1),
        
        selectInput(inputId = "chr",
          label = "Chromosome",
          choices = list("chr1" = "chr1", "chr2" = "chr2", "chr3" = "chr3", "chr4" = "chr4", "chr5" = "chr5", "chr6" = "chr6", "chr7" = "chr7", "chr8" = "chr8", "chr9" = "chr9", "chr10" = "chr10",
            "chr11" = "chr11", "chr12" = "chr12", "chr13" = "chr13", "chr14" = "chr14", "chr15" = "chr15", "chr16" = "chr16", "chr17" = "chr17", "chr18" = "chr18", "chr19" = "chr19", "chr20" = "chr20",
            "chr21" = "chr21", "chr22" = "chr22")
        )
      ),
      
      column(3,
        checkboxGroupInput(
          inputId = "MARS_time",
          label = "MARS Days",
          choices = list("00Hrs" = "00Hrs",
            "24Hrs" = "24Hrs",
            "48Hrs" = "48Hrs",
            "72Hrs" = "72Hrs",
            "96Hrs" = "96Hrs"),
          selected = "00Hrs",
          inline = TRUE),
        
        checkboxGroupInput(
          inputId = "ATAC_time",
          label = "ATAC Days",
          choices = list("00Hrs" = "00h",
            "24Hrs" = "24h",
            "48Hrs" = "48h"),
          selected = "00h",
          inline = TRUE),
      ),
      
      column(2,
        selectInput(inputId = "peaks_type",
          label = "ATAC peaks type",
          choices = list(TSS_1kb  = "TSS_1kb",
                      Intergenic = "intergenic",
                      CTCF = "CTCF",
                      Intergenic_CTCF = "CTCF_intergenic"),
        selected = "TSS_1kb" )
        
        ),
      
      column(2,
        actionButton(inputId = "update_button",
          label = "Click to update the plot",
          style='padding:50px 100px; font-size:200%')
      )
    )
  ),
  
  h3("Window of interest"),
  
  wellPanel(
    fluidRow(
      
      column(12,
        sliderInput(inputId = "position",
          label = "Position in chromosome in Mbp",
          min = 0,
          max = 250,
          step = 1,
          value = c(20,100),
          sep = TRUE,
          ticks = TRUE),
      )
    )
  ),
  
  mainPanel(width = 12,
    h3("MARS-seq plot :"),
    plotlyOutput(outputId = "MARS_plot"),
    h3("ATAC-seq plot :"),
    plotlyOutput(outputId = "ATAC_plot"),
    h3("MARSATAC plot :"),
    plotlyOutput(outputId = "MARSATAC_plot")
  )
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  UMI_data <- eventReactive(input$update_button,{UMI_mean_filtred_data %>%
      filter(donor == input$donor) %>%
      filter(chr == input$chr) %>%
      filter(condition %in% input$MARS_time)}
  )
  
  peaks_data <- eventReactive(input$update_button,{accessibility_peaks_data %>%
      filter(time %in% input$ATAC_time) %>%
      filter(chr == input$chr)%>%
      filter(type %in% input$peaks_type) %>%
      filter(yes_no == TRUE)}
  )
  
  output$MARS_plot <- renderPlotly({
    
    hover_text = paste0("Gene : ", UMI_data()$transcript_name_chr,
                        "\nNb cell : ", UMI_data()$sum_cell,
                        "\nVariance : ", round(UMI_data()$variance,1))
    
    plot <- ggplot() +
      geom_point(aes(x = UMI_data()$start_position, y = UMI_data()$avg_log2_UMI, 
                     color = UMI_data()$condition, text = hover_text),
                 shape = 1, size = 2) +
      xlim(input$position[1]*1e+06 , input$position[2]*1e+06) +
      labs(x = "chr position (bp)", y = "means(log2_UMI_sum) among cells", color = "MARS-seq time points")+
      ggtitle(paste("Chromosome =", input$chr, "| donor =", input$donor)) +
      theme(plot.title = element_text(size = 15, face = "bold"))
    
    ggplotly(plot, tooltip = "text") %>%
      config(displaylogo = FALSE) %>%
      layout(yaxis=list(fixedrange=TRUE)) %>%
      config(modeBarButtons = list(list("resetScale2d"), list("zoomOut2d"),list("zoomIn2d"), list("toImage")))


  })
  
  output$ATAC_plot <- renderPlotly({
    
    plot <- ggplot() +
      geom_rect(data = peaks_data() ,
        aes(xmin = peaks_data()$start, xmax=peaks_data()$end, ymin=1.5, ymax=max(UMI_data()$avg_log2_UMI), fill=peaks_data()$time),
        color = "grey",
        alpha = 0.8) +
      xlim(input$position[1]*1e+06 , input$position[2]*1e+06) +
      labs(x = "chr position (bp)", y = "means(log2_UMI_sum) among cells", fill = "ATAC-seq peak time :") +
      ggtitle(paste("Chromosome =", input$chr, "| donor =", input$donor, "| type of peak =", input$peaks_type)) +
      theme(plot.title = element_text(size = 15, face = "bold"))
    
    ggplotly(plot) %>% 
      config(displaylogo = FALSE) %>%
      layout(yaxis=list(fixedrange=TRUE)) %>%
      config(modeBarButtons = list(list("resetScale2d"), list("zoomOut2d"),list("zoomIn2d"), list("toImage")))

  })
  
  output$MARSATAC_plot <- renderPlotly({
    
   hover_text = paste0("Gene : ", UMI_data()$transcript_name_chr,
                        "\nNb cell : ", UMI_data()$sum_cell,
                        "\nVariance : ", round(UMI_data()$variance,1))
    
    plot <- ggplot() +
      geom_rect(data = peaks_data() ,
        aes(xmin = peaks_data()$start, xmax=peaks_data()$end, ymin=1.5, ymax=max(UMI_data()$avg_log2_UMI), fill=peaks_data()$time),
          color = "grey",
        alpha = 0.8) +
      xlim(input$position[1]*1e+06 , input$position[2]*1e+06) +
      labs(x = "chr position (bp)", y = "means(log2_UMI_sum) among cells", fill = "ATAC-seq peak time :")+
      # labs(x = "chr position (bp)", y = "means(log2_UMI_sum) among cells", fill = "ATAC-seq time points", color = "MARS-seq time points")+
      geom_point(data = UMI_data(),
        aes(x=UMI_data()$start_position, y=UMI_data()$avg_log2_UMI,
            color = UMI_data()$condition, text = hover_text), 
        shape = 1, size = 2) +  
      ggtitle(paste("Chromosome =", input$chr, "| donor =", input$donor, "| type of peak =", input$peaks_type)) +
      theme(plot.title = element_text(size = 15, face = "bold"))
    
    ggplotly(plot, tooltip = "text") %>%
     config(displaylogo = FALSE) %>%
      layout(yaxis=list(fixedrange=TRUE)) %>%
      config(modeBarButtons = list(list("resetScale2d"), list("zoomOut2d"),list("zoomIn2d"), list("toImage")))

    
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
