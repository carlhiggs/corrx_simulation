library(shiny)

# Define UI for slider demo application
shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("Power to detect a difference in Pearson correlations in Mz and Dz twins"),
  
  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(
    
    tags$h3("Scenario"),
            
    sliderInput("maxn", "Maximum N:", 
                min=30, max=1000, value=100, step=1),
    
    sliderInput("mzdz", "Ratio of Mz to Dz twin:", 
                min=0, max=10, value=1, step=0.01),
     
    sliderInput("r1", "Corr(x,y) in Mz twins:", 
                min=-1, max=1, value=0.6, step=0.01),

    sliderInput("r2", "Corr(x,y) in Dz twins:", 
                min=-1, max=1, value=0.2, step=0.01)
    
    # sliderInput("ordinate", "Normal ordinate:", 
    #             min=0, max=3, value=1.96, step=0.01)
    
    # selectInput("method", "Correlation method:",
    #                         c("Pearson"  = "pearson",
    #                           "Spearman" = "spearman")), 
    

  ),
  # Show a table summarizing the values entered
  mainPanel(
    plotOutput("graph1"),
    plotOutput("graph2"),
    tableOutput("datatable")
  )
 ))