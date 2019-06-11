#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

source('../workflow.R')

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Workflow data"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         fileInput('dataset_file',
                   'Dataset file'),
         tags$hr(),
         tags$h3(paste0('Parameters for Fisher scoring')),
         numericInput('fisher_vars', 'Amount of variables to keep:', min = 0, value = 5000),
         
         tags$hr(),
         
         tags$h3(paste0('Parameters for ReliefF')),
         checkboxInput('use_relieff', 'Perform a ReliefF step?', value = TRUE),
         numericInput('relieff_vars', 'Amount of variables to keep:', min = 0, value = 1000),
         selectInput('relieff_iters', 'Amount of iterations to perform:', c('Dataset size' = 0, 'ln(Dataset size)' = -1, 'sqrt(Dataset size)' = -2), selected = 0),
         selectInput('relieff_method', 'ReliefF algorithm to calculate the scores:', c('ReliefFequalK', 'ReliefFexpRank'), selected = 'ReliefFequalK'),
         
         tags$hr(),
         
         tags$h3(paste0('Parameters for sPLS-DA')),
         numericInput('components', 'Amount of PCA components to extract', value = 10, min = 2),
         numericInput('cv_folds', 'Amount of cross validation folds to perform', value = 5, min = 1),
         numericInput('cv_repeats', 'Amount of cross validation repetitions to perform', value = 10, min = 1),
         
         tags$hr(),
         #### LAST BUTTON ####
         actionButton('submit', 'Submit')
         ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput('fisher_plot'),
         textOutput('waiting_text'),
         textOutput('error_text')
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize = 100*1024^2)
  
  observeEvent(input$submit, {
    output$waiting_text = renderText('Please wait for the calculations to finish...')
    
    dataset_file = input$dataset_file
    dataset = load(dataset_file$datapath)
    if (!('classes' %in% dataset)){
      output$waiting_text = NULL
      output$error_text = renderText('ERROR: There is no classes array in the dataset! Please upload a valid dataset.')
      print(length(dataset))
    } else if (length(dataset) != 2){
      output$waiting_text = NULL
      output$error_text = renderText('ERROR: This dataset includes a wrong amount of data! Please upload a valid dataset.')
    } else {
      data_matrix = as.matrix.data.frame(get(dataset[3 - match('classes', dataset)]))
      positives = which(classes == classes[1])
      negatives = which(classes != classes[1])
      
      after_fisher = apply_fisher(data_matrix, positives, negatives, features_to_keep = input$fisher_vars, debug = FALSE)
      
      output$fisher_plot = renderPlot(plot(after_fisher$to_plot, ylab = 'Score'))
    }
  })
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

