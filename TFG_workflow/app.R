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
         #### Fisher ####
         tags$h3(paste0('Parameters for Fisher scoring')),
         numericInput('fisher_vars', 'Amount of variables to keep:', min = 1, value = 5000),
         actionButton('submit_fisher', 'Apply Fisher filter'),
         
         tags$hr(),
         #### ReliefF ####
         tags$h3(paste0('Parameters for ReliefF')),
         checkboxInput('use_relieff', 'Perform a ReliefF step?', value = TRUE),
         numericInput('relieff_vars', 'Amount of variables to keep:', min = 1, value = 1000),
         selectInput('relieff_iters', 'Amount of iterations to perform:', c('Dataset size' = 0, 'ln(Dataset size)' = -1, 'sqrt(Dataset size)' = -2), selected = 0),
         selectInput('relieff_method', 'ReliefF algorithm to calculate the scores:', c('ReliefFequalK', 'ReliefFexpRank'), selected = 'ReliefFequalK'),
         
         tags$hr(),
         #### PCA ####
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
         span(textOutput('error_text'), style='color:red'),
         textOutput('waiting_text'),
         plotOutput('fisher_plot')
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize = 100*1024^2)
  input_error = TRUE
  
  #### LOAD THE DATA AND PROCESS THE FIRST DATA ####
  observeEvent(input$dataset_file, {
    output$waiting_text = renderText('Please wait for the first calculations to finish...')
    
    dataset_file = input$dataset_file
    # First let's obtain the input data and check if it's well introduced
    if (endsWith(dataset_file$name, ".csv")){
      csv_data = read.csv(file = dataset_file$datapath, header = TRUE, as.is = TRUE)
      if(!('classes' %in% colnames(csv_data))){
        output$waiting_text = NULL
        output$error_text = renderText('ERROR: There is no classes column in the dataset! Please upload a valid dataset.')
      } else{
        input_error = FALSE
        output$error_text = NULL
        classes = csv_data$classes
        classes_position = match('classes', colnames(csv_data))
        data_matrix = as.matrix.data.frame(csv_data[-classes_position])
        positives = which(classes == classes[1])
        negatives = which(classes != classes[1])
      }
    } else if (endsWith(dataset_file$name, ".RData")){
      dataset = load(dataset_file$datapath)
      if(!('classes' %in% dataset)){
        output$waiting_text = NULL
        output$error_text = renderText('ERROR: There is no classes array in the dataset! Please upload a valid dataset.')
      } else if (length(dataset) != 2){
        output$waiting_text = NULL
        output$error_text = renderText('ERROR: This dataset includes a wrong amount of data! Please upload a valid dataset.')
      } else {
        input_error = FALSE
        output$error_text = NULL
        data_matrix = as.matrix.data.frame(get(dataset[3 - match('classes', dataset)]))
        positives = which(classes == classes[1])
        negatives = which(classes != classes[1])
      }
    } else {
      output$waiting_text = NULL
      output$error_text = renderText('ERROR: The file you uploaded is not valid. Make sure you upload a csv or RData file.')
    }
    
    if (!(input_error)){
      after_fisher = apply_fisher(data_matrix, positives, negatives, features_to_keep = input$fisher_vars, debug = FALSE)
      output$fisher_plot = renderPlot({
        plot(after_fisher$to_plot, ylab = 'Score')
        abline(v = input$fisher_vars, col = 'red', lwd = 2)
        })
    }
  })
  
  #### TO CHANGE THE INDICATOR LINE IN THE FISHER CHART ####
  observeEvent(input$fisher_vars, {
    if (!(input_error)){
      if (exists('after_fisher')){
        output$fisher_plot = renderPlot({
          plot(after_fisher$to_plot, ylab = 'Score')
          abline(v = input$fisher_vars, col = 'red', lwd = 2)
        })
      }
    }
  })
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

