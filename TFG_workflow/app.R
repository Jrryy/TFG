#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

if (!'rmarkdown' %in% installed.packages()){
  install.packages('rmarkdown')
}
library(rmarkdown)

if (!'png' %in% installed.packages()){
  install.packages('png')
}
library(png)

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
         numericInput('relieff_k', 'Amount of nearest neighbors to consider', min = 1, value = 10),
         selectInput('relieff_iters', 'Amount of iterations to perform:', c('Dataset size' = 0, 'ln(Dataset size)' = -1, 'sqrt(Dataset size)' = -2), selected = 0),
         selectInput('relieff_method', 'ReliefF algorithm to calculate the scores:', c('ReliefFequalK', 'ReliefFexpRank'), selected = 'ReliefFequalK'),
         checkboxInput('draw_relieff_heatmap', 'Draw heatmap of the selected variables after ReliefF?', value = TRUE),
         actionButton('submit_relieff', 'Apply ReliefF'),
         
         tags$hr(),
         #### PCA ####
         tags$h3(paste0('Parameters for sPLS-DA')),
         numericInput('components', 'Amount of PCA components to start calculating with', value = 5, min = 2),
         numericInput('cv_folds', 'Amount of cross validation folds to perform', value = 5, min = 1),
         numericInput('cv_repeats', 'Amount of cross validation repetitions to perform', value = 10, min = 1),
         checkboxInput('draw_indiv_plot', 'Draw a representation of the samples in the two first components?', value = TRUE),
         checkboxInput('draw_auroc', 'Draw ROC graphs?', value = FALSE),
         checkboxInput('draw_cim', 'Draw clustered image graph?', value = TRUE),
         checkboxInput('draw_loadings', 'Draw the contribution of each selected variable?', value = TRUE),
         actionButton('submit_splsda', 'Apply sPLS-DA'),
         
         tags$hr(),
         #### LAST BUTTON ####
         downloadButton('download', 'Save results as pdf')
         ),
      
      # Show a plot of the generated distribution
      mainPanel(
         conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                          tags$img(src = "loader.gif", style = 'position: absolute;'),
                          style = 'margin-left: auto; margin-right: auto; width: 50%;'),
         span(textOutput('error_text'), style='color:red'),
         plotOutput('pca_plot'),
         plotOutput('fisher_plot'),
         plotOutput('relieff_plot'),
         plotOutput('relieff_heatmap'),
         plotOutput('error_rate_plot'),
         plotOutput('tuning_plot'),
         plotOutput('individues_plot'),
         plotOutput('auroc_plot'),
         imageOutput('cim', height = '800px'),
         plotOutput('loadings')
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, clientData, session) {
  options(shiny.maxRequestSize = 100*1024^2)
  
  input_error = reactiveVal(TRUE)
  original_pca_data = reactiveVal(NULL)
  original_fisher_data = reactiveVal(NULL)
  after_fisher = reactiveVal(NULL)
  after_fisher_pca_data = reactiveVal(NULL)
  after_relieff = reactiveVal(NULL)
  after_relieff_pca_data = reactiveVal(NULL)
  after_plsda_perf = reactiveVal(NULL)
  after_splsda = reactiveVal(NULL)
  classes = reactiveVal(NULL)
  data_matrix = reactiveVal(NULL)
  positives = reactiveVal(NULL)
  negatives = reactiveVal(NULL)
  
  #### LOAD THE DATA AND PROCESS THE FIRST DATA ####
  observeEvent(input$dataset_file, {
    output$relieff_plot = NULL
    # output$waiting_text = renderText('Please wait for the first calculations to finish...')
    
    dataset_file = input$dataset_file
    # First let's obtain the input data and check if it's well introduced
    if (endsWith(dataset_file$name, ".csv")){
      csv_data = read.csv(file = dataset_file$datapath, header = TRUE, as.is = TRUE)
      if(!('classes' %in% colnames(csv_data))){
        # output$waiting_text = NULL
        output$error_text = renderText('ERROR: There is no classes column in the dataset! Please upload a valid dataset.')
        input_error(TRUE)
      } else{
        input_error(FALSE)
        output$error_text = NULL
        classes(csv_data$classes)
        classes_position = match('classes', colnames(csv_data))
        data_matrix(as.matrix.data.frame(csv_data[-classes_position]))
        classes_names = levels(factor(classes()))
        positives(which(classes == classes_names[1]))
        negatives(which(classes == classes_names[2]))
      }
    } else if (endsWith(dataset_file$name, ".RData")){
      dataset = load(dataset_file$datapath)
      if(!('classes' %in% dataset)){
        # output$waiting_text = NULL
        output$error_text = renderText('ERROR: There is no classes array in the dataset! Please upload a valid dataset.')
        input_error(TRUE)
      } else if (length(dataset) != 2){
        # output$waiting_text = NULL
        output$error_text = renderText('ERROR: This dataset includes a wrong amount of data! Please upload a valid dataset.')
        input_error(TRUE)
      } else {
        input_error(FALSE)
        output$error_text = NULL
        classes_index = match('classes', dataset)
        classes(get(dataset[classes_index]))
        data_matrix(as.matrix.data.frame(get(dataset[3 - classes_index])))
        classes_names = levels(factor(classes()))
        positives(which(classes() == classes_names[1]))
        negatives(which(classes() == classes_names[2]))
      }
    } else {
      # output$waiting_text = NULL
      output$error_text = renderText('ERROR: The file you uploaded is not valid. Make sure you upload a csv or RData file.')
      input_error(TRUE)
    }
    
    if (!(input_error())){
      after_splsda(NULL)
      updateNumericInput(session, 'fisher_vars', value = floor(ncol(data_matrix())/10))
      updateNumericInput(session, 'relieff_vars', value = floor(ncol(data_matrix())/100))
      original_pca_data(apply_pca(data_matrix(), classes(), debug = FALSE))
      output$pca_plot = renderPlot(plot(original_pca_data(), main = 'Principal components with all variables'))
      original_fisher_data(apply_fisher(data_matrix(), positives(), negatives(), features_to_keep = input$fisher_vars, debug = FALSE))
      # output$waiting_text = renderText('Drawing Fisher scores plot...')
      output$fisher_plot = renderPlot({
        plot(original_fisher_data()$to_plot, ylab = 'Score', main = 'Fisher scores')
        abline(v = input$fisher_vars, col = 'red', lwd = 2)
      })
    }
  })
  
  #### TO CHANGE THE INDICATOR LINE IN THE FISHER CHART ####
  # observeEvent(input$fisher_vars, {
  #   if (!(exists('input_error'))){
  #     # output$waiting_text = renderText('Drawing Fisher scores plot...')
  #     if (!(is.null(after_fisher))){
  #       output$fisher_plot = renderPlot({
  #         plot(after_fisher$to_plot, ylab = 'Score', main = 'Fisher scores')
  #         abline(v = input$fisher_vars, col = 'red', lwd = 2)
  #       })
  #       print('Am I here?')
  #       # output$waiting_text = NULL
  #     }
  #   }
  # })
  
  #### CALCULATE RELIEFF ####
  observeEvent(input$submit_fisher, {
    if (!(input_error())){
      after_splsda(NULL)
      after_fisher(apply_fisher(data_matrix(), positives(), negatives(), features_to_keep = input$fisher_vars, debug = FALSE))
      after_fisher_pca_data(apply_pca(after_fisher()$data, classes(), debug = FALSE))
      output$pca_plot = renderPlot(plot(after_fisher_pca_data(), main = 'Principal components with the selected variables after Fisher scoring'))
      if (input$use_relieff){
        after_relieff(apply_relieff(after_fisher()$data, classes(), input$relieff_k, features_to_keep = input$relieff_vars, iterations = input$relieff_iters, estimator = input$relieff_method, debug = FALSE))
        output$relieff_plot = renderPlot({
          plot(after_relieff()$sorted_attrs, ylab = 'Score', main = 'ReliefF scores')
          abline(v = input$relieff_vars, col = 'red', lwd = 2)
        })
        output$relieff_heatmap = NULL
        if (input$draw_relieff_heatmap){
          output$relieff_heatmap = renderPlot(heatmap(after_relieff()$data, main = paste('Heatmap of the', input$relieff_vars, 'variables after applying ReliefF')))
        } else {
          output$relieff_heatmap = NULL
        }
        output$error_text = NULL
      } else {
        after_relieff(after_fisher())
        output$relieff_plot = NULL
        output$error_test = NULL
      }
    }
  })
  
  #### CALCULATE PLS-DA PERF ####
  observeEvent(input$submit_relieff, {
    if (!(input_error())){
      # If there aren't errors in the input:
      if (is.null(after_relieff()) && input$use_relieff){
        # ERROR if there's no after_relieff data
        output$error_text = renderText('Please apply the Fisher filter before trying to apply ReliefF to your data.')
      } else {
        if (input$use_relieff){
          # If we want to use relieff, calculate all the related data
          if (input$fisher_vars < input$relieff_vars){
            output$error_text = renderText('ERROR: You can not keep more variables than you have. Review the amount of ReliefF variables.')
          } else{
            after_splsda(NULL)
            output$error_text = NULL
            after_relieff_pca_data(apply_pca(after_relieff()$data, classes(), debug = FALSE))
            output$pca_plot = renderPlot(plot(after_relieff_pca_data(), main = 'Principal components with the selected variables after ReliefF'))
            after_plsda_perf(apply_plsda_perf(after_relieff()$data, classes(), components = input$components, cv_folds = input$cv_folds, cv_repeats = input$cv_repeats, debug = FALSE))
          }
        } else {
          # Else, we will have to replace the relieff data with the data from Fisher
          after_plsda_perf(apply_plsda_perf(after_fisher()$data, classes(), components = input$components, cv_folds = input$cv_folds, cv_repeats = input$cv_repeats, debug = FALSE))
        }
    	  output$error_rate_plot = renderPlot({
    		  matplot(after_plsda_perf()$perf_plsda$error.rate$BER, type = 'l', lty = 1, col = color.mixo(1:3), main = 'Balanced Error Rate for amount of components', ylab = 'Balanced Error Rate')
    		  legend('topright', c('max.dist', 'centroids.dist', 'mahalanobis.dist'), lty = 1, col = color.mixo(1:3))
    	  })
        output$tuning_plot = renderPlot(plot(after_plsda_perf()$tune_splsda, col = color.jet(input$components), main = 'Error rates for number of components'))
      }
    }
  })
  
  #### CALCULATE sPLS-DA ####
  observeEvent(input$submit_splsda, {
    if (!(input_error())){
      if (is.null(after_plsda_perf())){
        output$error_text = renderText('Please apply the ReliefF step before trying to apply sPLS-DA.')
      } else {
        output$error_text = NULL
        final_components = after_plsda_perf()$final_ncomp
        features_to_keep = after_plsda_perf()$features_to_keep
        after_splsda(apply_splsda(after_relieff()$data, classes(), features_to_keep, components = final_components, cv_folds = input$cv_folds, cv_repeats = input$cv_repeats, debug = FALSE))
        if (input$draw_indiv_plot){
          output$individues_plot = renderPlot(plotIndiv(after_splsda()$splsda, comp = c(1, 2), ind.names = FALSE, legend = TRUE, ellipse = TRUE, title = 'sPLS-DA, final result, components 1 and 2'))
        } else {
          output$individues_plot = NULL
        }
        if (input$draw_auroc){
          output$auroc_plot = renderPlot(auroc(after_splsda()$splsda, roc.comp = 1))
        } else {
          output$auroc_plot = NULL
        }
        conds = levels(factor(classes()))
        cond.col = c('red', 'green')
        names(cond.col) = conds
        png("cim.png", height = 800, width = 800)
        cim(after_splsda()$splsda, row.sideColors = cond.col[classes()], legend = list())
        dev.off()
        if (input$draw_cim){
          output$cim = renderImage({
            list(src="cim.png")
            },
            deleteFile = FALSE)
        } else {
          output$cim = NULL
        }
        if (input$draw_loadings){
          output$loadings = renderPlot(plotLoadings(after_splsda()$splsda, contrib = 'max', method = 'mean'))
        } else {
          output$loading = NULL
        }
      }
    }
  })
  
  output$download = downloadHandler(
    filename = 'results.pdf',
    content = function(file) {
      params = list(input = input,
                    original_data = data_matrix(),
                    classes = classes(),
                    positives = positives(),
                    negatives = negatives(),
                    original_pca_data = original_pca_data(),
                    original_fisher_data = original_fisher_data(),
                    after_fisher_pca_data = after_fisher_pca_data(),
                    after_fisher_data = after_fisher(),
                    after_relieff_pca_data = after_relieff_pca_data(),
                    after_relieff_data = after_relieff(),
                    after_plsda_perf = after_plsda_perf(),
                    after_splsda = after_splsda())
      rmarkdown::render("../report.Rmd", output_file = file,
                        params = params, envir = new.env(parent = globalenv()))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

