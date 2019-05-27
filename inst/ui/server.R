if(FALSE){
  
  library(shiny)
  runApp()
  
}


library(bmdModeling)
library(shinysky) # for hotable
library(rmarkdown) # for automatic report
library(shinyjs) #for toggle

`%then%` <- shiny:::`%OR%`
track <- FALSE
dataDir <- system.file("extdata", package = "bmdModeling")

jsCode <- "shinyjs.pageCol = function(params){$('body').css('background', params);}"


serverFunction <- function(input, output, session){
  
  observe({
        
        if (is.null(input$debug_console))
          return(NULL)
        
        if (input$debug_console > 0) {
          
          options(browserNLdisabled = TRUE)
          saved_console <- ".RDuetConsole"
          if (file.exists(saved_console)) {load(saved_console)}
          isolate(browser())
          save(file = saved_console, list = ls(environment()))
          
        }
        
      })
  
  
  results <- reactiveValues()
  
  
#  output$downloadSink <- downloadHandler(
#      filename = paste0("consoleOutput.txt"),
#      content = function(file) {
#        file.copy(from = file.path(tempdir(), "consoleOutput.txt") , to = file)
#      })
  
  
  
  
  ## 1. Load Data ##
  
  results$initialData <- reactive({
        
        validate(need(input$sep != input$dec, 
                "Comma for both data separator and decimal separator is not allowed"))
        
#        if(input$example1 > 0) {
#          
#          load(file.path(dataDir, "das10.rda"))
#          return(das10)          
#          
#        } else if(input$example2 > 0){
#          
#          return(f.scan(file.path(dataDir, "methyleug.txt")))
#          
#        } else {
        
        inFile <- input$dataLoaded
        
        if (is.null(inFile))
          return(NULL)
        
        rawData <- read.table(inFile$datapath, header = TRUE, 
            sep = input$sep, quote = '"', stringsAsFactors = FALSE, 
            dec = input$dec)
        
        if (all(sapply(rawData, class) == "character"))
          dataType <- "proastData" else
          dataType <- "rawData"
        
        
        tryCatch({ 
              
              if (dataType == "rawData") {
                
                validate(need(ncol(rawData) > 1, "Data were not loaded correctly") %then%
                        need(nrow(rawData) > 1, "Data were not loaded correctly"))
                
                
                list(info = substring(inFile$name , first = 1, last = nchar(inFile$name) - 4), 
                    nvar = ncol(rawData),
                    varnames = colnames(rawData),
                    data = rawData,
                    dtype = rep(0, ncol(rawData)))
                
              } else {
                
                # also returns list with info (title), nvar, varnames, data, dtype
                tmpData <- f.scan(inFile$datapath, separator = input$sep)
                
                if (input$dec == ",")
                  tmpData$data <- as.data.frame(sapply(tmpData$data, 
                          function(x) as.numeric(sub(",", "\\.", as.character(x)))))
                
                return(tmpData)
                
              } 
              
              
            }, error = function(err) {
              
              return(err)
              
            })
        
#        }
        
      })
  
  
  observe({
        
        inFile <- input$dataLoaded
        
        if(!is.null(inFile)){
          newName <- inFile$name 
          
          extension <- substring(newName, first = (nchar(newName) - 2), 
              last = nchar(newName))
          
          if(extension == "csv")
            updateRadioButtons(session, "sep", selected = ',')
          
          if(extension == "txt")
            updateRadioButtons(session, "sep", selected = '\t')
          
        }
        
#        if(input$example1 > 0){
#          
#          updateSelectInput(session, inputId = "dtype", selected = 10)
#          
#        } else if(input$example2 > 0){
#          
#          updateSelectInput(session, inputId = "dtype", selected = 4)
#          
#        }
        
      })
  
  
  # Choices are named numbers
  results$varnamesChoices <- reactive({
        
        newNames <- (results$initialData())$varnames
        values <- 1:length(newNames)
        names(values) <- newNames
        return(values)
        
      })  
  
  
  # Select subset of the data
  output$subsetVariable <- renderUI({
        
        selectInput("subsetVariable", "Subset of the data according to:",
            choices = c("<none>" = "none", results$varnamesChoices()), selected = "none")
        
      })
  
  
  output$levelSubset <- renderUI({
        
        validate(need(input$dataLoaded, "No data loaded") %then%
                need(input$subsetVariable, "Please select subset variable"))
        
        if (input$subsetVariable == "none") {
          
          return(NULL)
          
        } else {
          
          values <- results$initialData()$data[, as.numeric(input$subsetVariable)]
          
          if (is.factor(values))
            choices <- levels(values) else 
            choices <- unique(values)
          
          
          return(selectInput("levelSubset", "keep value", choices = choices, 
                  multiple = TRUE))
        }
        
      })
  
  
  # Take subset if necessary
  results$dataLoaded <- reactive({
        
        validate(need(results$initialData(), "Please load data"))
        
        if (is.null(input$subsetVariable))
          return(NULL)
        
        if (input$subsetVariable != "none" & is.null(input$levelSubset))
          return(NULL)
        
        if (is(results$initialData(), "error"))
          return(NULL)
        
        
        if (input$subsetVariable != "none") {
          
          subsetData <- results$initialData()$data[
              (results$initialData()$data[, as.numeric(input$subsetVariable)]
                    %in% input$levelSubset),]
          
          dataComplete <- results$initialData()
          dataComplete$data <- subsetData
          
          return(dataComplete)
          
        } else {
          
          return(results$initialData())
          
        }
        
      })
  
  
  output$selectedResponses <- renderUI({
        
        isResponse <- results$initialData()$dtype != 0        
        
        selectInput("responses", "Which response(s) do you want to consider?", 
            choices = results$varnamesChoices(), 
            selected = results$varnamesChoices()[isResponse], multiple = TRUE)
        
      })
  
  output$nonResponses <- renderUI({
        
        selectedNonResponses <- !(results$varnamesChoices() %in% input$responses)
        if (any(selectedNonResponses))
          selectedNames <- paste(names(results$varnamesChoices())[selectedNonResponses],
              collapse = ", ") else
          selectedNames <- "<none>"
        
        list(
            tags$b("List of non-responses:"), 
            selectedNames,
            tags$br(),
            tags$br()
        )
        
      })
  
  
  observe({
        
        orderedResponses <- as.numeric(input$responses)[
            order(as.numeric(input$responses))]
        
        results$responseNames <- reactive(results$varnamesChoices()[orderedResponses])
        results$nonResponseNames <- reactive(results$varnamesChoices()[
                !(results$varnamesChoices() %in% orderedResponses)])
        
      })
  
  
  output$dataLoaded <- DT::renderDataTable({
        
        validate(need( !any(attr(results$initialData(), "class") == "error"), 
                    paste("Data were not loaded correctly: \n", results$initialData()$message)) %then%
                need(results$dataLoaded(), "No data selected"))
        
        DT::datatable(data = results$dataLoaded()$data, options = list(
                lengthMenu = list(c(15, 100, -1), c('15', '100', 'All')),
                pageLength = 15), rownames = FALSE)
        
      })
  
  
  
  
  ## 2. Fit Models ##
  
  ## 2.1 Define variables ##
  
  output$parameterQuestions <- renderUI({
        
        validate(need(results$nonResponseNames(), 
                    "Loaded data contains only response variables. Please define non-response variables") %then%
                need(results$responseNames(), 
                    "Loaded data contains no response variables. Please define response variables"))
        
        list(
            
            selectInput("xans",
                label = "Independent variable (e.g. dose)",
                choices = results$nonResponseNames()),
            
            selectInput("Vyans",
                label = "Response variable(s)",
                choices = results$responseNames(), 
                selected = results$responseNames()[1], multiple = TRUE)
        
        )
        
      })
  
  results$selectedResponses <- reactive(results$varnamesChoices()[as.numeric(input$Vyans)])
  
  
  observe({
        
        allTypes <- results$initialData()$dtype
        allChoices <- c("continuous" = 1, "quantal" = 4, "binary" = 2, "ordinal" = 3)
        
        bestChoice <- allChoices[(allChoices %in% allTypes[allTypes != 0])][1]
        
        validate(need(bestChoice, ""))
        
        if (bestChoice == 1) {
          
          updateSelectInput(session, inputId = "dtype", selected = 1)
          updateSelectInput(session, inputId = "dtype2", selected = "individual")
          updateCheckboxInput(session, inputId = "isLitter", value = FALSE)
          
        } else if (bestChoice == 10) {
          
          updateSelectInput(session, inputId = "dtype", selected = 1)
          updateSelectInput(session, inputId = "dtype2", selected = "summary")
          updateCheckboxInput(session, inputId = "isLitter", value = FALSE)
          
        } else if (bestChoice == 5) {
          
          updateSelectInput(session, inputId = "dtype", selected = 1)
          updateSelectInput(session, inputId = "dtype2", selected = "individual")
          updateCheckboxInput(session, inputId = "isLitter", value = TRUE)
          
        } else if (bestChoice == 4) {
          
          updateSelectInput(session, inputId = "dtype", selected = 4)
          updateCheckboxInput(session, inputId = "isLitter", value = FALSE)
          
        } else if (bestChoice == 6) {
          
          updateSelectInput(session, inputId = "dtype", selected = 4)
          updateCheckboxInput(session, inputId = "isLitter", value = TRUE)
          
        } else {
          
          updateSelectInput(session, inputId = "dtype", selected = bestChoice)
          
        }
        
      })
  
  
  output$parameterQuestions2 <- renderUI({
        
        createVnans <- function(iSelected) {       
          
          Vnans <- lapply(seq_along(input$Vyans), function(i) {
                
                selectizeInput(paste0("nans", i),
                    label = names(results$selectedResponses())[i],
                    choices = results$nonResponseNames(),
                    selected = results$nonResponseNames()[iSelected])
                
              })
          
          return( list(
                  tags$b("Sample size for response variable:"),
                  column(11, Vnans, offset = 1) 
              )
          )
          
        }
        
        validate(need(input$dtype, ""), 
            need(input$dtype2, ""), 
            need(input$Vyans, ""))
        
        
        if (input$dtype == "1" & input$dtype2 == "summary") { # Continuous summary data
          list(
              selectInput("sans", 
                  label = "Variation statistic",
                  choices = results$nonResponseNames(), 
                  selected = results$nonResponseNames()[2]),
              radioButtons("sd.se",
                  label = "Type of variation statistic",
                  choices = c("standard deviations" = 1, "standard errors" = 2)),
              createVnans(iSelected = 3)
          
          )
          
          
        } else if (results$dtype() %in% c(4, 6)) {
          
          createVnans(iSelected = 2)
          
        } else if (results$dtype() == 5){
          
          column(11, selectInput("nest.no", "with nested factor",
                  choices = results$nonResponseNames(),
                  selected = results$nonResponseNames()[2]), offset = 1)
          
        }
        
      })
  
  
  results$cont <- reactive({
        
        results$dtype() %in% c(1, 5, 10, 15, 25, 250, 26, 260)
        
      })
  
  
  output$modelChoices <- renderUI({
        
        if (input$singleModel == "no")
          return(NULL)
        
        selectInput("selectedModel", "Model", 
            choices = getModelNames(dtype = results$dtype()))
        
      })
  
  output$blankLine <- renderUI({
        
        validate(need(FALSE, " "))
        
      })
  
  
  ## 2.2 Models to fit ##
  output$selectedModels <- renderUI({
        
        if (!results$dtype() %in% c(2, 4) | input$singleModel == "yes")
          return(NULL)
        
        fittedModels <- getModelNames(dtype = results$dtype())[-c(1,2)]        
        
        list(
            
            checkboxInput("performMA", " Perform model averaging", value = TRUE),
            
            conditionalPanel("input.performMA == true",
                
                actionLink("advancedMA", "Advanced settings", icon = icon("cog")),
                
                
#                        radioButtons(inputId = "defaultModelsMA",
#                            label = "Selected models for model averaging",
#                            choices = c("Default (converged models)" = "yes", 
#                                "Advanced" = "no")),
                
                fluidRow(
                    column(11, conditionalPanel("input.advancedMA % 2 == 1",  # Advanced settings
                            
                            list(
                                checkboxInput("doNaiveApproach", "Include results for averaged bmd (naive approach)", value = FALSE),
                                strong("Selected models for model averaging"),
                                lapply(seq_along(results$selectedResponses()), function(iResponse) {
                                      
                                      selectInput(inputId = paste0("modelsMA", iResponse), 
                                          label = names(results$selectedResponses())[iResponse], 
                                          choices = fittedModels,
                                          selected = fittedModels,
                                          multiple = TRUE)})
                            
                            )
                        
                        ), offset = 1)
                )
            )
        )
        
      })
  
  
  
  ## 2.3 General parameters ##
  
  output$ces.ans <- renderUI({
        
        if (is.null(results$dtype()))
          return(NULL)
        
        # for quantal response only
        if (results$dtype() %in% c(2, 4, 6)) {
          
          selectInput("ces.ans", label = "Benchmark criterion",  
              choices = c("ED50" = 1, 
                  "Additional risk, i.e. P[BMD] - P[0]" = 2,
                  "Extra risk, i.e. (P[BMD]-P[0])/(1-P[0])" = 3),
              #                  "CED for latent variable" = 4))
              selected = 3) 
          
          
        } else if (input$dtype == '3') {
          
          yNames <- names(results$varnamesChoices())[as.numeric(input$Vyans)]
          yValues <- results$dataModified()$data[, as.numeric(input$Vyans), drop = FALSE]
          nLevels <- apply(yValues, 2, max)
          
          allAnswers <- list()
          
          for (i in 1:length(input$Vyans)) {
            
            allAnswers[[i]] <- column(11, 
                selectizeInput(paste0("CES.cat", i),
                    label = yNames[i], choices = unique(yValues[,i])), 
                offset = 1)
            
          }
          
          list(
              tags$b("Severity category associated with CED"),
              allAnswers
          )
          
        }
        
      })
  
  
  output$ces <- renderUI({
        
        if (is.null(results$dtype()))
          return(NULL)
        
        # Non-quantal response
        if (!results$dtype() %in% c(2, 4, 6)) {
          return(
              list(
                  numericInput("CES", label = "Value for CES", value = 0.05),
                  actionLink(inputId = "helpCES", label = "CES/BMR",
                      icon = icon("info-circle")),
                  
                  conditionalPanel("input.helpCES % 2 == 1",
                      p(em("Critical Effect Size, i.e. the Benchmark Response (BMR), 
                                  defined as the estimated difference in response compared with the background response."))
                  ),
                  
                  conditionalPanel("input.helpCES % 2 == 0",
                      tags$br()
                  )
              )
          )
        }
        
        if(is.null(input$ces.ans))
          return(NULL)
        
        #Switch text depending on shinyInput$ces.ans:
        cesLabel <- switch(input$ces.ans, 
            '1' = "Value for CES (positive)", 
            '2' = "Value for the BMR, in terms of additional risk",
            '3' = "Value for the BMR, in terms of extra risk") 
#            '4' = "Value for the CES, defined for latent variable")
        
        if(input$ces.ans == "1") # No impact of CES value
          return(NULL) else       
          return(numericInput("CES", label = cesLabel, value = 0.10))
        
        
      })
  
  
  results$dataModified <- reactive({
        
        validate(need(results$dataLoaded(), "Please load data") %then%
                need(input$xans, "Please define variables"))
        
        # Remove missing values
        dataComplete <- f.remove.NAs(variableIndices = as.numeric(c(input$xans, 
                    input$Vyans, input$sans, input$nans)), 
            originalData = results$dataLoaded()$data, track = track)
        
        # index of the variables used for parameters (a, b, variance (theta), c, d)
        allFactors <- as.numeric(c(input$fct1.no, input$fct2.no, input$fct3.no, 
                input$fct4.no, input$fct5.no, input$factors))
        allFactors <- allFactors[allFactors != 0]
        
        if(length(allFactors) != 0){
          
          dataComplete <- f.remove.NAs(variableIndices = allFactors, 
              originalData = dataComplete, track = track)
          
        }
        
        # Check whether x and y values are numeric vectors
        
        validate(need(input$xans, "Please define independent variable") %then%
                need(input$Vyans, "Please define response variable(s)"))
        
        f.check.nonneg.num(dataFrame = dataComplete[, as.numeric(input$xans)], track = track)
        isProblem <- f.check.nonneg.num(dataFrame = dataComplete[, as.numeric(input$Vyans)], track = track)
        
        validate(need(!any(isProblem), 
                "No log-transformation can be performed due to negative response values."))
        
        
        list(info = (results$dataLoaded())$info, 
            nvar = (results$dataLoaded())$nvar,
            varnames = (results$dataLoaded())$varnames, 
            data = dataComplete, 
            dtype = (results$dataLoaded())$dtype)        
        
      })
  
  
  
  
  ## 2.4 Model parameters ##
  
  # TODO if loglik = NA, give warning and let choose other startValues?
  results$initialModel <- reactive({
        
        if (input$singleModel == "no")
          return(NULL)   
        
        validate(need(input$selectedModel, "No model selected"))
        
        initShinyInput <- results$shinyInput()
        initShinyInput$yans <- initShinyInput$Vyans[1]
        initShinyInput$Vyans <- NULL
        initShinyInput$nans <- initShinyInput$Vnans[1]
        initShinyInput$Vnans <- NULL
        initShinyInput$CES.cat <- initShinyInput$allCES.cat[1]
        
        initShinyInput$startValues <- NULL
        initShinyInput$parameterConstraints <- NULL
        
        initData <- initShinyInput$data
        initShinyInput$data <- NULL
        
        fittedModel <- fitSingleModel(data = initData, 
            shinyInput = initShinyInput,
            selectedModel = as.numeric(input$selectedModel),
            noFit = TRUE)
        
        validate(need(!is(fittedModel, "error"), 
                paste("Error: Please check whether you specified correct values for the input fields.",
                    fittedModel$message)))        
        
        return(fittedModel)
        
      })
  
  output$factors <- renderUI({
        
        selectInput("factors", label = "Covariates",
            choices = results$nonResponseNames(), multiple = TRUE)
        
      })
  
  results$allParameters <- reactive({
        
        if (results$cont()) {
          
          allParameters <- c("Background response parameter (a)" = 1, 
              "Potency parameter (BMD/CED)" = 2, "Variance (var)" = 3, 
              "Maximum response parameter (c)" = 4, 
              "'Steepness' parameter (d)" = 5)
          
        } else if (input$dtype == '3') {
          
          allParameters <- c("Background response parameter (a)" = 1, 
              "Potency parameter (BMD/CED)" = 2, 
              "Variance (var)" = 3, 
              "Maximum response parameter (c)" = 4, 
              "'Steepness' parameter (d)" = 5)
          
        } else {
          
          allParameters <- c("Background response parameter (a)" = 1, 
              "Potency parameter (BMD/CED)" = 2, 
              "'Steepness' parameter (c)" = 4)
          
        }
        
        
        return(allParameters)
        
      })
  
  
  output$includeCovariates <- renderUI({
        
        # Show only present parameters for single model; always include the third
#      if(!is.null(results$presentParameters()))
#        allParameters <- allParameters[sapply(allParameters, function(iParameter)
#                  any(grepl(names(allParameters)[iParameter], results$presentParameters()))
#            )]
        
        
        parameterFields <- lapply(1:5, function(iParameter) {
              
              if(!iParameter %in% results$allParameters())
                return(NULL)
              
              if (iParameter < 3) # Only for these parameters factors are selected by default
                selectInput(paste0("fct", iParameter, ".no"), 
                    names(results$allParameters())[iParameter],
                    choices = results$nonResponseNames(),
                    selected = results$varnamesChoices()[as.numeric(input$factors)],
                    multiple = TRUE) else
                selectInput(paste0("fct", iParameter, ".no"), 
                    names(results$allParameters())[results$allParameters() == iParameter],
                    choices = results$nonResponseNames(),
                    multiple = TRUE)
              
            })
        
        
        if (is.null(input$factors))
          textDefault <- NULL else
          textDefault <- em("By default no covariate is selected for", 
              paste(names(results$allParameters())[-c(1, 2)], collapse = ", "),
              "and the input field is blank. \n")
        
        
        return(
            list(
                textDefault,
                tags$br(),
                tags$b("Covariate with respect to"),
                column(11, parameterFields[[3]], offset = 1),
                column(11, list(tags$em("Scale parameters"), 
                        tags$br(), parameterFields[1:2]), offset = 1),
                column(11, list(tags$em("Shape parameters"),
                        tags$br(), parameterFields[4:5]), offset = 1)
            )
        )  
        
      })
  
  
  output$extraCovariates <- renderUI({
        
        if(!is.null(input$fct3.no) & !results$cont()) {
          
          numericInput("fct3.ref", "Reference value for study duration", 
              value = 0)
          
        } else {
          
          NULL
          
        }
        
      })
  
  observe({
        
        validate(need(results$initialModel(), 
                "Please wait... a single model is being fit"))
        
        presentParameters <- f.text.par(results$initialModel())
        
        validate(need(length(presentParameters) == length(results$initialModel()$lb), 
                "Please wait... a single model is being fit"),
            need(length(presentParameters) == length(results$initialModel()$par.start),
                "Please wait... a single model is being fit"))
        
        initDataFrame <- data.frame(parameterName = presentParameters,
            lowerBound = results$initialModel()$lb, 
            upperBound = results$initialModel()$ub,
            startValues = results$initialModel()$par.start, 
            stringsAsFactors = FALSE)
        
        validate(need(initDataFrame, "Default parameter constraints are unknown") %then%
                need(ncol(initDataFrame) == 4, 
                    "No valid result for initial parameter values and constraints"))
        
        colnames(initDataFrame) <- c("Parameter", "Lower bound", "Upper bound", "Start values")
        
        
        results$parameterValues <- initDataFrame
        
      })
  
  
  output$parameterValues <- renderHotable({
        
        results$parameterValues
        
      }, readOnly = c(TRUE, FALSE, FALSE, FALSE))
  
  
  output$setParameterValues <- renderUI({
        
        extraWarning <- NA
        # Check for one level factors & provide warning
        allFactors <- as.numeric(c(input$fct1.no, input$fct2.no, input$fct3.no, 
                input$fct4.no, input$fct5.no))
        allFactors <- unique(allFactors[allFactors != 0])
        
        if (length(allFactors) > 0) {
          
          oneLevel <- apply(results$dataModified()$data[, allFactors, drop = FALSE], 2, 
              function(x) length(unique(x)) == 1)
          
          if (any(oneLevel))
            extraWarning<- paste("The covariate(s)", 
                paste(names(results$varnamesChoices())[allFactors[oneLevel]], collapse = ","),
                "has/have only one level.")
          
        }
        
        validate(need(is.na(extraWarning), extraWarning),
            need(input$singleModel == "yes", 
                "Default settings for estimation of model parameters can only be changed when a single model is being fit")
        )
        
        modelChoices <- getModelNames(dtype = results$dtype())
        
        # For full model and continuous response: parameters are fixed (except first, var)
        validate(need(!is.na(results$cont()), ""), 
            need(!is.null(input$selectedModel), ""))
        if(results$cont() & as.numeric(input$selectedModel) == 2)
          warningText <- "For full model, all model parameters except 'var' are fixed."
        else warningText <- NULL
        
        list(
            tags$b("Estimation of model parameters:",
                names(modelChoices)[modelChoices == as.numeric(input$selectedModel)], "model"),
            warningText,
            tags$br(),
            hotable("parameterValues"),
            tags$em("Note: Define infinite constraints with Inf or -Inf. Define fixed parameters with equal lower and upper bound.")
        )
        
      })
  
  
  
  ## 3. Apply proast ##
  
  # Define proast parameters
  results$dtype <- reactive({
        
        validate(need(input$dtype, "Please select type of response in Data tab page"))
        
        dtype <- as.numeric(input$dtype)
        
        if (dtype == 1) {
          
          if (input$dtype2 == "summary") dtype <- 10
          if (input$isLitter) dtype <- 5
          
        } else if (dtype == 4) {
          
          if (input$isLitter) dtype <- 6
          
        }
        
        
#        # No log-transform of response
#        if(!is.null(input$transformation)){
#          
#          switch(input$transformation,
#              
#              "<none>" = {
#                
#                dtype[dtype == 1] <- 25
#                dtype[dtype == 10] <- 250
#                
#              },
#              
#              "square root" = {
#                
#                dtype[dtype == 1] <- 26
#                dtype[dtype == 10] <- 260
#                
#              })
#          
#        }
        
        return(dtype)
        
      })
  
  
  results$shinyInput <- reactive({
        
        if (is.null(input$CES))
          CES <- 0.05 else
          CES <- input$CES
        
        
        Vnans <- c()
        allCES.cat <- c()
        
        for (iName in 1:length(input$Vyans)){
          
          if(!is.null(input[[paste0("nans", iName)]]))
            Vnans[iName] <- as.numeric(input[[paste0("nans", iName)]])
          if(!is.null(input[[paste0("CES.cat", iName)]]))
            allCES.cat[iName] <- as.numeric(input[[paste0("CES.cat", iName)]])
          
        } 
        
        
        parameterConstraints <- NULL
        startValues <- NULL
        
        if (input$singleModel == "yes" & !is.null(input$parameterValues)) {
          
          parameterSettings <- hot.to.df(input$parameterValues)
          parameterConstraints <- data.frame(lowerBound = as.numeric(parameterSettings[, 2]),
              upperBound = as.numeric(parameterSettings[, 3]))
          startValues <- as.numeric(parameterSettings[, 4])
          
        }
        
        ## Define covariate for parameters
        fct1.no <- 0
        fct2.no <- 0
        fct3.no <- 0
        fct4.no <- 0
        fct5.no <- 0
        
        if (input$advancedSettings %% 2 == 0 & !is.null(input$factors)) { # No advanced settings
          
          fct1.no <- as.numeric(input$factors)
          fct2.no <- as.numeric(input$factors)
          
        } 
        
        if (input$advancedSettings %% 2 != 0) {
          
          if (!is.null(input$fct1.no))
            fct1.no <- as.numeric(input$fct1.no)
          
          if (!is.null(input$fct2.no))
            fct2.no <- as.numeric(input$fct2.no)
          
          if (!is.null(input$fct3.no))
            fct3.no <- as.numeric(input$fct3.no)
          
          if (!is.null(input$fct4.no))
            fct4.no <- as.numeric(input$fct4.no)
          
          if (!is.null(input$fct5.no))
            fct5.no <- as.numeric(input$fct5.no)
          
        }
        
        
        # Create multi-factor data and shinyInput
        currentData <- results$dataModified()$data
        lastIndex <- ncol(currentData)
        
        for (i in 1:5) {
          
          if (length(get(paste0("fct", i, ".no"))) > 1) {
            
            selectedFactors <- currentData[, get(paste0("fct", i, ".no"))]
            names(selectedFactors) <- names(currentData)[get(paste0("fct", i, ".no"))]
            currentData <- cbind(currentData, 
                interaction(selectedFactors, drop = TRUE, sep = "_"))
            names(currentData)[lastIndex + 1] <- paste(names(currentData)[get(paste0("fct", i, ".no"))], collapse = "_")
            assign(paste0("fct", i, ".no"), lastIndex + 1)
            lastIndex <- lastIndex + 1
            
          }
          
        }
        
        newData <- results$dataModified()
        newData$nvar <- ncol(currentData)
        newData$varnames <- names(currentData)
        newData$data <- currentData
        newData$dtype <- c(newData$dtype, 
            rep(0, ncol(currentData) - ncol(results$dataModified()$data)))
        
        
        inputValues <- c(
            
            list(
                data = newData,
                dtype = results$dtype(), 
                xans = as.numeric(input$xans), 
                Vyans = as.numeric(input$Vyans),
                CES = CES,
                cont = results$cont(),
                conf.lev = input$conf.lev,
                parameterConstraints = parameterConstraints,
                startValues = startValues,
                #include covariates
                fct1.no = fct1.no,
                fct2.no = fct2.no,
                fct3.no = fct3.no,
                fct4.no = fct4.no,
                fct5.no = fct5.no,
                fct3.ref = input$fct3.ref
            ),
            
            # nSamples not required for binary model
            if(results$dtype() != 2)
              list(Vnans = Vnans),
            
            # continuous summary models
            if(results$dtype() %in% c(10, 250, 260))
              list(
                  sans = as.numeric(input$sans), 
                  sd.se = as.numeric(input$sd.se)
              ),
            
            # type of benchmark response, quantal and binary models
            if (results$dtype() %in% c(2, 4, 6) && !is.null(input$ces.ans))
              list(ces.ans = as.numeric(input$ces.ans)),
            
            # severity category for each response y; ordinal data
            if (input$dtype == '3')
              list(allCES.cat = allCES.cat),
            
            # number of bootstrap runs for CI of CES with continuous clustered
            if (results$dtype() == 5)
              list(nruns = input$nBootstraps, nest.no = as.numeric(input$nest.no))
        
        )
        
        inputValues
        
      })
  
  
  ## 3.1  Fit proast models ##
  
  results$fittedModels <- eventReactive(input$submit, {
        
        withProgress(message = "Fitting Models", value = 0, {
              
#        sink(file.path(tempdir(), "consoleOutput.txt"))
#        cat("*** Fitting the models *** \n \n")
              
              returnObject <- lapply(seq_along(results$selectedResponses()), function(iResponse) {
                    
                    incProgress(1/(length(results$selectedResponses())), 
                        detail = paste0("for response '", names(results$selectedResponses()[iResponse]), "'"))      
                    
                    currentShinyInput <- results$shinyInput()
                    currentShinyInput$yans <- currentShinyInput$Vyans[iResponse]
                    currentShinyInput$Vyans <- NULL
                    currentShinyInput$nans <- currentShinyInput$Vnans[iResponse]
                    currentShinyInput$Vnans <- NULL
                    currentShinyInput$CES.cat <- currentShinyInput$allCES.cat[iResponse]
                    
                    currentData <- currentShinyInput$data
                    currentShinyInput$data <- NULL
                    
                    if (input$singleModel == "no") {
                      
                      fitAllModels(data = currentData, 
                          shinyInput = currentShinyInput, fitCovariateCombinations = TRUE)
                      
                    } else {
                      
                      fitSingleModel(data = currentData, 
                          shinyInput = currentShinyInput,
                          selectedModel = as.numeric(input$selectedModel))
                      
                    }
                    
                  })
              
#        sink()
              
              return(returnObject)
              
            })
        
      })
  
  
  # Create summary of all fitted models
  summarizeResults <- function(){
    
    if (!is.null(results$fittedModels())) {
      
      results$summaryTables <- reactive({
            
            lapply(results$fittedModels(), function(fittedModelsResponse) {
                  
                  if (input$singleModel == "yes") {
                    
                    validate(need( !any(is(fittedModelsResponse, "error")), 
                            paste("Error(s) in calculation: \n", 
                                fittedModelsResponse$message))) 
                    
                  } else {
                    
                    validate(need(!all(sapply(fittedModelsResponse, function(x) is(x, "error"))),
                            paste0("Error(s) in calculation: \n",
                                paste(sapply(fittedModelsResponse, function(x) x$message), collapse = "\n"))))
                    
                  }
                  
                  summaryResults <- summaryModels(savedResults = fittedModelsResponse)
                  
                  return(summaryResults)
                  
                })
            
          })
      
      
      sapply(seq_along(results$fittedModels()), function(iResponse) {
            
            summaryResult <- results$summaryTables()[[iResponse]]
            
            if (input$allCovariates) {
              
              tmpTable <- summaryResult$summaryTable
              
            } else {
              
              tmpTable <- filterBestCovariates(summaryResult$summaryTable)
              
            }
            
            colnames(tmpTable) <- attr(tmpTable, "columnNames")
            
            if (all(nchar(tmpTable[, "Included covariate(s)"]) < 1))
              tmpTable <- tmpTable[, names(tmpTable) != "Included covariate(s)"]
            
            responseName <- summaryResult$extraInfo$responseName
            
            
            output[[paste0("downloadTable", responseName)]] <- downloadHandler(
                filename = paste0("bmdTable_", responseName, ".csv"),
                content = function(file) {
                  write.csv(tmpTable, file)
                })
            
            stopAnalysis <- FALSE
            if (!is.null(summaryResult$extraInfo$aicWarning))
              if (grepl("null model", summaryResult$extraInfo$aicWarning))
                stopAnalysis <- TRUE
            
            
            if (input$singleModel == "no") {
              
              if (!stopAnalysis) {
                
                bestModelIndex <- summaryResult$extraInfo$bestModelIndex
                
                sapply(seq_along(bestModelIndex), function(iBest){
                      
                      output[[paste0("downloadPlot", responseName, iBest)]] <- downloadHandler(
                          filename = paste0("bmdPlot_", responseName, iBest, ".png"),
                          content = function(file) {
                            png(file)
                            f.plot.all(results$fittedModels()[[iResponse]][[bestModelIndex[iBest]]])
                            dev.off()
                          })
                    })
                
              } else {
                
                output[[paste0("downloadPlot", responseName)]] <- downloadHandler(
                    filename = paste0("bmdPlot_", responseName, ".png"),
                    content = function(file) {
                      png(file)
                      f.plot.all(results$fittedModels()[[iResponse]][[2]])
                      dev.off()
                    })
                
              }
              
            } else {
              
              output[[paste0("downloadPlot", responseName)]] <- downloadHandler(
                  filename = paste0("bmdPlot_", responseName, ".png"),
                  content = function(file) {
                    png(file)
                    f.plot.all(results$fittedModels()[[iResponse]])
                    dev.off()
                  })    
              
            }
            
          })
      
      
      results$warningText <- reactive({
            
            lapply(results$fittedModels(), function(fittedModelsResponse) {
                  
                  if (results$dtype() == 3)
                    paste0("Note: Estimated values for BMD, BMDL and BMDU are reported for the chosen severity category: ",
                        fittedModelsResponse[[length(fittedModelsResponse)]]$CES.cat, ".") else
                    NULL
                  
                  
                })
            
          })
      
      results$errorModels <- reactive({
            
            lapply(results$fittedModels(), function(fittedModelsResponse) {
                  
                  if (input$singleModel == "no") {
                    
                    errorModels <- lapply(fittedModelsResponse, function(iModel){
                          
                          if (is(iModel, "error")) 
                            tags$li(iModel$message) else 
                            NULL
                          
                        }) 
                    
                    if (!all(sapply(errorModels, is.null)))
                      return(tags$p("Error(s) in calculation:", 
                              tags$ul(unique(errorModels)))) else 
                      return(NULL)
                    
#                      # If none of the fitted models is without error
#                      validate(need(length(errorModels) != length(fittedModelsResponse),
#                              errorModels))
                    
                  }
                  
                })
            
          })
      
    }
    
  }
  
  
  
  ## 3.2 Model averaging ##
  
  results$averagedModels <- eventReactive(input$submit, {
        
        if (is.null(input$performMA))
          return(NULL)
        
        if (!input$performMA)
          return(NULL)
        
        if (!(results$dtype() %in% c(2, 4)) | input$singleModel == "yes")
          return(NULL)
        
        withProgress(message = "Model Averaging", value = 0, {
              
#        sink(file = file.path(tempdir(), "consoleOutput.txt"), append = TRUE)
#        cat("\n\n*** Model averaging *** \n\n")
              
              returnObject <- lapply(seq_along(results$selectedResponses()), function(iResponse) {
                    
                    incProgress(1/(length(results$selectedResponses())), 
                        detail = paste0("for response '", names(results$selectedResponses()[iResponse]), "'"))
                    
                    fittedModels <- getModelNames(dtype = results$dtype())
                    
                    tmpTable <- results$summaryTables()[[iResponse]]$summaryTable
                    aicNull <- tmpTable[which(tmpTable$model == "Null"), "aic"]
                    bestSummaryTable <- filterBestCovariates(summaryTable = tmpTable)
                    
                    # Default selected models for MA
                    if (input$advancedMA %% 2 == 0) {
                      
                      isAccepted <- bestSummaryTable$converged == "yes"
                      isFitted <- !(bestSummaryTable$model %in% c("Null", "Full", "Exp model 4"))
                      
                      selectedModels <- bestSummaryTable$model[isAccepted & isFitted]
                      
                      # Advanced selected models for MA  
                    } else {
                      
                      selectedModels <- names(fittedModels)[fittedModels %in% 
                              as.numeric(input[[paste0("modelsMA", iResponse)]])]
                      
                    }
                    
                    # No model averaging if no model is better than null model
                    if (!is.null(results$summaryTables()[[iResponse]]$extraInfo$aicWarning))
                      if (grepl(pattern = "null model", x = results$summaryTables()[[iResponse]]$extraInfo$aicWarning))
                        return(NULL)
                    
                    currentShinyInput <- results$shinyInput()
                    currentShinyInput$yans <- currentShinyInput$Vyans[iResponse]
                    currentShinyInput$Vyans <- NULL
                    currentShinyInput$nans <- currentShinyInput$Vnans[iResponse]
                    currentShinyInput$Vnans <- NULL
                    currentShinyInput$CES.cat <- currentShinyInput$allCES.cat[iResponse]
                    
                    currentData <- currentShinyInput$data
                    currentShinyInput$data <- NULL
                    
                    # Consider only selected models
                    selectedModelsOrdered <- bestSummaryTable$model[bestSummaryTable$model %in% selectedModels]
                    selectedParameterNames <- bestSummaryTable$parameterNames[
                        bestSummaryTable$model %in% selectedModels]
                    
                    selectedIndices <- matchWithResults(savedResults = results$fittedModels()[[iResponse]],
                        modelNames = selectedModelsOrdered, 
                        parameterNames = selectedParameterNames)
                    
                    modelResults <- results$fittedModels()[[iResponse]][selectedIndices]
                    
                    
                    # Calculate weights per response
                    weights <- calculateWeights(aicValues = bestSummaryTable$aic[
                            bestSummaryTable$model %in% selectedModels])
                    
                    validate(need(length(modelResults) == length(weights), 
                            "Number of fitted models does not equal number of weights"))
                    
                    
                    # Naive approach
                    if (input$doNaiveApproach) {
                      
                      maSimpleBmd <- optimizeBmd(weights = weights, 
                          modelResults = modelResults, 
                          naiveApproach = TRUE)
                      
                      withProgress(message = "Estimating naive model-averaged CIs", value = 0, {
                            
                            maSimpleBootstrap <- bootstrapBmd(proastData = currentData, 
                                weights = weights, 
                                modelResults = modelResults, 
                                shinyInput = currentShinyInput, 
                                aicNull = aicNull,
                                nBootstraps = input$nBootstraps,
                                naiveApproach = TRUE,
                                showProgress = TRUE)
                            
                          })
                      
                      maSimpleBounds <- findBootstrapBounds(bootstrapBmd = maSimpleBootstrap$bootstrapBmd,
                          confidenceLevel = currentShinyInput$conf.lev)
                      
                      bmdTable1 <- cbind(maSimpleBmd, maSimpleBounds)
                      colnames(bmdTable1) <- c("BMD", "BMDL", "BMDU")
                      
                    } else {
                      
                      bmdTable1 <- NULL
                      
                    } 
                    
                    # Non-naive appraoch
                    maBmd <- optimizeBmd(weights = weights, 
                        modelResults = modelResults)
                    
                    withProgress(message = "Estimating model-averaged CIs", value = 0, {
                          
                          maBootstrap <- bootstrapBmd(proastData = currentData, 
                              weights = weights, 
                              modelResults = modelResults,
                              shinyInput = currentShinyInput, 
                              aicNull = aicNull,
                              nBootstraps = input$nBootstraps,
                              showProgress = TRUE)
                          
                        })
                    
                    maBounds <- findBootstrapBounds(bootstrapBmd = maBootstrap$bootstrapBmd,
                        confidenceLevel = currentShinyInput$conf.lev)
                    
                    bmdTable2 <- cbind(maBmd, maBounds)
                    colnames(bmdTable2) <- c("BMD", "BMDL", "BMDU")
                    
                    bmdTable <- list(bmdTable1, bmdTable2)
                    attr(bmdTable, "responseName") <- results$summaryTables()[[
                        iResponse]]$extraInfo$responseName
                    
                    weightsData <- t(data.frame(weights))
                    colnames(weightsData) <- selectedModelsOrdered
                    rownames(weightsData) <- "Estimated model weights"
                    
                    
                    # Make plot
                    png(file.path(tempdir(), paste0("maBmdPlot_", iResponse, ".png")),
                        width = 600, height = 400)
                    plotAverageModel(proastData = currentData,
                        xans = currentShinyInput$xans, 
                        yans = currentShinyInput$yans, 
                        nans = currentShinyInput$nans, 
                        bmd = maBmd, 
                        modelResults = modelResults, 
                        bootstrapBmd = maBootstrap$bootstrapBmd, 
                        bootstrapModelResults = maBootstrap$modelResults,
                        confidenceLevel = currentShinyInput$conf.lev)
                    dev.off()
                    
                    
                    return(
                        list(bmdTable = bmdTable,
                            modelResults = modelResults,
                            bootstrapBmd = maBootstrap$bootstrapBmd,
                            bootstrapModelResults = maBootstrap$modelResults,
                            weights = weightsData
                        )
                    )
                    
                  })
              
#        sink()
              
              
              return(returnObject)
              
            })
        
      })
  
  
  # Summary of model-averaged results
  summarizeMA <- function() {
    
    sapply(seq_along(results$selectedResponses()), function(iResponse) {
          
          if (is.null(results$averagedModels()[[iResponse]]))
            return()
          
          iTables <- results$averagedModels()[[iResponse]]$bmdTable
          responseName <- attr(iTables, "responseName")
          bmdValues <- iTables[[2]]$BMD
          names(bmdValues) <- rownames(iTables[[2]])
          
          
          output[[paste0("downloadTableMA1", responseName)]] <- downloadHandler(
              filename = paste0("maBmdTable1_", responseName, ".csv"),
              content = function(file) {
                write.csv(iTables[[1]], file)
              })
          
          output[[paste0("downloadTableMA2", responseName)]] <- downloadHandler(
              filename = paste0("maBmdTable2_", responseName, ".csv"),
              content = function(file) {
                write.csv(iTables[[2]], file)
              })  
          
          output[[paste0("downloadPlotMA", responseName)]] <- downloadHandler(
              filename = paste0("maBmdPlot_", responseName, ".png"),
              content = function(file) {
                file.copy(file.path(tempdir(), paste0("maBmdPlot_", iResponse, ".png")), file)
              })   
          
          return(NULL)
          
        })
    
  }
  
  
  
  results$forestData <- reactive({
        
        tryCatch({
              
              # Remove previous plot
              if (file.exists(file.path(tempdir(), "forestPlot.png")))
                unlink(file.path(tempdir(), "forestPlot.png"))
              
              firstResult <- results$summaryTables()[[1]]
              
              # Results for MA are available
              if (is.null(input$performMA))
                doMA <- FALSE else if (!input$performMA)
                doMA <- FALSE else
                doMA <- TRUE
              
              if (input$singleModel == "yes")
                doMA <- FALSE
              
              
              if (doMA)
                covariateNames <- rownames(results$averagedModels()[[1]]$bmdTable[[2]]) else
                covariateNames <- sub(pattern = "bmd.", replacement = "", 
                    names(firstResult$summaryTable)[
                        grepl(pattern = "bmd.", names(firstResult$summaryTable), fixed = TRUE)])
              
              if (length(covariateNames) < 2)
                covariateNames <- NULL
              
              if (doMA) {
                
                bmdValues <- lapply(seq_along(results$selectedResponses()), function(iResponse) {
                      
                      iResult <- results$averagedModels()[[iResponse]]
                      
                      if (is.null(iResult))
                        return(NULL)
                      
                      data.frame(response = names(results$selectedResponses())[iResponse], 
                          bmd = t(iResult$bmdTable[[2]]$BMD),
                          bmdl = t(iResult$bmdTable[[2]]$BMDL), 
                          bmdu = t(iResult$bmdTable[[2]]$BMDU))
                      
                    })
                
              } else {
                
                bmdValues <- lapply(seq_along(results$selectedResponses()), function(iResponse) {
                      
                      iResult <- results$summaryTables()[[iResponse]]
                      
                      if (!is.null(iResult$extraInfo$aicWarning))
                        if (grepl("null model", iResult$extraInfo$aicWarning))
                          return(NULL)
                      
                      if (input$singleModel == "yes") {
                        
                        if (is.null(covariateNames)) {
                          
                          bmdValues <- data.frame(response = iResult$extraInfo$responseName,
                              bmd = iResult$summaryTable$bmd,
                              bmdl = iResult$summaryTable$bmdl, 
                              bmdu = iResult$summaryTable$bmdu)                    
                          
                        } else {
                          
                          bmdValues <- data.frame(response = iResult$extraInfo$responseName,
                              iResult$summaryTable[grepl("bmd.", names(iResult$summaryTable), fixed = TRUE)],
                              iResult$summaryTable[grepl("bmdl.", names(iResult$summaryTable), fixed = TRUE)], 
                              iResult$summaryTable[grepl("bmdu.", names(iResult$summaryTable), fixed = TRUE)])
                          
                        }
                        
                      } else {
                        
                        bmdValues <- data.frame(response = iResult$extraInfo$responseName,
                            as.list(iResult$extraInfo$minBmdl), 
                            as.list(iResult$extraInfo$maxBmdu))
                        
                        if (is.null(covariateNames))
                          names(bmdValues) <- c("response", "bmdl", "bmdu")
                        
                        bmdValues
                        
                      }
                      
                    })
                
              }
              
              bmdData <- do.call(rbind, bmdValues)
              if (any(is.na(bmdData)))
                stop("Missing values in estimated BMD values")
              
              return(list(
                      bmdData = bmdData,
                      covariateNames = covariateNames
                  ))
              
            }, error = function(err) 
              validate(need(FALSE, paste("Error in forestplot:", err$message)))
        )
        
      })
  
  
  # Make download button
  output$downloadForestPlot <- downloadHandler(
      filename = "forestPlot.png",
      content = function(file) {
        file.copy(file.path(tempdir(), "forestPlot.png"), file)
      })   
  
  
  
  ## 3.3 Results to show in UI ##
  
  output$summaryAnalysis <- renderUI({
        
        input$submit
        
        isolate({
              
              summarizeResults()
              
              if (!all(sapply(results$averagedModels(), is.null)))
                summarizeMA()
              
              # Create forest plot
              png(file.path(tempdir(), "forestPlot.png"),
                  width = 1000, height = 400)
              makeForestPlot(bmdData = results$forestData()$bmdData, 
                  groupNames = results$forestData()$covariateNames)
              dev.off()
              
              
              responseTabs <- lapply(seq_along(results$selectedResponses()), function(iResponse) {
                    
                    summaryResult <- results$summaryTables()[[iResponse]]
                    
                    if (input$allCovariates) {
                      
                      tmpTable <- summaryResult$summaryTable
                      
                    } else {
                      
                      tmpTable <- filterBestCovariates(summaryResult$summaryTable)
                      
                    }
                    
                    extraInfo <- summaryResult$extraInfo
                    colnames(tmpTable) <- attr(tmpTable, "columnNames")
                    
                    if (all(nchar(tmpTable[, "Included covariate(s)"]) < 1) & !is.null(input$factors)) {
                      
                      tmpTable <- tmpTable[, names(tmpTable) != "Included covariate(s)"]
                      warningCovariates <- "None of the models had a better fit when including the covariate(s)"
                      
                    } else warningCovariates <- NULL
                    
                    
                    # All fitted models
                    if (input$singleModel == "no") {
                      
                      bestModel <- apply(extraInfo$bestModel, 1, function(x){
                            if (x["parameterNames"] == "") 
                              x["model"] else
                              paste(x, collapse = " with covariate for ") 
                          })
                      bestModelIndex <- extraInfo$bestModelIndex
                      
                      
                      topResults <- list(
                          h4("Fitted Models"),
                          tags$em(results$errorModels()[[iResponse]]),
                          tags$em(results$warningText()[[iResponse]]),
                          tags$em(warningCovariates),
                          renderTable(tmpTable),
                          downloadButton(paste0("downloadTable", extraInfo$responseName), "Download")
                      )
                      
                      
                      # BMD without model averaging
                      if (is.null(results$averagedModels()[[iResponse]])) {
                        
                        stoppedAnalysis <- FALSE
                        if (!is.null(extraInfo$aicWarning))
                          if (grepl("null model", extraInfo$aicWarning))
                            stoppedAnalysis <- TRUE
                        
                        if (!stoppedAnalysis) {
                          
                          midResults <- list(
                              p("Lowest BMDL(s): ", 
                                  paste(extraInfo$minBmdl, collapse = "; ")),
                              p("Highest BMDU(s): ", 
                                  paste(extraInfo$maxBmdu, collapse = "; ")),
                              p("Best model(s): ", paste(bestModel, collapse = "; ")),
                              lapply(seq_along(bestModelIndex), function(iBest){
                                    list(
                                        renderPlot(f.plot.all(
                                                results$fittedModels()[[iResponse]][[bestModelIndex[iBest]]]),
                                            width = 600),
                                        downloadButton(paste0("downloadPlot", extraInfo$responseName, iBest), "Download")
                                    )
                                  })
                          )
                          
                        } else {
                          
                          midResults <- list(
                              tags$br(),
                              tags$br(),
                              renderPlot(f.plot.all(results$fittedModels()[[iResponse]][[2]]),
                                  width = 600),
                              downloadButton(paste0("downloadPlot", extraInfo$responseName), "Download")
                          )
                          
                        }
                        
                        bottomResults <- NULL
                        
                        
                        # BMD with model averaging
                      } else { 
                        
                        bmdTables <- results$averagedModels()[[iResponse]]$bmdTable
                        weightsTable <- results$averagedModels()[[iResponse]]$weights
                        
                        bmdValues <- bmdTables[[2]]$BMD
                        names(bmdValues) <- rownames(bmdTables[[2]])
                        
                        responseName <- attr(bmdTables, "responseName")
                        
                        if (input$doNaiveApproach) {
                          
                          printBmdTables <- fluidRow(
                              column(6,
                                  h5("Approach 1: Averaged bmd"),
                                  renderTable(bmdTables[[1]], rownames = TRUE),
                                  downloadButton(paste0("downloadTableMA1", responseName), "Download")
                              
                              ),
                              column(6, 
                                  h5("Approach 2: Averaged response"),
                                  renderTable(bmdTables[[2]], rownames = TRUE),
                                  downloadButton(paste0("downloadTableMA2", responseName), "Download")
                              )
                          )
                          
                        } else {
                          
                          printBmdTables <- list(
                              renderTable(bmdTables[[2]], rownames = TRUE),
                              downloadButton(paste0("downloadTableMA2", responseName), "Download")
                          )
                          
                        }
                        
                        
                        midResults <- NULL                       
                        bottomResults <- list(
                            
                            tags$br(),
                            tags$br(),
                            h4("Model Averaging"),
                            
                            fluidRow(
                                column(6, 
                                    renderTable(weightsTable, rownames = TRUE),
                                    printBmdTables
                                ),
                                
                                column(6, 
                                    renderImage({
                                          list(src = file.path(tempdir(), paste0("maBmdPlot_", iResponse, ".png")))
                                        }, deleteFile = FALSE),
                                    downloadButton(paste0("downloadPlotMA", responseName), "Download")
                                )
                            )
                        )
                        
                      } 
                      
                      tabPanel(title = extraInfo$responseName, 
                          list(
                              tags$br(),
                              if (!is.null(extraInfo$aicWarning))
                                tags$span(style = "color:red", 
                                    tags$p(tags$strong("Please consult a BMD specialist:"), extraInfo$aicWarning)),
                              # Results of all fitted models
                              topResults, midResults, bottomResults,
                              tags$hr()
                          )
                      )
                      
                      
                    } else {
                      
                      tabPanel(title = extraInfo$responseName, 
                          list(
                              tags$br(),
                              
                              tags$em(results$warningText()[[iResponse]]),
                              tags$em(warningCovariates),
                              renderTable(tmpTable),
                              results$fittedModels()[[iResponse]]$warningCovariates,
                              
                              renderPlot(f.plot.all(results$fittedModels()[[iResponse]]), width = 600), 
                              downloadButton(paste0("downloadPlot", extraInfo$responseName), "Download"),
                              tags$hr()                      
                          
                          )
                      )
                    }
                    
                  })
              
              myTabs <- c(responseTabs,
                  list(tabPanel(title = "Summary",
                          list(
                              tags$br(),
                              if (file.exists(file.path(tempdir(), "forestPlot.png")))
                                downloadButton("downloadForestPlot", "Download"),
                              tags$br(),
                              renderImage({
                                    list(src = file.path(tempdir(), "forestPlot.png"),
                                        alt = "No plot available")
                                  }, deleteFile = FALSE)                              
                          )
                      )
                  )
              )
              
              
              js$scrollTo(id = "#summaryAnalysis")
              
              
              do.call(tabsetPanel, myTabs)
              
              
            })
        
      })
  
  
  
# From proast manual:
#The log-likelihood profile is plotted, which was used to calculate the confidence interval for the CED. The horizontal line indicates the 90% (two-sided) confidence level, and its intersections with the log-likelihood profile provide the confidence limits. This plot is given as a check: it should be a smoothly increasing and then decreasing function. Irregularities in this function indicate estimation problems. 
  
# Show submit button
  output$fitModel <- renderUI({
        
        allWarnings <- list()
        
        validate(need(results$dataModified(), "Please load data and specify all parameters") %then%
                need(input$Vyans, "Please define response variable(s)"))
        
        doInvalidate <- results$shinyInput()
        
        # Need number of bootstrap runs for CI of BMD when dtype == 5
        if (results$dtype() %in% c(2, 4, 5))
          validate(need(input$nBootstraps, "Please provide numeric value for the number of bootstrap runs") %then%
                  need(input$nBootstraps > 0, "Please provide positive value for the number of bootstrap runs"))
        
        
        # Check start values for model parameters
        if(input$singleModel == "yes") {
          validate(need(results$initialModel(), "Initial model is being fit") %then%
                  need(!is(results$initialModel(), "error"), 
                      paste("Error: \n", results$initialModel()$message)))
          
          if(results$initialModel()$errorAdjustStartValues)
            allWarnings <- c(allWarnings, list( 
                    "\nWarning: Please choose other start values for the model parameters:
                        The log-likelihood could not be calculated."))
          
        }
        
        # Check range of dose variable
        rawX <- results$dataModified()$data[, as.numeric(input$xans)]
        
        if(is.factor(rawX)){
          
          allWarnings <- c(allWarnings, list("\nWarning: Please select a numeric independent variable."))
          
#        } else {
#          
#          xValues <- rawX/input$sf.x
#          xName <- names(results$varnamesChoices())[as.numeric(input$xans)]
#          xRange <- range(xValues[xValues != 0])
#          
#          if(xRange[2] > 100){
#            
#            allWarnings <- c(allWarnings, list(paste("\nWarning: (Nonzero)", xName, "values range from", xRange[1],
#                        "to", xRange[2], "\nit is recommended to scale",
#                        xName, "to prevent numerical problems.")))
#            
#          } 
        }
        
        # Check response trend and compare to input for CES
        if (results$dtype() %in% c(4, 6) & !is.null(results$shinyInput()$Vnans)) {
          
          responseValues <- results$dataModified()$data[, as.numeric(input$Vyans), drop = FALSE]
          nResponses <- results$dataModified()$data[, as.numeric(results$shinyInput()$Vnans), drop = FALSE]
          
          validate(need(ncol(responseValues) == ncol(nResponses), "Please define sample size for each response"))
          allRawY <- responseValues/nResponses
          
        } else {
          
          allRawY <- results$dataModified()$data[, as.numeric(input$Vyans), drop = FALSE]
          
        }
        
        allWarnings <- tryCatch({
              
              coefficients <- sapply(allRawY, function(rawY) lm(rawY ~ rawX)[[1]])
              maxResponse <- coefficients[1,] + coefficients[2,] * max(rawX)
              baseResponse <- coefficients[1,]
              
              if (results$dtype() %in% c(2, 4, 6) & !is.null(input$ces.ans)) {
                
                bmrFunction <- paste0("bmr", input$ces.ans)
                
                if (bmrFunction == "bmr1")
                  CES <- 0.50 else
                  CES <- input$CES
                
                maxCes <- do.call(bmrFunction, list(baseResponse = baseResponse, 
                        maxResponse = maxResponse)) * sign(coefficients[2,])
                
                if (all(coefficients[2,] < 0))
                  allWarnings <- c(allWarnings, list("Note: Response values are expected to increase monotonically.",
                          "It is recommended to convert the response values (e.g. 1 - response values) before performing the analysis."))
                
                if (any(maxCes < input$CES))
                  allWarnings <- c(allWarnings, list("Note: Given the data, we expect the BMD to exceed the range of observed dose values.",
                          "You may want to decrease the value for the BMR."))
                
              } else {
                
                difference <- maxResponse - baseResponse
                
                isTrue <- ifelse(coefficients[2,] > 0,
                    (1 + input$CES) * baseResponse < maxResponse,
                    (1 - input$CES) * baseResponse > maxResponse)
                
                if (!all(isTrue))
                  allWarnings <- c(allWarnings, list("Given the data, the value for CES might not be attained.",
                          "You may want to decrease the value for CES."))          
                
              }
              
              allWarnings
              
            }, error = function(err){
              
              allWarnings <- c(allWarnings, list("The data trend could not be assessed.",
                      "You may want to adapt the input fields in 'Define variables'.")) 
              
              return(allWarnings)
              
            })
        
        
        # Number of rows excluded due to missing values
        rowsMissing <- nrow(results$dataLoaded()$data) - nrow(results$dataModified()$data)
        
        if(rowsMissing > 0){
          
          allWarnings <- c(allWarnings, list(paste("\nWarning: There are", rowsMissing, 
                      "rows with missing values removed from the data.")))
        }
        
        
        if(length(allWarnings) > 0)
          allWarnings <- c(allWarnings, list(tags$br(), tags$br()))
        
        
        list(
            
            allWarnings,
            actionButton(inputId = "submit", label = "Fit Model(s)", 
                styleclass = "primary", size = "mini")
        
        )
        
      })
  
  
  
  ## 4. Make report ##
  
  output$downloadReport <- downloadHandler(
      filename = 'reportBmd.docx',
      
      content = function(file) {
        
        tmpDir <- tempdir()
        
#        filesToCopy <- list.files("www")
#        sapply(filesToCopy, function(x)
#              file.copy(from = file.path("www", x), 
#                  to = file.path(tmpDir, x),
#                  overwrite = TRUE)
#        )
        
        file.copy(normalizePath('www/reportBmd.Rmd'), 
            file.path(tmpDir, 'reportBmd.Rmd'), overwrite = TRUE)
        file.copy(normalizePath('www/flowChart.png'), 
            file.path(tmpDir, 'flowChart.png'), overwrite = TRUE)
        file.copy(normalizePath('www/Docx-landscape'), 
            file.path(tmpDir, 'Docx-landscape'), overwrite = TRUE)
        
        out <- render(file.path(tmpDir, 'reportBmd.Rmd'))
        file.rename(out, file)
        
      }
  )
  
  
  
}