library(shinysky)
library(shinyjs)

jsScroll <- "shinyjs.scrollTo = function(params){$('html, body').animate({scrollTop: $(params['id']).offset().top}, 2000);}"

fluidPage(
    
    useShinyjs(),
    extendShinyjs(text = jsScroll),
    
    tags$br(),
    actionButton(inputId = "debug_console", label = "Connect with console"),
    tags$br(),
    tags$br(),
    
    downloadButton("downloadReport", label = "Download report"),
    
#    fluidRow(
#        column(12, downloadButton("downloadReport", label = "Download report"), 
#        downloadButton("downloadSink", label = "Download log file"))
#    ),
    
    helpText(
        h5(a(href="https://support-efsa.openanalytics.eu/projects/bmd/wiki", 
                target="_blank", "About"), align = "right"),
        h5(a(href="https://support-efsa.openanalytics.eu/projects/bmd/issues/new/", 
                target="_blank", "Report new issue"), align = "right")
    ),
    
    titlePanel(title = div(img(src = "EFSA_logo.JPG", 
                float = "top", height = "60px", hspace = "50px"),
            "Benchmark Dose Modelling"), 
        windowTitle = "Benchmark Dose Modelling"),
    
    tags$br(),
    
    tabsetPanel(
        
        tabPanel("Data",
            
            tags$br(),
            
            fluidRow(
                
                column(6, 
                    
                    wellPanel(
                        
                        fileInput('dataLoaded', '',
                            accept = c('text/csv', 
                                'text/comma-separated-values,text/plain', 
                                '.csv')),
                        
                        actionLink(inputId = "helpProast", label = "Data format",
                            icon = icon("info-circle")),
                        
                        conditionalPanel("input.helpProast % 2 == 1",
                            em("Option 1: Proast data"),
                            p(strong("Line 1:"), "one-word title, used in plots (no minus signs or spaces allowed)",
                                br(), strong("Line 2:"), "the number of columns of the data matrix",
                                br(), strong("Line 3:"), "a code for the data type in that particular column, choose between",
                                tags$ul(
                                    tags$li("0: Non-response"),
                                    tags$li("1: Continuous response"),
                                    tags$li("2: Binary response"),
                                    tags$li("3: Ordinal response"),
                                    tags$li("4: Quantal response"),
                                    tags$li("5: Nested continuous response"),
                                    tags$li("6: Nested quantal response"),
                                    tags$li("10: Mean (continuous) response")
                                ),
                                strong("Line 4:"), "one-word title for each column (no minus signs or spaces allowed)",
                                br(), strong("Line 5 onwards:"), "the data matrix"),
                            em("Option 2: Raw data"),
                            p(strong("Line 1:"), "one-word title for each column",
                                br(), strong("Line 2 onwards:"), "the data matrix")
                        ),
                        
                        conditionalPanel("input.helpProast % 2 == 0",
                            tags$br()
                        ),
                        
                        
                        fluidRow(
                            column(6, 
                                radioButtons('sep', 'Data separator',
                                    c("Tab (.txt file)" = '\t',
                                        "Space (.txt file)" = ' ',
                                        "Comma (.csv file)" = ',',
                                        "Semicolon (.csv file)"= ';')
                                )
                            ),
                            column(6, 
                                radioButtons('dec', 'Decimal separator',
                                    c("Point" = '.', "Comma" = ','))
                            )
                        ),
                        
                        fluidRow(
                            column(6, uiOutput("subsetVariable")),
                            column(6, uiOutput("levelSubset"))
                        ),
                        
                        uiOutput('selectedResponses'),
                        uiOutput('nonResponses'),
                        
                        selectInput("dtype", label = "Type of response",
                            choices = c("continuous" = 1,  
                                "quantal" = 4, "binary" = 2, "ordinal" = 3)),
                        
                        conditionalPanel("input.dtype == '1'",
                            selectInput("dtype2", "", 
                                choices = c("individual data" = "individual",
                                    "summary data" = "summary"))                            
                        ),
                        
                        conditionalPanel("input.dtype == '1' | input.dtype == '4'",
                            checkboxInput("isLitter", "Litter effect")
                        ),
                        
                        
                        conditionalPanel("input.dtype != '2' & input.dtype != '4'", 
                            
                            em("Model averaging is only available for quantal response data.")
                        
                        )
                    )
                
                ),
                
                column(6,
                    
                    DT::dataTableOutput('dataLoaded')
                
                )
            
            )
        
        ),
        
        tabPanel("Fit Models",
            
            tags$br(),
            
            fluidRow(
                
                column(6, 
                    
                    wellPanel(
                        
                        h4("Define variables"),
                        
                        # Main questions
                        uiOutput("parameterQuestions"),
                        uiOutput("parameterQuestions2"),
                        uiOutput("blankLine"),
                        
                        h4("Models to fit"),
                        
                        fluidRow(
                            column(6,
                                selectInput("singleModel", "",
                                    choices = c("Set of models" = "no", "Single model" = "yes"))
                            ),
                            column(6,
                                uiOutput("modelChoices")
                            )
                        )
                    
                    ),
                    
                    uiOutput("selectedModels")
                
                ),
                
                column(6, 
                    
                    wellPanel(
                        
                        tabsetPanel(
                            
                            tabPanel(title = "Additional settings",
                                
                                uiOutput("ces.ans"),
                                uiOutput("ces"),
                                
                                numericInput("conf.lev", "Confidence level for the BMD confidence intervals", 
                                    value = 0.9),
                                
                                conditionalPanel("input.dtype == '2' | (input.dtype == '4' & input.isLitter == false) | input.dtype == '5'",
                                    numericInput("nBootstraps", "The number of bootstrap runs for calculating BMD confidence intervals",
                                        value = 200)
                                )
                            
                            ), 
                            
                            tabPanel(title = "Model parameters",
                                
                                uiOutput("factors"),
                                conditionalPanel("input.singleModel == 'no'", 
                                    checkboxInput("allCovariates", "Show results for all covariate combinations",
                                        value = FALSE)
                                ),
                                actionLink("advancedSettings", "Advanced settings", icon = icon("cog")),
                                
                                conditionalPanel("input.advancedSettings % 2 == 1",
                                    
                                    uiOutput("includeCovariates"),
                                    uiOutput("extraCovariates"),
                                    
                                    uiOutput("setParameterValues")
                                )
                            
                            )
                        )
                    )
                )
            ),
            
            h6("Output generated with Proast, version 61.3", align = "right"),
            tags$hr(),   
            
            uiOutput("fitModel"),
            tags$br(),
            uiOutput("summaryAnalysis")
            
        )    
    )    

)
