#' Test function for defining parameter constraints
#' @param shinyInput list with input from Shiny
#' @param data dataset
#' @param modelAns integer, type of response fitted
#' @param doPrint logical, if TRUE, additional information concerning the test is printed
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @import testthat
#' @return no returned value
#' @export
testConstraints <- function(shinyInput, data, modelAns, doPrint, track){
  
  dtype <- shinyInput$dtype
  shinyInput$ces.ans <- 1
  
  test_that("Adapt default parameter constraints", {
        
        result <- lapply(modelAns, function(model.ans){
              
              cat("** model.ans =", model.ans, "** \n")
              shinyInput$model.ans <- model.ans
              
              # Initial run with default parameter constraints
              shinyInput$parameterConstraints <- NULL
              shinyInput$main.ans <- NULL
              
              fittedModel <- f.proast(odt = data, shinyInput = shinyInput, track = track)
              
              if(doPrint)
                cat(" * Initial parameter constraints \n",
                    f.text.par(fittedModel), "\n", signif(fittedModel$lb, 5), "\n",
                    signif(fittedModel$ub, 5), "\n \n")
              
              # Adapt default parameter constraints
              shinyInput$parameterConstraints <- data.frame(
                  lowerBound = fittedModel$lb, 
                  upperBound = fittedModel$ub, 
                  stringsAsFactors = FALSE)
              
              if(shinyInput$cont){
                
                shinyInput$parameterConstraints$upperBound[1] <- Inf
                if(model.ans != 1)
                  shinyInput$parameterConstraints$lowerBound[2] <- -Inf
                
              } else {
                
                shinyInput$parameterConstraints$lowerBound[1] <- -Inf
                if(model.ans != 1)
                  shinyInput$parameterConstraints$upperBound[2] <- Inf
                
              }
              
              
              shinyInput$main.ans <- 4
              
              finalModel <- f.proast(odt = data, shinyInput = shinyInput, track = track)
              tmp <- f.plot.all(finalModel, track = track)
              
              if(doPrint)
                cat(" * New parameter constraints \n",
                    f.text.par(finalModel), "\n", signif(finalModel$lb, 5), "\n",
                    signif(finalModel$ub, 5), "\n",
                    "Loglik: ", signif(finalModel$loglik, 5), "\n", 
                    "MLE: ", signif(finalModel$MLE), "\n \n")
              
              if (model.ans != 11) { #fixed parameters (and thus no constraints) for full model
                
                expect_equal(finalModel$lb, shinyInput$parameterConstraints$lowerBound)
                expect_equal(finalModel$ub, shinyInput$parameterConstraints$upperBound)
                
              }
              
              return(finalModel)
              
            })
        
      })
  
}

#' Test function for defining starting values
#' @param shinyInput list with input from Shiny
#' @param data dataset
#' @param modelAns integer, type of response fitted
#' @param startValues list with starting values for each of the response fitted
#' This should be named with the modelAns as character; default value is NULL
#' @param results list with expected results, default value is NULL;
#' should contain the following items
#' \itemize{
#' 	\item{par.start: }{list with starting values for the parameters, of length \code{modelAns}}
#' 	\item{loglik: }{vector with values for the log likelihood, of length \code{modelAns}}
#' 	\item{CED: }{(only for continous model) vector with CED values, of length \code{modelAns}}
#' 	\item{CED.matr: }{(only for quantal model) vector with CED values, 
#' of length (\code{modelAns} - 2) (without null and full models)}
#' 	\item{conf.int: }{vector with confidence intervals, with 2 rows (lower/upper bounds)
#' 	and with number of columns, for:
#' 	\itemize{
#' 		\item{continuous model}{length \code{modelAns}}
#' 		\item{quantal model}{length (\code{modelAns} - 2) (without null and full models)}
#' 	}}
#' }
#' @param doPrint logical, if TRUE, additional information concerning the test is printed
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @import testthat
#' @return no returned value
#' @export
testStartingValues <- function(shinyInput, data, modelAns, 
    startValues = NULL, results = NULL, doPrint, track){
  
  dtype <- shinyInput$dtype
  shinyInput$ces.ans <- 1
  nameCED <- ifelse(shinyInput$cont, "CED", "CED.matr")
  
  
  test_that(paste("Adapt default start values of model parameters, dtype =", dtype), {
        
        maxIndexFullModel <- ifelse(shinyInput$cont, 13, 16)
        if(shinyInput$dtype == 3) maxIndexFullModel <- 13
        
        # Adapt default start values
        if(is.null(startValues)) {
          
          cat("** Generate initial start values **\n")
          
          initialResult <- lapply(modelAns, function(model.ans){
                
                cat("* model.ans =", model.ans, "\n")
                shinyInput$model.ans <- model.ans
                
                if (model.ans < maxIndexFullModel)
                  shinyInput$main.ans <- 4 else if(dtype == 5)
                  shinyInput$main.ans <- c(4, 6, 7) else
                  shinyInput$main.ans <- c(4, 6)
                
                if (shinyInput$dtype == 3 & model.ans == 14) 
                  shinyInput$main.ans <- 4
                
                
                # Initial run with default start values
                shinyInput$startValues <- NULL
                
                fittedModel <- f.proast(odt = data, shinyInput = shinyInput, track = track)
                
                if(doPrint)
                  cat(f.text.par(fittedModel), "\n", 
                      signif(fittedModel$par.start, 5), "\n \n")
                
                return(fittedModel)
                
              })
          
          startValues <- lapply(initialResult, function(x) x$par.start)
          names(startValues) <- as.character(modelAns)
          
          if(is.null(results))
            results <- list(
                par.start = lapply(initialResult, function(x) x$par.start),
                loglik = sapply(initialResult, function(x) x$loglik),
                CED.tmp = sapply(initialResult[-c(1,2)], function(x) x[[nameCED]]),
                conf.int =  sapply(initialResult[-c(1,2)], function(x) x$conf.int)
            )
          
          names(results)[3] <- nameCED          
          
        }
        
        
        cat("** Analysis with new start values **\n")
        
        result <- lapply(modelAns, function(model.ans){
              
              cat("* model.ans =", model.ans, "\n")
              shinyInput$model.ans <- model.ans
              
              if (model.ans < maxIndexFullModel)
                shinyInput$main.ans <- 4 else if(dtype == 5)
                shinyInput$main.ans <- c(4, 6, 7) else 
                shinyInput$main.ans <- c(4, 6)
              
              if (shinyInput$dtype == 3 & model.ans == 14) 
                shinyInput$main.ans <- 4
              
              shinyInput$startValues <- startValues[[as.character(model.ans)]]
              
              finalModel <- f.proast(odt = data, shinyInput = shinyInput, track = track)
              
              if(doPrint)
                cat(f.text.par(finalModel), "\n", signif(finalModel$par.start, 5), "\n",
                    "Loglik: ", signif(finalModel$loglik, 5), "\n", 
                    "MLE: ", signif(finalModel$MLE), "\n \n")
              
              # check plot
              tmp <- f.plot.all(finalModel)
              
              return(finalModel)
              
            })
        
        ## checks
        if(!is.null(results)){
          
          # specified starting values of parameters
          expect_equal(sapply(result, function(x) signif(x$par.start, 4)), 
              sapply(results$par.start, function(x) signif(x, 4)))
          
          expect_equal(sapply(result, function(x) signif(x$loglik, 4)), 
              signif(results$loglik, 4))
          
          # For CED do not consider null and full model
          resultUsed <- result[-c(1,2)]
          
          expect_equal(sapply(resultUsed, function(x) signif(x[[nameCED]], 4)), 
              signif(results[[nameCED]], 4))
          
          # This test is not ok for dtype = 6
          if (dtype != 6)
          expect_equal(sapply(resultUsed, function(x) signif(x$conf.int, 4)), 
              signif(results$conf.int, 4))
          
        } else {
          
          resultUsed <- result[-c(1,2)]
          cat("conf.int:", pasteToVector(sapply(resultUsed, function(x) signif(x$conf.int, 4))))
          
        }
        
      })   
  
}

#' Test function for including covariates
#' @param shinyInput list with input from Shiny
#' @param shinyInputCovariates list with input from Shiny for parameters
#' concerning the covariates
#' @param data dataset
#' @param modelAns integer, type of response fitted
#' @param results list with expected results, should contain the following items
#' \itemize{
#' 	\item{par.start: }{list with starting values for the parameters, of length \code{modelAns}}
#' 	\item{loglik: }{vector with values for the log likelihood, of length \code{modelAns}}
#' 	\item{CED: }{vector with CED values, of length \code{modelAns}}
#' 	\item{conf.int: }{vector with confidence intervals, with 2 rows (lower/upper bounds)
#' 	and length \code{modelAns} columns}
#' }
#' @param doPrint logical, if TRUE, additional information concerning the test is printed
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @import testthat
#' @return no returned value
#' @export
testCovariates <- function(shinyInput, shinyInputCovariates,
    data, modelAns, results = NULL, doPrint, track){
  
  dtype <- shinyInput$dtype
  
  shinyInputUsed <- c(shinyInput, shinyInputCovariates)
  
  maxIndexFullModel <- ifelse(shinyInput$cont, 13, 16)
  if(shinyInput$dtype == 3) maxIndexFullModel <- 13
  
  test_that("Run f.proast() - fit model and plot for all models & calculate CED", {
        
        result <- lapply(modelAns, function(model.ans){
              
              cat("** model.ans =", model.ans, "** \n")              
              shinyInputUsed$model.ans <- model.ans
              
              if (model.ans < maxIndexFullModel)
                shinyInputUsed$main.ans <- 4 else if(dtype == 5)
                shinyInputUsed$main.ans <- c(4, 6, 7) else
                shinyInputUsed$main.ans <- c(4, 6)
              
              if (shinyInput$dtype == 3 & model.ans == 14) 
                shinyInputUsed$main.ans <- 4
              
              if(shinyInput$cont){
                
                tmpResult <- f.proast(odt = data, shinyInput = shinyInputUsed, track = track)
                tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
                
                return(tmpResult)
                
              } else {
                
                allResults <- lapply(1:4, function(ces.ans){
                      #ces.ans only for non-continuous response types
                      
                      if(ces.ans < 4){
                        
                        shinyInputUsed$ces.ans <- ces.ans		
                        
                        if(doPrint) cat("ces.ans = ", ces.ans, "\n")
                        
                        if (ces.ans %in% c(2,3) & model.ans == 16 & shinyInputCovariates$fct4.no != 0) {
                          
                          if(doPrint) cat("Error in proast61.3 \n")   
                          
                        } else if (ces.ans > 1 & model.ans == 25 & dtype != 3){
                          #Covariates not implemented for probit model
                          
                          if(doPrint) cat("Error in proast61.3 \n")
                          
                        } else {                            
                          
                          tmpResult <- f.proast(odt = data, shinyInput = shinyInputUsed, track = track)
                          tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
                          
                          return(tmpResult) 
                          
                        }
                        
                      } else { # issue with ces.ans == 4
                        
                        if(doPrint) cat("ces.ans = ", ces.ans, ": Model not tested, bug in the code \n")
                        
                      }
                      
                    })
                
                
                return(allResults[[1]])
                
              }
              
            })
        
        ## checks
        
        # specified covariate name
        specifiedCovariateNames <- grep("fct[[:digit:]].no", names(shinyInputCovariates), value = TRUE)
        # has some values for covariates been specified in input parameters
        if(length(specifiedCovariateNames) > 0){
          # are there different than 0
          specifiedCovariate <- names(which(sapply(specifiedCovariateNames,	function(name) 
                        any(shinyInputCovariates[[name]] != 0))))
          if(length(specifiedCovariate) > 0){
            tmp <- sapply(specifiedCovariate, function(name){
                  resultsNameCovariate <- sapply(result, function(x) x[[sub("no", "name", name)]])
                  # filter null values (e.g. if model doesn't contain the specified parameter)
                  isNullCovariate <- sapply(resultsNameCovariate, is.null)
                  if(any(isNullCovariate))
                    cat(paste("NOTE: Covariate for parameter", 
#                            sub("fct([[:digit:]]).no", "\\1", name),
                            name,
                            "has not been taken into account (NULL) for model(s):",
                            paste(modelAns[which(isNullCovariate)], collapse = ", ")), "\n\n")
                  if(any(!isNullCovariate)){
                    resultsNameCovariate <- unlist(resultsNameCovariate[!isNullCovariate])
                    
                    if(any(grepl("interaction", resultsNameCovariate))){
                      
                      expect_equal(resultsNameCovariate,
                          rep(paste0("interaction_",
                                  paste0(data$varnames[shinyInputCovariates[[name]]], collapse = "*")), 
                              length(resultsNameCovariate)))
                      
                    } else {
                      
                      expect_equal(resultsNameCovariate,
                          rep(data$varnames[shinyInputCovariates[[name]]], length(resultsNameCovariate)))
                      
                    }
                    
                  }
                })
          }
        }
        
        # values for log-likelihood
        if(!is.null(results))
          expect_equal(sapply(result, function(x) signif(x$loglik, 5)), results$loglik) else
        if(doPrint) cat("loglik:", pasteToVector(sapply(result, function(x) signif(x$loglik, 5))), "\n")
        
        
        # for CED: don't consider null and full model for quantal model
        resultUsed <- result[-c(1,2)]
        
        # CED
        nameCED <- ifelse(shinyInput$cont, "CED.origScale", "CED.matr")
        if(dtype == 5) nameCED <- "CED"
        # use unique for quantal model (var and c (methyleug))
        if(!is.null(results))
          expect_equal(sapply(resultUsed, function(x) signif(unique(x[[nameCED]]), 5)), results[[nameCED]]) else 
        if(doPrint) cat("CED:", pasteToVector(unlist(sapply(resultUsed, function(x) signif(unique(x[[nameCED]]), 5)))), "\n")
        
        # confidence intervals
        if(!is.null(results))
          expect_equal(sapply(resultUsed, function(x) signif(x$conf.int, 5)), results$conf.int) else
        if(doPrint) cat("conf.int:", pasteToVector(unlist(sapply(resultUsed, function(x) signif(x$conf.int, 5)))), "\n")
        
      })
  
  
  ## check if summary function works in performAnalysis.R
  test_that("Using summary function for set of models", {
        
        fittedModels <- fitAllModels(data = data, shinyInput = shinyInputUsed, 
            track = track)
        summaryResult <- summaryModels(savedResults = fittedModels)
        tmpTable <- summaryResult$summaryTable
        if(doPrint)	print(tmpTable)
        
        tmpPlot <- f.plot.all(fittedModels[[summaryResult$extraInfo$bestModelIndex[1]]], track = track)
        
      })
  
  test_that("Using summary function for single model", {
        
        tmp <- sapply(modelAns, function(iModel){
              
              if(doPrint) print(iModel)
              
              fittedModel <- fitSingleModel(data = data, shinyInput = shinyInputUsed,
                  selectedModel = iModel, track = track)
              
              if(is(fittedModel, "error")){
                
                if(doPrint) print(paste("Error for model", iModel))
                
              } else {
                
                tmpTable <- summaryModels(savedResults = fittedModel)$summaryTable
                if(doPrint)	print(tmpTable)
                
                tmpPlot <- f.plot.all(fittedModel, track = track)
                
              }
              
            })
        
      })
  
}

