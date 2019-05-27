#' Fit single model for the continuous or quantal response
#' @param data list, as returned from the f.scan() function;
#' @param shinyInput list, all parameter values as defined in the shiny UI
#' @param selectedModel, integer represents the model to be fit; run function
#' getModelNames to obtain all models and their coding
#' @param noFit boolean, if TRUE the model is not fitted, but default parameter
#' names and constraints can be obtained; default value is FALSE
#' @param estimateCed boolean, if TRUE the ced (i.e. bmd) is estimated along with
#' confidence interval using profile likelihood; default value is TRUE 
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list with for each model a list with all results obtained during 
#' the analysis with f.proast(); attribute "modelNames" contains the full names
#' of all models that were fitted; if error occurs return error message 
#' @return a list with all results obtained during the analysis with f.proast();
#' attribute "modelName" contains the full name of fitted model; 
#' if error occurs return error message 
#' @export
fitSingleModel <- function(data, shinyInput, selectedModel, noFit = FALSE, 
    estimateCed = TRUE, track = FALSE) {
  
  tryCatch({
        
        currentShinyInput <- shinyInput
        currentShinyInput$model.ans <- selectedModel
        currentShinyInput$main.ans <- 4
        
        modelChoices <- getModelNames(dtype = shinyInput$dtype)
        
        if (noFit) {
          
          currentShinyInput$parameterConstraints <- NULL
          
          fittedModel <- f.proast(odt = data, shinyInput = currentShinyInput, 
              track = track)
          
        } else {
          
          #Replace missing constraints by infinite
          if(!is.null(currentShinyInput$parameterConstraints)){
            
            currentShinyInput$parameterConstraints[is.na(currentShinyInput$parameterConstraints[,1]), 1] <- -Inf
            currentShinyInput$parameterConstraints[is.na(currentShinyInput$parameterConstraints[,2]), 2] <- Inf
            
          }
          
          if (selectedModel %in% modelChoices[1:2]) {
            
            fittedModel <- f.proast(odt = data, shinyInput = currentShinyInput, 
                track = track)
            
            # no BMD for null and full models
            fittedModel$CED <- NA
            fittedModel$conf.int <- rep(NA, 2)
            
          } else {
            
            if (estimateCed) {
              
              if(shinyInput$dtype == 5)
                currentShinyInput$main.ans <- c(4, 6, 7) else 
                currentShinyInput$main.ans <- c(4, 6)
              
            }
            
            fittedModel <- f.proast(odt = data, shinyInput = currentShinyInput, 
                track = track)
            
            if (is.null(fittedModel$CED)) {
              
              # for other models, labelled 'CED.matr' -> necessary?
              fittedModel$CED <- fittedModel$CED.matr 
              
            }
            
          }  
          
        }
        
        selectedName <- names(modelChoices)[modelChoices == selectedModel]
        
        attr(fittedModel, "modelName") <- selectedName
        
        
        return(fittedModel)
        
      }, error = function(err) {
        
        modelChoices <- getModelNames(dtype = shinyInput$dtype)
        selectedName <- names(modelChoices)[modelChoices == selectedModel]
        
        err$message <- paste(selectedName, "model:", err$message)
        
        return(err)
        
      })   
  
}




#' Fit all covariate combinations for model parameters in given model family
#' @param data list, as returned from the f.scan() function;
#' @param shinyInput list, all parameter values as defined in the shiny UI
#' @param selectedModel, integer represents the model to be fit; run function
#' getModelNames to obtain all models and their coding
#' @param noFit boolean, if TRUE the model is not fitted, but default parameter
#' names and constraints can be obtained; default value is FALSE
#' @param estimateCed boolean, if TRUE the ced (i.e. bmd) is estimated along with
#' confidence interval using profile likelihood; default value is TRUE 
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list with savedResults, for each model a list with all results 
#' obtained during the analysis with fitSingleModel() 
#' @export
fitAllCovariates <- function(data, shinyInput, selectedModel, 
    noFit = FALSE, estimateCed = TRUE, track = FALSE) {
  
  
  modelParameters <- getModelParameters(selectedModel = selectedModel,
      dtype = shinyInput$dtype)
  
  allCombinations <- expand.grid(fct1.no = unique(c(0, shinyInput$fct1.no)), 
      fct2.no = unique(c(0, shinyInput$fct2.no)), 
      fct3.no = unique(c(0, shinyInput$fct3.no)),
      fct4.no = unique(c(0, shinyInput$fct4.no)), 
      fct5.no = unique(c(0, shinyInput$fct5.no)))
  
  
  if (all(is.na(modelParameters))) {
    
    selectedCombinations <- allCombinations[1, , drop = FALSE]
    
  } else {
    
    selectedCombinations <- unique(allCombinations[, modelParameters, drop = FALSE])
    
  }
  
  savedResults <- lapply(seq_len(nrow(selectedCombinations)), function(i){
        
        shinyInput[names(selectedCombinations)] <- selectedCombinations[i, ]
        
        fitSingleModel(data = data, shinyInput = shinyInput,
            selectedModel = selectedModel, noFit = noFit, 
            estimateCed = estimateCed, track = track)
        
      })
  
  
  return(savedResults)
  
}




#' Fit all model families for the continuous or quantal response
#' @param data list, as returned from the f.scan() function;
#' @param shinyInput list, all parameter values as defined in the shiny UI
#' @param fitCovariateCombinations boolean, if TRUE all possible combinations of 
#' the included covariates for parameters are fitted; default value is FALSE
#' @param noFit boolean, if TRUE the model is not fitted, but default parameter
#' names and constraints can be obtained; default value is FALSE
#' @param estimateCed boolean, if TRUE the ced (i.e. bmd) is estimated along with
#' confidence interval using profile likelihood; default value is TRUE 
#' @param modelIndices numeric vector, indices of the models to be fit; if NULL
#' values are obtained from getModelNames(); default value is NULL
#' @param track logical, if TRUE (FALSE by default) print the name of the function which is currently being run
#' @return list with for each model a list with all results obtained during 
#' the analysis with fitSingleModel()
#' @export
fitAllModels <- function(data, shinyInput, fitCovariateCombinations = FALSE,
    noFit = FALSE, estimateCed = TRUE, 
    modelIndices = NULL, track = FALSE) {
  
  if (is.null(modelIndices))
    modelIndices <- getModelNames(dtype = shinyInput$dtype)
  
  if (!fitCovariateCombinations) {
    
    savedResults <- lapply(modelIndices, function(iModel){
          
          fitSingleModel(data = data, shinyInput = shinyInput,
              selectedModel = iModel, noFit = noFit, estimateCed = estimateCed, 
              track = track)
          
        })
    
  } else {
    
    savedResults <- lapply(modelIndices, function(iModel){
          
          fitAllCovariates(data = data, shinyInput = shinyInput,
              selectedModel = iModel, noFit = noFit, estimateCed = estimateCed,
              track = track)
          
        })
    
    savedResults <- do.call(c, savedResults)
    
  }
  
  
  return(savedResults)
  
}



#' Bind results of one/several fitted model(s) into one summary table
#' @param savedResults list as returned by the function fitAllModels() or 
#' fitAllCovariates() 
#' @return data frame containing for each model 
#' \itemize{
#'  \item{model: }{character, model name}
#'  \item{parameterNames: }{character, names of the parameters for which a covariate was included}
#'  \item{npar: }{integer, number of parameters in the model}
#'  \item{loglik: }{numeric, estimated log likelihood}
#'  \item{aic: }{numeric, estimated AIC}
#'  \item{bmd, bmdl, bmdu: }{per covariate value, the estimate, lower and upper
#'  bound of the benchmark dose}
#'  \item{converged: }{boolean, did fitting procedure ended with convergence}
#' }
#' @export
bindModelResults <- function(savedResults) {
  
  # A. Multiple models
  if (length(savedResults[[1]]) > 1) {
    
    # Exclude models with error from table
    withoutError <- !sapply(savedResults, function(x) is(x, "error"))
    usedResults <- savedResults[withoutError]
    
    isFitted <- !sapply(usedResults, function(x) attr(x, "modelName")) %in% c('Null', 'Full', 'Exp model 4')
    fittedResults <- usedResults[isFitted]
    
    parameterCombinations <- t(sapply(fittedResults, function(x) 
              c(x$fct1.no, x$fct2.no, x$fct3.no, x$fct4.no, x$fct5.no)))
    fittedParameterNames <- sapply(seq_len(nrow(parameterCombinations)), function(i) {
          
          includedParameters <- getModelParameters(selectedModel = NA,
              dtype = fittedResults[[1]]$dtype)
          parameterNames <- names(includedParameters)[parameterCombinations[i, includedParameters] != 0]
          if (length(parameterNames) == 0) parameterNames <- ""
          return(paste(parameterNames, collapse = ", "))
          
        })
    
    multiLevels <- sapply(usedResults, function(x) length(x$levelNames))
    
    if (any(multiLevels > 1)) {
      
      maxLevels <- which.max(multiLevels)
      levelNames <- lapply(usedResults, function(x) x$levelNames)[[maxLevels]]
      
      
      expandBmds <- function(variable = c("bmd", "bmdl", "bmdu"), fittedResults) {
        
        fittedValues <- t(sapply(fittedResults, function(x) {
                  
                  # Match levels and results with x$levelNames
                  tmpValues <- as.vector(switch(variable,
                      "bmd" = {
                        if(x$dtype == 3)
                          x$CED.matr[, x$CES.cat] * x$sf.x else
                          x$CED * x$sf.x
                      }, 
                      "bmdl" = x$conf.int[,1],
                      "bmdu" = x$conf.int[,2]))
                  
                  
                  if (length(tmpValues) < length(levelNames))
                    tmpValues <- rep(tmpValues, length.out = length(levelNames))
                  
                  if (length(tmpValues) > length(levelNames))
                    tmpValues <- unique(tmpValues)[seq_along(levelNames)]
                  
                  
                  return(tmpValues)
                  
                }))
        
        newValues <- matrix(data = NA, nrow = length(usedResults), 
            ncol = ncol(fittedValues))
        newValues[isFitted,] <- fittedValues
        
        return(newValues)
        
      }
      
      bmd <- expandBmds(variable = "bmd", fittedResults = fittedResults)
      bmdl <- expandBmds(variable = "bmdl", fittedResults = fittedResults)
      bmdu <- expandBmds(variable = "bmdu", fittedResults = fittedResults)
      
      bmdValues <- as.data.frame(matrix(rbind(bmd, bmdl, bmdu), ncol = 3*length(levelNames)))
      colnames(bmdValues) <- paste0(c("bmd.", "bmdl.", "bmdu."), rep(levelNames, each = 3))
      
      columnNames <- c("Model", "Included covariate(s)", 
          "Number of parameters", "Log-likelihood", "AIC", 
          paste(c("BMD -", "BMDL -", "BMDU -"), rep(levelNames, each = 3)), 
          "Converged")
      
    } else {
      
      bmd <- sapply(usedResults, function(x) {
            if(x$dtype == 3)
              x$CED.matr[1, x$CES.cat] * x$sf.x
            else x$CED[1] * x$sf.x
          }) 
      
      bmdl <- sapply(usedResults, function(x) x$conf.int[1])
      bmdu <- sapply(usedResults, function(x) x$conf.int[2])
      
      bmdValues <- data.frame(bmd = bmd, bmdl = bmdl, bmdu = bmdu)
      
      columnNames <- c("Model", "Included covariate(s)", 
          "Number of parameters", "Log-likelihood", "AIC", 
          "BMD", "BMDL", "BMDU", "Converged")
      
    }
    
    parameterNames <- rep("", length(usedResults))
    parameterNames[isFitted] <- fittedParameterNames
    
    summaryTable <- cbind(data.frame(
            model = sapply(usedResults, function(x) attr(x, "modelName")),
            parameterNames = parameterNames,
            npar = sapply(usedResults, function(x) x$npar),
            loglik = sapply(usedResults, function(x) x$loglik),
            aic = sapply(usedResults, function(x) 
                  2 * x$npar - 2 * x$loglik), 
            row.names = NULL, stringsAsFactors = FALSE),
        bmdValues)
    
    summaryTable$converged <- as.logical(sapply(usedResults, function(x) x$converged))
    
    attr(summaryTable, "columnNames") <- columnNames
    row.names(summaryTable) <- seq_along(savedResults)[withoutError]
    
    # B. Single model
  } else {
    
    levelNames <- savedResults$levelNames
    
    summaryTable <- with(savedResults, {
          
          multiLevels <- length(levelNames)
          
          if (attr(savedResults, "modelName") %in% c('Null', 'Full', 'Exp model 4')) {
            
            bmdValues <- data.frame(bmd = NA, bmdl = NA, bmdu = NA)
            fittedParameterNames <- ""
            
          } else {
            
            if (multiLevels > 1) {
              
              expandBmds <- function(variable = c("bmd", "bmdl", "bmdu"), savedResults){
                
                # Match levels and results with x$levelNames
                tmpValues <- switch(variable,
                    "bmd" = {
                      if(dtype == 3)
                        bmd <- CED.matr[, CES.cat] * sf.x else 
                        bmd <- CED * sf.x
                    }, 
                    "bmdl" = conf.int[,1],
                    "bmdu" = conf.int[,2])
                
                if (length(tmpValues) < length(levelNames))
                  tmpValues <- rep(tmpValues, length.out = length(levelNames))
                
                if (length(tmpValues) > length(levelNames))
                  tmpValues <- unique(tmpValues)
                
                return(as.numeric(tmpValues))
                
              }
              
              bmd <- expandBmds(variable = "bmd", savedResults = savedResults)
              bmdl <- expandBmds(variable = "bmdl", savedResults = savedResults)
              bmdu <- expandBmds(variable = "bmdu", savedResults = savedResults)
              
              bmdValues <- as.data.frame(matrix(rbind(bmd, bmdl, bmdu), ncol = 3*length(levelNames)))
              colnames(bmdValues) <- paste0(c("bmd.", "bmdl.", "bmdu."), rep(levelNames, each = 3))
              
            } else {
              
              if (dtype == 3)
                bmd <- CED.matr[1, CES.cat] * sf.x else 
                bmd <- CED[1] * sf.x
              
              bmdl <- conf.int[1]
              bmdu <- conf.int[2]
              
              bmdValues <- data.frame(bmd = bmd, bmdl = bmdl, bmdu = bmdu)
              
            }
            
            parameterCombinations <- c(fct1.no, fct2.no, fct3.no, fct4.no, fct5.no)
            parameterNames <- names(getModelParameters(selectedModel = NA,
                    dtype = dtype))[parameterCombinations != 0]
            if (length(parameterNames) == 0) parameterNames <- ""
            fittedParameterNames <- paste(parameterNames, collapse = ", ")
            
          }
          
          summaryTable <- cbind(data.frame(
                  model = attr(savedResults, "modelName"),
                  parameterNames = fittedParameterNames,
                  npar = npar,
                  loglik = loglik,
                  aic = 2 * npar - 2 * loglik, 
                  row.names = NULL, stringsAsFactors = FALSE),     
              bmdValues)
          summaryTable$converged <- as.logical(converged)
          
          if (multiLevels > 1 & !(summaryTable$model %in% c('Null', 'Full', 'Exp model 4'))) {
            
            attr(summaryTable, "columnNames") <- c("Model", "Included covariate(s)", 
                "Number of parameters", "Log-likelihood", "AIC", 
                paste(c("BMD -", "BMDL -", "BMDU -"), rep(levelNames, each = 3)),
                "Converged")
            
          } else {
            
            attr(summaryTable, "columnNames") <- c("Model", "Included covariate(s)", 
                "Number of parameters", "Log-likelihood", "AIC", 
                "BMD", "BMDL", "BMDU", "Converged")
            
          }
          
          
          return(summaryTable)
          
        })
    
  }
  
  
  return(summaryTable)
  
}



#' Retain model with smallest AIC per family when considering covariates 
#' @param summaryTable data frame as returned by the function bindModelResults()
#' @return data frame, subset of summaryTable with one model per family, which 
#' has the smallest AIC; attribute columnNames same as original summaryTable 
#' @importFrom plyr ddply
#' @export
filterBestCovariates <- function(summaryTable) {
  
  # Return original table if not various covariate combinations
  bestCovariates <- ddply(summaryTable, "model", function(subFamily) {
        
        if (nrow(subFamily) == 1) return(subFamily)
        
        subFamily[which.min(subFamily$aic), ]
        
      })
  
  rownames(bestCovariates) <- bestCovariates$model
  bestCovariates <- bestCovariates[unique(summaryTable$model),]
  rownames(bestCovariates) <- NULL
  attr(bestCovariates, "columnNames") <- attr(summaryTable, "columnNames") 
  
  
  return(bestCovariates)
  
} 



#' Find the best model given a table of fitted models results
#' @param summaryTable data frame as returned by the function bindModelResults()
#' @return if one model in summaryTable NULL, else list with 
#' \itemize{
#' 	\item{aicWarning: }{character, warning in procedure to find best model}
#'  \item{acceptedAic: }{boolean, is fitted model accepted based on AIC}
#'  \item{minBmdl: }{numeric vector, minimum of bmdl values (among accepted models)}
#'  \item{maxBmdu: }{numeric vector, maximum of bmdu values (among accepted models)}
#'  \item{minAic: }{numeric, minimum of aic values (among accepted models)}
#'  \item{bestModel: }{character vector, name and included parameters of best (with minAic) model}
#'  \item{bestBmd: }{numeric vector, bmd values of best model}
#' }
#' @export
findBestModel <- function(summaryTable) {
  
  # A. Single model
  if (nrow(summaryTable) == 1) {
    
    return(NULL)
    
    # B. Multiple models    
  } else {
    
    isFitted <- !(summaryTable$model %in% c('Null', 'Full', 'Exp model 4'))
    aicFitted <- summaryTable$aic[isFitted]
    
    minAic <- min(aicFitted, na.rm = TRUE)
    accepted <- rep(NA, nrow(summaryTable))
    aicWarning <- NULL
    
    
    # Consider smallest aic for comparisons with Null and Full model
    if (any(summaryTable$model == "Null"))
      aicNull <- min(summaryTable$aic[(summaryTable$model == "Null")]) else
      aicNull <- Inf
    if (any(summaryTable$model %in% c("Full", 'Exp model 4')))
      aicFull <- min(summaryTable$aic[(summaryTable$model %in% c("Full", 'Exp model 4'))]) else
      aicFull <- Inf
    
    
    if (minAic > (aicNull - 2) & length(aicNull) > 0) {
      
      aicWarning <- "None of the fitted models is better than the null model: All fitted models' AIC values are larger than null model's AIC - 2"
      accepted[isFitted] <- FALSE
      
    } else {
      
      if (!any(summaryTable$converged)) {
        
        aicWarning <- "None of the fitted models converged"
        accepted[isFitted] <- FALSE
        
      } else if (minAic > (aicFull + 2) & length(aicFull) > 0) {
        
        aicWarning <- "None of the fitted models is at least as good as the full model: All fitted models' AIC values are larger than full model's AIC + 2"
        accepted[isFitted] <- FALSE
        
      } else accepted[isFitted] <- (aicFitted <= (minAic + 2))
      
    }
    
    if (!is.null(aicWarning))
      warning(aicWarning)
    
    isMultiLevel <- (ncol(summaryTable) > 9)
    
    # Find minBmdl and maxBmdu
    bmdl <- summaryTable[accepted & !is.na(accepted), 
        grepl("bmdl", colnames(summaryTable))]
    bmdu <- summaryTable[accepted & !is.na(accepted), 
        grepl("bmdu", colnames(summaryTable))]
    
    findBmd <- function(bmdValues, FUN = c(min, max)) {
      
      if (all(is.na(bmdValues))) {
        
        toReturn <- NA
        
      } else {
        
        if (isMultiLevel) {
          
          if (sum(accepted, na.rm = TRUE) == 1) 
            toReturn <- bmdValues else 
            toReturn <- apply(bmdValues, 2, FUN, na.rm = TRUE)
          
        } else toReturn <- do.call(FUN, list(bmdValues, na.rm = TRUE))
        
      }  
      
      return(toReturn)
    }
    
    minBmdl <- findBmd(bmdValues = bmdl, FUN = min)
    maxBmdu <- findBmd(bmdValues = bmdu, FUN = max)
    
    bestModel <- summaryTable[which(summaryTable$aic == minAic)[1], 
        c("model", "parameterNames")]
    
    bestBmdColumns <- grepl("bmd", colnames(summaryTable)) & 
        !grepl("bmdl", colnames(summaryTable)) & !grepl("bmdu", colnames(summaryTable))
    bestBmd <- summaryTable[which(summaryTable$aic == minAic), bestBmdColumns]
    
    
    extraInfo <- list(aicWarning = aicWarning, acceptedAic = accepted,
        minBmdl = minBmdl, maxBmdu = maxBmdu, 
        minAic = minAic, bestModel = bestModel, bestBmd = bestBmd)
    
    
    return(extraInfo)
    
  }
  
}


#' Summarize the results of all fitted models for continuous or quantal response
#' @param savedResults list as returned by the function fitAllModels() or 
#' fitAllCovariates() 
#' @return list with summaryTable and extraInfo, respectively as returned by 
#' bindModelResults() and findBestModel()
#' @export
summaryModels <- function(savedResults) {
  
  summaryTable <- bindModelResults(savedResults = savedResults)
  extraInfo <- findBestModel(summaryTable = summaryTable)
  
  
  summaryTable$converged <- ifelse(summaryTable$converged, "yes", "no")
  
  if (!is.null(extraInfo)) {
    
    summaryTable$accepted <- ifelse(is.na(extraInfo$acceptedAic), "",
        ifelse(extraInfo$acceptedAic, "yes", "no"))
    
    attr(summaryTable, "columnNames") <- c(attr(summaryTable, "columnNames"), "Accepted AIC")
    
    extraInfo$bestModelIndex <- which(extraInfo$minAic == 
            sapply(savedResults, function(x) 2 * x$npar - 2 * x$loglik))[1]
    extraInfo$responseName <- savedResults[[1]]$y.leg
    
  } else {
    
    extraInfo <- list(responseName = savedResults$y.leg)
    
  }
  
  
  return(list(summaryTable = summaryTable, extraInfo = extraInfo))
  
}




#' Forest plot for estimated bmd, bmdl and bmdu values per response and covariate 
#' @param bmdData data.frame with response name, estimated bmd, bmdl and bmdu 
#' values
#' @param groupNames character vector, names of the covariate groups;
#' default value is NULL 
#' @return no object is returned, show plot in current device
#' @importFrom RColorBrewer brewer.pal
#' @importFrom forestplot forestplot fpColors fpTxtGp
#' @importFrom grid gpar
#' @export
makeForestPlot <- function(bmdData, groupNames = NULL) {
  
  if (is.null(bmdData))
    return(NULL)
  
  
  if (!is.null(groupNames)) {
    
    nGroups <- length(groupNames)
    colors <- brewer.pal(max(3, nGroups), "Dark2")[1:nGroups]
    
    bmdColumns <- grepl("bmd.", names(bmdData), fixed = TRUE)
    bmdlColumns <- grepl("bmdl.", names(bmdData), fixed = TRUE)
    bmduColumns <- grepl("bmdu.", names(bmdData), fixed = TRUE)
    
  } else {
    
    nGroups <- 1
    colors <- "black"
    
    bmdColumns <- which(names(bmdData) == "bmd")
    bmdlColumns <- which(names(bmdData) == "bmdl")
    bmduColumns <- which(names(bmdData) == "bmdu")
    
  }
  
  if (!any(bmdColumns)) {
    
    bmdColumns <- bmdlColumns
    boxsize <- 0
    
    tableText <- cbind(
        c("Response", as.character(bmdData$response)),
        c("BMDL", apply(signif(bmdData[, bmdlColumns, drop = FALSE], 5), 1, 
                function(x) paste(x, collapse = "\n"))),
        c("BMDU", apply(signif(bmdData[, bmduColumns, drop = FALSE], 5), 1, 
                function(x) paste(x, collapse = "\n"))))
    
  } else {
    
    boxsize <- 0.1
    
    tableText <- cbind(
        c("Response", as.character(bmdData$response)),
        c("BMDL", apply(signif(bmdData[, bmdlColumns, drop = FALSE], 5), 1, 
                function(x) paste(x, collapse = "\n"))),
        c("BMD", apply(signif(bmdData[, bmdColumns, drop = FALSE], 5), 1,
                function(x) paste(x, collapse = "\n"))), 
        c("BMDU", apply(signif(bmdData[, bmduColumns, drop = FALSE], 5), 1, 
                function(x) paste(x, collapse = "\n"))))
    
  }
  
  dataRange <- signif(range.default(bmdData[, -1], finite = TRUE)*c(0.9, 1.1), 3)
  
  
  # Example forestplot from https://www.r-bloggers.com/forest-plot-with-horizontal-bands/
  plot.new()
  forestplot(tableText,
#      xticks = seq(0, roundUp(max(bmdData[, bmduColumns], na.rm = TRUE)*1.1), length.out = 4),
      boxsize = boxsize,
      clip = dataRange,
      mean = rbind(rep(NA, nGroups), bmdData[, bmdColumns, drop = FALSE]),
      lower = rbind(rep(NA, nGroups), bmdData[, bmdlColumns, drop = FALSE]),
      upper = rbind(rep(NA, nGroups), bmdData[, bmduColumns, drop = FALSE]),
      lwd.ci = 2, ci.vertices = TRUE, 
      col = fpColors(box = colors, lines = colors),
      txt_gp = fpTxtGp(ticks = gpar(cex = 1.2), cex = 1.5))
  if (!is.null(groupNames)) {
    par(xpd = TRUE)
    legend("topright", groupNames, col = colors, bty = "n", pch = 15, cex = 1.5)
  }  
  
  #   toPlot <- plot_ly(data = bmdData, x = ~bmd, y = ~response,
  ##      if (!is.null(bmdData$covariate))
  ##        color = ~covariate, 
#      color = ~covariate,
#      type = "scatter", mode = "markers",
#      hoverinfo = "text", 
#      text = ~paste0(bmdText["bmd"], ": ", bmd,
#          "</br>", bmdText["bmdl"], ": ", bmdl,
#          "</br>", bmdText["bmdu"], ": ", bmdu),
#      error_x = list(type = "data", symmetric = FALSE, array = ~bmdu,
#          arrayminus = ~bmdl))
#  
#  toPlot %>% 
#      layout(xaxis = list(title = "BMD values"), 
#          yaxis = list(title = ""))  
  
}


#' Index in Proast and full names of the fitted models
#' @param dtype integer, indicating the response type
#' @return named vector with model indices for specific data type and names the 
#' model names
#' @export
getModelNames <- function(dtype){
  
  if (dtype %in% c(1, 5, 10, 15, 25, 250, 26, 260)) {
    
    modelNames <- c("Null" = 1, "Full" = 11, "Exp model 3" = 13,
        "Exp model 5" = 15, "Hill model 3" = 23, "Hill model 5" = 25)
    
  } else if (dtype == 3){
    
    modelNames <- c("Null" = 1, "Exp model 4" = 14, "Exp model 3" = 13,
        "Exp model 5" = 15, "Hill model 3" = 23, "Hill model 5" = 25)
    
  } else {
    
    modelNames <- c("Null" = 1, "Full" = 14, "Logistic" = 26,
        "Probit" = 25, "Log-logistic" = 18, "Log-probit" = 21,
        "Weibull" = 19, "Gamma" = 24, "Two-stage" = 16)
    
  }
  
  
  return(modelNames)
  
}



#' Obtain the available model parameters for a given model family
#' @param selectedModel, integer represents the model to be fit; run function
#' getModelNames to obtain all models and their coding; default value is NA
#' @param dtype integer, indicating the response type
#' @return vector with names the model parameters
#' @export
getModelParameters <- function(selectedModel = NA, dtype) {
  
  if (dtype %in% c(1, 5, 10, 15, 25, 250, 26, 260, 3)) 
    codes <- c("a" = 1, "b" = 2, "var" = 3, "c" = 4, "d" = 5) else 
    codes <- c("a" = 1, "b" = 2, "c" = 4)
  
  
  if (is.na(selectedModel)) {
    
    return(codes)
    
  } else {
    
    if (dtype %in% c(1, 5, 10, 15, 25, 250, 26, 260)) {
      
      if (selectedModel %in% c(1, 11)) 
        modelParameters <- NA else if (selectedModel %in% c(13, 23))
        modelParameters <- c("a", "b", "var", "d") else
        modelParameters <- c("a", "b", "var", "c", "d")
      
    } else if (dtype == 3) {
      
      if (selectedModel %in% c(1, 11)) 
        modelParameters <- NA else if (selectedModel %in% c(13, 23))
        modelParameters <- c("a", "b", "var", "d") else if(selectedModel == 14)
        modelParameters <- c("a", "b", "var", "c") else
        modelParameters <- c("a", "b", "var", "c", "d")
      
    } else {
      
      if (selectedModel %in% c(1, 14)) 
        modelParameters <- NA else if (selectedModel %in% c(25, 26, 18))
        # TODO due to error in proast we cannot include factor for parameter c
        modelParameters <- c("a", "b") else
        modelParameters <- c("a", "b", "c")
      
    }
    
    return(codes[modelParameters])
    
  }
  
}



#' Find which fitted models' results match the listed model and parameter names  
#' @param savedResults list as returned by the function fitAllModels() or 
#' fitAllCovariates() 
#' @param modelNames character vector, names of the models to be selected; 
#' see getModelNames() for available options 
#' @param parameterNames character vector, same length as \code{modelNames},
#' names of the model parameters to be selected; see getModelParameters() for
#' available options; default value is NULL
#' @return numeric vector, indices of the matching fitted results in the 
#' list \code{savedResults}
#' @export
matchWithResults <- function(savedResults, modelNames, parameterNames = NULL) {
  
  if (length(savedResults[[1]]) > 1) {    
    # Exclude models with error from table
    withoutError <- !sapply(savedResults, function(x) is(x, "error"))
    usedResults <- savedResults[withoutError]
    
    isFitted <- !sapply(usedResults, function(x) attr(x, "modelName")) %in% c('Null', 'Full', 'Exp model 4')
    fittedResults <- usedResults[isFitted]
    
    allModelNames <- rep("", length(savedResults)) 
    allModelNames[withoutError] <- sapply(usedResults, function(x) attr(x, "modelName"))
    
    
    if (!is.null(parameterNames)) {
      
      if (length(modelNames) != length(parameterNames))
        stop("The lengths of modelNames and parameterNames need to be equal.")
      
      parameterCombinations <- t(sapply(fittedResults, function(x) 
                c(x$fct1.no, x$fct2.no, x$fct3.no, x$fct4.no, x$fct5.no)))
      fittedParameterNames <- sapply(seq_len(nrow(parameterCombinations)), function(i) {
            
            includedParameters <- getModelParameters(selectedModel = NA,
                dtype = fittedResults[[1]]$dtype)
            parameterNames <- names(includedParameters)[parameterCombinations[i, includedParameters] != 0]
            if (length(parameterNames) == 0) parameterNames <- ""
            return(paste(parameterNames, collapse = ", "))
            
          })
      
      subParameterNames <- rep("", length(usedResults))
      subParameterNames[isFitted] <- fittedParameterNames
      allParameterNames <- rep("", length(savedResults))
      allParameterNames[withoutError] <- subParameterNames  
      
      bestIndices <- which(paste(allModelNames, allParameterNames) %in% 
              paste(modelNames, parameterNames))
      
    } else {
      
      bestIndices <- which(allModelNames %in% modelNames)      
      
    }
    
  } else {
    
    bestIndices <- 1  
    
  }
  
  
  return(bestIndices)
  
}



#' Print for debugging
#' @param x R object that will be printed
#' @return NULL, print output in the console
#' @export
printer <- function(x){
  
  cat("MV", deparse(substitute(x)), "\n")
  print(x)
  
}

#' Paste elements of vector into string vector (for testing)
#' @param x vector
#' @return string of the form "c(<elements of x>)"
#' @export
pasteToVector <- function(x) {
  
  if (is.character(x))
    elements <- paste(paste0("'", x, "'"), collapse = ", ")
  else elements <- paste(x, collapse = ", ")
  
  paste0("c(", elements, ")")
  
}


