#' Calculate weights based on fitted models' aic
#' @param aicValues numeric vector, each model's estimated aic value 
#' @return numeric vector, estimated weight for each model.
#' Same order as aicValues but models with missing aicValues are excluded
#' @export
calculateWeights <- function(aicValues){
  
  aicValues <- aicValues[!is.na(aicValues)]
  newAic <- aicValues - min(aicValues)  # To prevent numerical problems
  weights <- exp(-0.5*newAic) / sum(exp(-0.5*newAic))
  
  
  return(weights)
  
}



#' Calculate the weighted average of estimated response values
#' @param weights numeric vector, estimated weights as returned by calculateWeights() 
#' @param dose numeric vector, dose values at which response values are estimated
#' @param modelResults list, with results for each model, same length as weights. 
#' For each model a list with at least model.ans, regr.par, CES and ces.ans; 
#' these are by default included in result from f.proast(). Eventually 
#' contains also fct1 and fct2 if factors are included for the model parameters.
#' @param groupIndex integer, index for covariate combination of interest;
#' default value is 1; if NA average response for all groups is returned
#' @return numeric (vector), averaged response value at given dose for group(s)
#' of interest
#' @export
averageResponse <- function(weights, dose, modelResults, groupIndex = 1) {
  
  responseValues <- lapply(modelResults, function(iModel) {
        
        if (is.null(iModel$fct1))
          fct1 <- 1 else
          fct1 <- iModel$fct1
        
        if (all(fct1 == 1)) fct1 <- 1
        
        if (is.null(iModel$fct2))
          fct2 <- 1 else
          fct2 <- iModel$fct2
        
        if (all(fct2 == 1)) fct2 <- 1
        
        # Only unique combinations should be remained
        allFactors <- cbind(fct1, fct2)
        uniqueFactors <- allFactors[!duplicated(allFactors), , drop = FALSE]
        
        fct1 <- uniqueFactors[, 1]
        fct2 <- uniqueFactors[, 2]
        
        # Define names
        if (all(fct1 == 1)) 
          namesFct1 <- " " else 
          namesFct1 <- paste0(names(iModel$data.0)[iModel$fct1.no], ".", iModel$fct1.txt[fct1]) 
        
        if (all(fct2 == 1))
          namesFct2 <- " " else
          namesFct2 <- paste0(names(iModel$data.0)[iModel$fct2.no], ".", iModel$fct2.txt[fct2]) 
        
        allResponses <- data.frame(matrix(
                f.expect.bin(model.ans = iModel$model.ans, 
                    x = rep(dose, times = nrow(uniqueFactors)), 
                    regr.par = iModel$regr.par, 
                    fct1 = rep(fct1, each = length(dose)), 
                    fct2 = rep(fct2, each = length(dose)),
                    fct1.full = iModel$fct1.full, fct2.full = iModel$fct2.full,
                    CES = iModel$CES, ces.ans = iModel$ces.ans, 
                    track = FALSE), ncol = nrow(uniqueFactors)))
        
        responseNames <- apply(cbind(namesFct1, namesFct2),
            1, function(x) paste(unique(x), collapse = " & "))
        responseNames <- gsub("  & | &  ", "", responseNames)
        
        if (all(responseNames %in% c("", " ")))
          responseNames <- NULL else
          names(allResponses) <- responseNames
        
#          allResponses <- unlist(allResponses)
#          names(allResponses) <- NULL
#          
#          return(allResponses)
#          
#        } else { 
#          
         
          
          return(allResponses)  
          
#        }
        
      })
  
  
  weightedResponses <- lapply(seq_along(weights), function(i)
        responseValues[[i]] * weights[i])
  responseColumns <- sapply(weightedResponses, ncol)
  
    # Multiple covariate groups
  if (any(responseColumns > 1)) {
    
    if (length(unique(responseColumns)) > 1) {
      
      maxColumns <- max(responseColumns)
      toRepeat <- which(responseColumns < maxColumns)
      
      for (i in toRepeat) {
      
        weightedResponses[[i]] <- do.call(cbind, rep(weightedResponses[[i]], 
            each = maxColumns/responseColumns[i]))
    
      }
      
      averagedResponse <- Reduce('+', weightedResponses)
      
    }
      
    
    if (!is.na(groupIndex)) {
      
      averagedResponse <- averagedResponse[, groupIndex, drop = FALSE]
      
    }    
    
    # Single covariate group
  } else {
    
    weightedResponses <- do.call(cbind, weightedResponses)

    if (length(dose) > 1)
      averagedResponse <- apply(weightedResponses, 1, sum, na.rm = TRUE) else
      averagedResponse <- sum(weightedResponses, na.rm = TRUE)
    
  }
  
  
  return(averagedResponse)
  
}



#' Target function for optimization if ces.ans = 1
#' @param weights numeric vector, estimated weights as returned by calculateWeights() 
#' @param dose numeric, value for the dose at which response values are estimated
#' @param modelResults list, with results for each model, same length as weights. 
#' For each model a list with at least model.ans, regr.par, CES and ces.ans; 
#' these are by default included in result from f.proast(). Eventually 
#' contains also fct1 and fct2 if factors are included for the model parameters.
#' @param groupIndex integer, index for covariate combination of interest;
#' default value is 1; if NA results for all groups are returned
#' @param baseResponse numeric, estimated response value for dose = 0;
#' redundant but retained for convenience with bmr2() and bmr3()
#' @param maxResponse numeric, response value for which bmr should be calculated;
#' default value is NULL 
#' @return numeric (vector), estimated response value for a given dose and 
#' group(s) of interest; same object as returned by averageResponse()
#' @export
bmr1 <- function(weights, dose, modelResults, groupIndex = 1, baseResponse, 
    maxResponse = NULL) {
  
  if (!is.null(maxResponse)) {  
    
    bmr <- maxResponse
    
  } else {
    
    bmr <- averageResponse(weights = weights, dose = dose, 
        modelResults = modelResults, groupIndex = groupIndex)
    
  }
  
  return(bmr)
  
}

#' Target function for optimization if ces.ans = 2
#' @param weights numeric vector, estimated weights as returned by calculateWeights() 
#' @param dose numeric, value for the dose at which response values are estimated
#' @param modelResults list, with results for each model, same length as weights. 
#' For each model a list with at least model.ans, regr.par, CES and ces.ans; 
#' these are by default included in result from f.proast(). Eventually 
#' contains also fct1 and fct2 if factors are included for the model parameters.
#' @param groupIndex integer, index for covariate combination of interest;
#' default value is 1; if missing average response for all groups is returned
#' @param baseResponse numeric, estimated response value for dose = 0
#' @param maxResponse numeric, response value for which bmr should be calculated;
#' default value is NULL 
#' @return numeric (vector), estimated added risk for group(s) of interest
#' @export
bmr2 <- function(weights, dose, modelResults, groupIndex = 1, baseResponse, 
    maxResponse = NULL){
  
  if (!is.null(maxResponse)) {
    
    bmr <- maxResponse - baseResponse
    
  } else {
    
    bmr <- averageResponse(weights = weights, dose = dose, 
        modelResults = modelResults, groupIndex = groupIndex) - baseResponse
    
  }
  
  return(bmr)
  
} 

#' Target function for optimization if ces.ans = 3 
#' @param weights numeric vector, estimated weights as returned by calculateWeights() 
#' @param dose numeric, value for the dose at which response values are estimated
#' @param modelResults list, with results for each model, same length as weights. 
#' For each model a list with at least model.ans, regr.par, CES and ces.ans; 
#' these are by default included in result from f.proast(). Eventually 
#' contains also fct1 and fct2 if factors are included for the model parameters.
#' @param groupIndex integer, index for covariate combination of interest;
#' default value is 1; if missing average response for all groups is returned
#' @param baseResponse numeric, estimated response value for dose = 0
#' @param maxResponse numeric, response value for which bmr should be calculated;
#' default value is NULL 
#' @return numeric (vector), estimated extra risk for group(s) of interest
#' @export
bmr3 <- function(weights, dose, modelResults, groupIndex = 1, baseResponse,
    maxResponse = NULL) {
  
  if (!is.null(maxResponse)) {
    
    bmr <- (maxResponse - baseResponse) / (1 - baseResponse)
    
  } else {
    
    bmr <- (averageResponse(weights = weights, dose = dose, 
              modelResults = modelResults, groupIndex = groupIndex) - 
          baseResponse) / (1 - baseResponse)
    
  }
  
  return(bmr)
  
} 

#' Estimate the model-averaged BMD using numeric optimization techniques 
#' @param weights numeric vector, estimated weights as returned by calculateWeights() 
#' @param modelResults list, with results for each model, same length as weights. 
#' For each model a list with at least npar, loglik, x, y, model.ans, regr.par, 
#' CES and ces.ans; these are by default included in result from f.proast(). 
#' Eventually contains also fct1 and fct2 if factors are included for the model parameters.
#' @param nIterations integer, number of iterations for numeric optimization;
#' default value is 400
#' @param naiveApproach boolean, TRUE if the model-averaged BMD is estimated as the 
#' weighted average of bmd values, FALSE if the model-averaged BMD is estimated
#' based on weighted average of response values; default value is FALSE
#' @return numeric, the model-averaged BMD; error is returned if the numeric 
#' procedure did not converge after \code{nIterations} iterations
#' @export
optimizeBmd <- function(weights, modelResults, nIterations = 400,
    naiveApproach = FALSE){
  
  
  withError <- sapply(modelResults, function(x) is(x, "error"))
  usedResults <- modelResults[!withError]
  
  if (length(usedResults) != length(weights))
    warning("The number of fitted models without error does not equal the number of weights")
  
  CES <- unique(sapply(usedResults, function(iModel) iModel$CES))
  if (length(CES) != 1)
    stop("All modelResults should be based on same value for CES, currently", 
        paste(CES, collapse = ", "))
  
  ces.ans <- unique(sapply(usedResults, function(iModel) iModel$ces.ans))
  if (length(ces.ans) != 1)
    stop("All modelResults should be based on same value for ces.ans, currently", 
        paste(ces.ans, collapse = ", "))
  
  # Function for bmr that should be optimized
  bmrFunction <- paste0("bmr", ces.ans)
  if (ces.ans == 1) CES <- 0.50  #target: Predicted response value equals 0.50  
  
  
  increase <- sign(lm(usedResults[[1]]$y ~ usedResults[[1]]$x)[[1]][2])
  
  if (naiveApproach) {
    
    multiLevels <- sapply(modelResults, function(x) length(x$levelNames))
    
    summaryTable <- bindModelResults(savedResults = modelResults)
    
    if (any(multiLevels > 1)) {
      
      bmdValues <- summaryTable[, grepl("bmd\\.", colnames(summaryTable))]
      names(bmdValues) <- gsub("bmd\\.", "", names(bmdValues))
      return(apply(bmdValues * weights[!withError], 2, sum, na.rm = TRUE))      
      
    } else {
      
      bmdValues <- summaryTable$bmd    
      return(sum(bmdValues * weights[!withError], na.rm = TRUE))      
      
    }
    
  } else {
    
    # Averaged response for dose = 0
    baseResponses <- averageResponse(weights = weights, dose = 0, 
        modelResults = usedResults, groupIndex = NA)
    
    bmdValues <- sapply(seq_along(baseResponses), function(iGroup) {
          
          # Bracketing
          tmpDose <- 0.5
          multiplier <- 1
          bounded <- FALSE
          
          for (i in 1:20) {
            
            multiplier <- multiplier*2
            tmpDose <- 0.5*multiplier
            
            bmr <- do.call(bmrFunction, list(weights = weights, dose = tmpDose, 
                    modelResults = usedResults, baseResponse = baseResponses[iGroup], 
                    groupIndex = iGroup)) * increase
            
            if (bmr >= CES){
              
              bounded <- TRUE
              break
              
            }       
            
          }
          
          if (!bounded)
            return(NA)
          
          
          # Bisection 
          top <- log(tmpDose)
          bottom <- log(2.2e-308)
          mid <- (top + bottom)/2
          
          iterations <- 1 
          
          #  First, make sure the BMD is above bottom (we already know it is below top)
          #  if BMD is below bottom, just return 0, as bottom is already as small as
          #  we can represent.
          bmr <- do.call(bmrFunction, list(weights = weights, dose = exp(bottom), 
                  modelResults = usedResults, baseResponse = baseResponses[iGroup], 
                  groupIndex = iGroup)) * increase
          
          if (bmr > CES) {
            
            return(0)
            
          } else {
            
            while ( abs(CES - bmr) >= 1.0e-10 & abs(top - bottom) >= 1.0e-10 &
                iterations <= nIterations) {
              
              tmpDose <- exp(mid) 
              bmr <- do.call(bmrFunction, list(weights = weights, dose = tmpDose, 
                      modelResults = usedResults, baseResponse = baseResponses[iGroup],
                      groupIndex = iGroup)) * increase
              
              if (bmr > CES)
                top <- mid else 
                bottom <- mid
              
              mid <- (top + bottom)/2 
              iterations <- iterations + 1
              
            }
            
            if (iterations > nIterations) 
              stop("The bmd could not be estimated after ", nIterations, " iterations.") else
              return(tmpDose) 
            
          }
          
        })
    
    names(bmdValues) <- names(baseResponses)
    
    return(bmdValues)
    
  }
  
}


#' Plot for the model-averaged response values and bmd
#' @param proastData list, data in proast format as returned by f.scan()
#' @param xans integer, indicates the column in proastData$data  
#' containing x-values (dose)
#' @param yans integer, indicates the column in proastData$data 
#' containing y-values (response)
#' @param nans integer, indicates the column in proastData$data 
#' containing sample size; if irrelevant set to 0; default value is 0
#' @param bmd numeric value, the estimated model-averaged BMD as returned by 
#' optimizeBmd() 
#' @param modelResults list, with results for fitted model for bmd; contains 
#' at least npar, loglik, model.ans, regr.par, CES and ces.ans; 
#' these are by default included in result from f.proast(); Eventually contains 
#' also fct1 and fct2 if factors are included for the model parameters
#' @param bootstrapBmd numeric vector, estimated bmd for each of the bootstrap data 
#' sets as returned by bootstrapBmd(); default value is NULL
#' @param bootstrapModelResults list with for each of the bootstrap data sets 
#' results of the fitted model; should contain at least see param modelResults above;
#' as returned by bootstrapBmd(); default value is NULL
#' @param confidenceLevel numeric, defines constructed confidence interval 
#' for bmd, e.g. if 0.9 then lower and upper bound of 90% CI for bmd are estimated;
#' default value is 0.9
#' @param naiveApproach boolean, TRUE if the model-averaged BMD is estimated as the 
#' weighted average of bmd values, FALSE if the model-averaged BMD is estimated
#' based on weighted average of response values; default value is FALSE
#' @return no return value; plot is written to the current device
#' @importFrom RColorBrewer brewer.pal
#' @export
plotAverageModel <- function(proastData, xans, yans, nans = 0, bmd, modelResults,
    bootstrapBmd = NULL, bootstrapModelResults = NULL, confidenceLevel = 0.9,
    naiveApproach = FALSE){
  
  xValues <- proastData$data[,xans]
  nEvents <- proastData$data[,yans]
  
  if (is.null(nans)) 
    nObs <- rep(1, length(nEvents)) else if (nans == 0) 
    nObs <- rep(1, length(nEvents)) else
    nObs <- proastData$data[,nans]
  
  yValues <- nEvents/nObs
  
  # Confidence limits of observed data
  confLimits <- sapply(seq_along(nObs), function(iGroup) {
        
        f.LL.bin(k = nEvents[iGroup], n = nObs[iGroup])
        
      })
  
  
  plotDose <- seq(min(xValues), max(xValues), length.out = 100)
  
  # Estimate response for all bootstrap data sets and original data
  tmpModelResults <- c(list(modelResults), bootstrapModelResults)
  
  plotResponses <- lapply(tmpModelResults, function(jFit) {
        
        # Calculate weights 
        aicValues <- as.numeric(sapply(jFit, function(iResult)            
                  2 * iResult$npar - 2 * iResult$loglik))
        weights <- calculateWeights(aicValues = aicValues)
        withError <- sapply(jFit, function(x) is(x, "error"))
        
        averageResponse(weights = weights, dose = plotDose, 
            modelResults = jFit[!withError], groupIndex = NA)
        
      })
  
  yLimits <- round(range(c(plotResponses, confLimits)), 2)
  
  # Plot empty
  plot(x = NULL, y = NULL, xlim = range(xValues), ylim = yLimits,
      xlab = proastData$varnames[xans], 
      ylab = proastData$varnames[yans], main = "Averaged response model")
  
  # Add fitted model lines
  # Multigroups?
  if (!is.null(ncol(plotResponses[[1]]))) {
    
    multiGroups <- TRUE
    colors <- brewer.pal(max(3, ncol(plotResponses[[1]])), "Dark2")
    
  } else multiGroups <- FALSE
  
  tmp <- sapply(length(plotResponses):1, function(i) {
        
        if (multiGroups) {
          
          sapply(seq_len(ncol(plotResponses[[i]])), function(j)
                lines(plotDose, plotResponses[[i]][,j], 
                    col = ifelse(i == 1, "black", colors[j])))
          
        } else {
          
          lines(plotDose, plotResponses[[i]], col = ifelse(i == 1, "black", "gray"))
          
        }
        
      })
  
  
  # Add reference lines for bmd, bmdl and bmdu
  if (!is.null(bootstrapBmd)) {
    
    bounds <- findBootstrapBounds(bootstrapBmd = bootstrapBmd, confidenceLevel = confidenceLevel)
    bmdValues <- c(list(bmd), bounds)
    
  } else {
    
    bmdValues <- list(bmd)
    
  }
  
  tmp <- sapply(length(bmdValues):1, function(i) {
        
        if (all(is.na(bmdValues[[i]])))
          return()
        
        sapply(seq_along(bmdValues[[i]]), function(iGroup) {
              
              if (i == 1) {
                
                crucialModelResults <- modelResults
                color <- "black"
                
              } else if (i == 2) {
                
                if (length(bounds$bmdl) == 1) {
                  
                  crucialModelResults <- bootstrapModelResults[[
                      which(as.numeric(bootstrapBmd) == bounds$bmdl)[1]]]
                  color <- "gray"
                  
                } else {
                  
                  crucialModelResults <- bootstrapModelResults[[
                      which(bootstrapBmd[, iGroup] == bounds$bmdl[iGroup])[1]]]
                  color <- colors[iGroup]
                  
                }
                
                
              } else {
                
                if (length(bounds$bmdu) == 1) {                  
                  
                  crucialModelResults <- bootstrapModelResults[[
                      which(bootstrapBmd == bounds$bmdu)[1]]]
                  color <- "gray"
                  
                } else {
                  
                  crucialModelResults <- bootstrapModelResults[[
                      which(bootstrapBmd[, iGroup] == bounds$bmdu[iGroup])[1]]]
                  color <- colors[iGroup]
                  
                }
                
              }
              
              # Calculate weights 
              aicValues <- as.numeric(sapply(crucialModelResults, function(iResult)            
                        2 * iResult$npar - 2 * iResult$loglik))
              weights <- calculateWeights(aicValues = aicValues)
              withError <- sapply(crucialModelResults, function(x) is(x, "error"))
              
              bmdResponse <- averageResponse(weights = weights, 
                  dose = bmdValues[[i]][iGroup], modelResults = crucialModelResults[!withError], 
                  groupIndex = iGroup)
              
              lines(c(min(xValues), bmdValues[[i]][iGroup]), 
                  rep(bmdResponse, 2), 
                  lty = 2, col = color)  # horizontal left-right
              lines(rep(bmdValues[[i]][iGroup], 2), 
                  c(bmdResponse, yLimits[1]), 
                  lty = 2, col = color)  # vertical up-down
            })
        
      })
  
  
  # Plot observed data
  points(xValues, yValues)
  
  # Plot confidence limits for observed data
  arrows(x0 = xValues, y0 = confLimits[1,], x1 = xValues, y1 = confLimits[2,],
      length = 0.05, angle = 90, code = 3)
  
  # Add legend if > 1 group
  if (multiGroups)
    legend("bottomright", legend = names(bmd), col = colors, lty = 1,
        bty = "n")
  
  
}



#' Estimate lower and upper bound for model-averaged bmd using parametric bootstrap 
#' @param proastData list, data in proast format as returned by f.scan()
#' @param weights numeric vector, estimated weights as returned by calculateWeights() 
#' @param modelResults list, with results for each model, same length as weights. 
#' For each model a list with at least npar, loglik,  model.ans, regr.par, CES and ces.ans; 
#' these are by default included in result from f.proast(). Eventually 
#' contains also fct1 and fct2 if factors are included for the model parameters.
#' @param shinyInput list with necessary parameters used for fitted models in 
#' modelResults 
#' @param naiveApproach boolean, TRUE if the model-averaged BMD is estimated as the 
#' weighted average of bmd values, FALSE if the model-averaged BMD is estimated
#' based on weighted average of response values; default value is FALSE
#' @param aicNull numeric, aic value for null model as criterion for accepting 
#' bootstrap data, if NA all bootstrap data are accepted; default value is NA
#' @param nBootstraps integer, the number of bootstrap data sets to generate;
#' default value is 200
#' @param seed integer, allows reproducing results; default value is 1
#' @param showProgress boolean, whether progress bar should be shown in shiny
#' application; important: only use this option when function is called from 
#' within shiny application; default value is FALSE
#' @return list with modelResults and bootstrapBmd. The modelResults contain 
#' modelResults for each bootstrap data set; bootstrapBmd is data frame with
#' all estimated bmd values per group
#' @importFrom shiny incProgress
#' @export
bootstrapBmd <- function(proastData, weights, modelResults, shinyInput, 
    naiveApproach = FALSE, aicNull = NA, nBootstraps = 200, seed = 1,
    showProgress = FALSE) {
  
  set.seed(seed)
  
  # To define groups
  if (is.null(shinyInput$fct1))
    namesFct1 <- " "  else if (shinyInput$fct1 == 0) 
    namesFct1 <- " " else 
    namesFct1 <- paste0(names(proastData$data)[shinyInput$fct1], ".", proastData$data[, shinyInput$fct1]) 
  
  if (is.null(shinyInput$fct1)) 
    namesFct2 <- " " else if (shinyInput$fct2 == 0)
    namesFct2 <- " " else
    namesFct2 <- paste0(names(proastData$data)[shinyInput$fct2], ".", proastData$data[, shinyInput$fct2]) 
  
  groupNames <- apply(cbind(namesFct1, namesFct2),
      1, function(x) paste(unique(x), collapse = " & "))
  groupNames <- gsub("  & | &  ", "", groupNames)
  
  
  # Generate bootstrap data
#  allGroups <- proastData$data[, c(shinyInput$fct1, shinyInput$fct2)]
#  uniqueGroups <- allGroups[!duplicated(allGroups), , drop = FALSE]
  
  # Parameters for condition to accept bootstrap data
  summaryTable <- bindModelResults(modelResults)
  bestModel <- findBestModel(summaryTable)$bestModel
  allModelNames <- getModelNames(dtype = shinyInput$dtype)
  selectedModel <- as.numeric(allModelNames[names(allModelNames) == 
              as.character(bestModel$model)])
  
  selectedShinyInput <- shinyInput
  allParameters <- getModelParameters(selectedModel = selectedModel,
      dtype = selectedShinyInput$dtype)
  bestParameters <- unlist(strsplit(bestModel$parameterNames, ", "))
  selectedParameters <- allParameters[which(names(allParameters) %in% bestParameters)]
  selectedCombination <- data.frame(fct1.no = 0, fct2.no = 0, fct3.no = 0,
      fct4.no = 0, fct5.no = 0)
  
  if (any(bestParameters != "")) {
    
    selectedCombination[, selectedParameters] <- 
        c(selectedShinyInput$fct1.no, selectedShinyInput$fct2.no,
            selectedShinyInput$fct3.no, selectedShinyInput$fct4.no,
            selectedShinyInput$fct5.no)[selectedParameters]
    
  } 
  
  selectedShinyInput[names(selectedCombination)] <- selectedCombination
  
  
  k <- 1
  nRemoved <- 0
  bootstrapData <- list()
  
  while(k <= nBootstraps) {
    
    dose <- proastData$data[, shinyInput$xans]
    averageProb <- averageResponse(weights = weights, dose = dose, 
        modelResults = modelResults, groupIndex = NA)
    
    if (is.null(dim(averageProb))) {
    
      if (is.null(shinyInput$nans))
        nObs <- rep(1, length(averageProb)) else if(shinyInput$nans == 0)
        nObs <- rep(1, length(averageProb)) else
        nObs <- as.numeric(proastData$data[, shinyInput$nans])
      
      proastData$data[, shinyInput$yans] <- sapply(seq_len(length(averageProb)), function(i)
            rbinom(1, nObs[i], averageProb[i]))
      
    } else {
      
      if (is.null(shinyInput$nans))
        nObs <- rep(1, nrow(averageProb)) else if(shinyInput$nans == 0)
        nObs <- rep(1, nrow(averageProb)) else
        nObs <- as.numeric(proastData$data[, shinyInput$nans])
      
      proastData$data[, shinyInput$yans] <- sapply(seq_len(nrow(averageProb)), function(i){
            rbinom(1, nObs[i], averageProb[i, groupNames[i]])
          })
      
    }
    
    
    if (is.na(aicNull)) {
      
      bootstrapData[[k]] <- proastData
      k <- k + 1
      
    } else {
      
      tmpResult <- fitSingleModel(data = proastData, 
          shinyInput = selectedShinyInput, selectedModel = selectedModel, 
          estimateCed = FALSE)
      newAic <- 2 * tmpResult$npar - 2 * tmpResult$loglik 
      
      if (newAic > (aicNull - 2)) {
        
        nRemoved <- nRemoved + 1
        next
        
      } else {
        
        bootstrapData[[k]] <- proastData
        k <- k + 1
        
      }
      
    }
    
  }
  
  cat("Number of removed bootstrap data sets:", nRemoved, "\n")
  
  # Investigate bootstrap data
  if (FALSE) {
    
    newResponses <- lapply(bootstrapData, function(x) x$data[, shinyInput$yans])
    origResponses <- proastData$data[, shinyInput$yans]
    plot(x = NULL, y= NULL, xlim = range(proastData$data[, shinyInput$xans]),
        ylim = range(newResponses, origResponses))
    sapply(newResponses, function(currentResponses)
          points(proastData$data[, shinyInput$xans], currentResponses, col = "grey"))
    points(proastData$data[, shinyInput$xans], origResponses, pch = 19)
    
  }
  
  
  modelIndices <- sapply(modelResults, function(x) x$model.ans)
  
  
  selectedShinyInput <- shinyInput
  fittedModels <- list()
  
  for (k in seq_along(bootstrapData)) {
    
    if (showProgress) {
      
      incProgress(1/nBootstraps, 
          detail = paste("bootstrap", k, "out of", nBootstraps))
      
    }
    
    # Fit all models for each bootstrap data set
    iData <- bootstrapData[[k]]
    
    fittedModels[[k]] <- lapply(modelResults, function(iModel){
          
          selectedShinyInput[c("fct1.no", "fct2.no", "fct3.no", "fct4.no", "fct5.no")] <- 
              iModel[c("fct1.no", "fct2.no", "fct3.no", "fct4.no", "fct5.no")]
          
          fitSingleModel(data = iData, shinyInput = selectedShinyInput,
              selectedModel = iModel$model.ans, estimateCed = naiveApproach)
          
        })
    
    # Estimate bmd for each bootstrap data set
    # Calculate weights 
    aicValues <- as.numeric(sapply(fittedModels[[k]], function(iResult)            
              2 * iResult$npar - 2 * iResult$loglik))
    weights <- calculateWeights(aicValues = aicValues)
    
    # Calculate BMD
    if (k == 1) {
      
      tmpBmd <- optimizeBmd(weights = weights, modelResults = fittedModels[[k]], 
          naiveApproach = naiveApproach)
      bootstrapBmd <- matrix(NA, nrow = nBootstraps, ncol = max(1, length(tmpBmd))) 
      bootstrapBmd[k, ] <- tmpBmd
      colnames(bootstrapBmd) <- names(tmpBmd)
      
    } else {
      
      bootstrapBmd[k,] <- optimizeBmd(weights = weights, modelResults = fittedModels[[k]], 
          naiveApproach = naiveApproach)
      
    }
    
  }
  
  
  return(
      list(
          modelResults = fittedModels,
          bootstrapBmd = bootstrapBmd
      )
  )
  
}


#' Find the lower and upper bound for bmd based on parametric bootstrap results
#' @param bootstrapBmd data.frame, bmd values based on parametric bootstrap; 
#' as returned by bootstrapBmd()
#' @param confidenceLevel numeric, defines constructed confidence interval 
#' for bmd, e.g. if 0.9 then lower and upper bound of 90% CI for bmd are estimated;
#' default value is 0.9
#' @return numeric vector, estimated lower and upper bound for model-averaged bmd  
#' @export
findBootstrapBounds <- function(bootstrapBmd, confidenceLevel = 0.9){
  
  if (ncol(bootstrapBmd) > 1) {
    
    bmdl <- c()
    bmdu <- c()
    
    for (iCol in 1:ncol(bootstrapBmd)) {
      
      bmdValues <- bootstrapBmd[, iCol]
      bmdValues <- bmdValues[!is.na(bmdValues)]
      nBootstraps <- length(bmdValues)
      orderedValues <- bmdValues[order(bmdValues)]
      
      bmdl <- c(bmdl, orderedValues[max(1, round((1 - confidenceLevel)/2*nBootstraps))])
      bmdu <- c(bmdu, orderedValues[min(nBootstraps, round((1 - (1 - confidenceLevel)/2)*nBootstraps))])
      
    }
    
    bmdBounds <- data.frame(bmdl = bmdl, bmdu = bmdu)
    rownames(bmdBounds) <- colnames(bootstrapBmd)
    
  } else {
    
    bootstrapBmd <- bootstrapBmd[!is.na(bootstrapBmd)]
    nBootstraps <- length(bootstrapBmd)
    
    orderedValues <- bootstrapBmd[order(bootstrapBmd)]
    bmdl <- orderedValues[max(1, round((1 - confidenceLevel)/2*nBootstraps))]
    bmdu <- orderedValues[min(nBootstraps, round((1 - (1 - confidenceLevel)/2)*nBootstraps))]
    
    bmdBounds <- data.frame(bmdl = bmdl, bmdu = bmdu)
    
  }
  
  
  return(bmdBounds)
  
}