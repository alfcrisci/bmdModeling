library(bmdModeling)
library(testthat)

dataDir <- system.file("extdata", package = "bmdModeling")

track <- FALSE


context("Proast - Quantal data (methyleug)")

metry <- f.scan(file.path(dataDir, "methyleug.txt"))

shinyInput <- list(
    # quantal
    dtype = 4, 
    cont = FALSE,
    # xans = independent variable = dose
    xans = 1, 
    # response: number of responding animals: forest.kk, liver.kk
    yans = 2, # also 3
    # associated sample size
    nans = 4, # sample.size
    # need for f.mm6.cat
    CES = 0.1,
    ces.ans = 3,
    # confidence level
    conf.lev = 0.9
)


test_that("No model averaging, 1 response", {
      
      savedResults <- fitAllModels(data = metry, shinyInput = shinyInput, 
          track = track)
      summaryResult <- summaryModels(savedResults = savedResults)
      
      bmdData <- data.frame(response = summaryResult$extraInfo$responseName, 
#          bmd = summaryResult$summaryTable$bmd[summaryResult$extraInfo$bestModelIndex],
          bmdl = summaryResult$extraInfo$minBmdl, 
          bmdu = summaryResult$extraInfo$maxBmdu)
      bmdData
      #    response     bmdl     bmdu
      # 1 forest.kk 1.538738 3.512874
      
      makeForestPlot(bmdData = bmdData)
      
    })


test_that("No model averaging, multiple responses", {
      
      bmdValues <- lapply(c(2, 3), function(yans) {
            
            shinyInput$yans <- yans
            savedResults <- fitAllModels(data = metry, shinyInput = shinyInput, 
                track = track)
            summaryResult <- summaryModels(savedResults = savedResults)
            
            data.frame(response = summaryResult$extraInfo$responseName, 
#                bmd = summaryResult$summaryTable$bmd[summaryResult$extraInfo$bestModelIndex],
                bmdl = summaryResult$extraInfo$minBmdl, 
                bmdu = summaryResult$extraInfo$maxBmdu)
            
          })
      
      bmdData <- do.call(rbind, bmdValues)
      bmdData
      #    response     bmdl     bmdu
      # 1 forest.kk 1.538738 3.512874
      # 2  liver.kk 2.558796 5.931440

      
#      bmdText <- c(bmd = "Best model BMD", bmdl = "Lowest BMDL", 
#          bmdu = "Highest BMDU")
      
      makeForestPlot(bmdData = bmdData)
      
    })


test_that("With model averaging", {
      
      bmdValues <- lapply(c(2, 3), function(yans) {
            
            shinyInput$yans <- yans
            fittedModels <- fitAllModels(data = metry, shinyInput = shinyInput, 
                track = track)
            
            aicValues <- sapply(fittedModels, function(iResult)            
                  2 * iResult$npar - 2 * iResult$loglik)
            estimatedWeights <- calculateWeights(aicValues = aicValues)
            bmd <- optimizeBmd(weights = estimatedWeights, 
                modelResults = fittedModels)
            
            nullModel <- fitSingleModel(data = metry, shinyInput = shinyInput, 
                selectedModel = 1)
            aicNull <- 2 * nullModel$npar - 2 * nullModel$loglik 
            
            bootstrapResults <- bootstrapBmd(proastData = metry, 
                weights = estimatedWeights, 
                modelResults = fittedModels, 
                shinyInput = shinyInput, aicNull = aicNull, nBootstraps = 20)
            bounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
            
            data.frame(response = metry$varnames[yans], 
                bmd = bmd,
                bmdl = bounds["bmdl"], 
                bmdu = bounds["bmdu"])
            
          })
      
      bmdData <- do.call(rbind, bmdValues)
      bmdData
      #    response      bmd     bmdl     bmdu
      # 1 forest.kk 2.756860 2.383461 8.156733
      # 2  liver.kk 4.602876 3.791038 6.492105

#      bmdText <- c(bmd = "Model-averaged BMD", bmdl = "Model-averaged BMDL", 
#          bmdu = "Model-averaged BMDU")
      
      makeForestPlot(bmdData = bmdData)
      
    })


test_that("With all covariate combinations", {
      
      shinyInput$fct1.no <- 5
      shinyInput$fct2.no <- 5
      
      bmdValues <- lapply(c(2, 3), function(yans) {
            
            shinyInput$yans <- yans
            savedResults <- fitAllModels(data = metry, shinyInput = shinyInput, 
                fitCovariateCombinations = TRUE)
            summaryResult <- summaryModels(savedResults = savedResults)
            covariateNames <<- sub(pattern = "bmd.", replacement = "", 
                names(summaryResult$summaryTable)[grepl(pattern = "bmd.", names(summaryResult$summaryTable), fixed = TRUE)])
            
            data.frame(response = summaryResult$extraInfo$responseName,
#                summaryResult$extraInfo$bestBmd,
                as.list(summaryResult$extraInfo$minBmdl), 
                as.list(summaryResult$extraInfo$maxBmdu))            
            
          })
      
      bmdData <- do.call(rbind, bmdValues)
      bmdData
      #                                   response bmdl.sex_1 bmdl.sex_2 bmdu.sex_1
      # summaryResult.extraInfo.minBmdl  forest.kk   1.136994   1.690616   3.016317
      # summaryResult.extraInfo.minBmdl1  liver.kk   4.559233   1.629791  14.036617
      #                                  bmdu.sex_2
      # summaryResult.extraInfo.minBmdl    4.587193
      # summaryResult.extraInfo.minBmdl1   5.320930

      
      makeForestPlot(bmdData = bmdData, groupNames = covariateNames)
      
    })


test_that("With model averaging and with covariates", {
      
      shinyInput$fct1.no <- 5
      shinyInput$fct2.no <- 5
      
      bmdValues <- lapply(c(2, 3), function(yans) {
            
            shinyInput$yans <- yans
            savedResults <- fitAllModels(data = metry, shinyInput = shinyInput, 
                track = track, fitCovariateCombinations = TRUE)
            
            summaryTable <- bindModelResults(savedResults = savedResults)
            bestSummaryTable <- filterBestCovariates(summaryTable = summaryTable)
            aicNull <- summaryTable[1, "aic"]
            
            bestIndices <- matchWithResults(savedResults = savedResults,
                modelNames = bestSummaryTable$model[-c(1,2)], 
                parameterNames = bestSummaryTable$parameterNames[-c(1,2)])
            
            estimatedWeights <- calculateWeights(aicValues = bestSummaryTable$aic[-c(1,2)])

            maBmd <- optimizeBmd(weights = estimatedWeights, 
                modelResults = savedResults[bestIndices])
            
            bootstrapResults <- bootstrapBmd(proastData = metry, 
                weights = estimatedWeights, 
                modelResults = savedResults[bestIndices], 
                shinyInput = shinyInput, aicNull = aicNull, nBootstraps = 20)
            maBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
            
            data.frame(response = metry$varnames[yans], 
                bmd = as.list(maBmd),
                bmdl = t(maBounds[, "bmdl"]), 
                bmdu = t(maBounds[, "bmdu"]))
            
          })
      
      bmdData <- do.call(rbind, bmdValues)
      bmdData
      #    response bmd.sex_1 bmd.sex_2   bmdl.1   bmdl.2    bmdu.1   bmdu.2
      # 1 forest.kk  2.170676  3.062535 1.773318 2.420598  2.885813 4.182478
      # 2  liver.kk  7.280856  3.074729 5.953864 1.992275 11.166782 4.036783

      
#      bmdText <- c(bmd = "Model-averaged BMD", bmdl = "Model-averaged BMDL", 
#          bmdu = "Model-averaged BMDU")
      
      makeForestPlot(bmdData = bmdData, groupNames = covariateNames)
      
    })



test_that("Single model with covariates", {
      
      shinyInput$fct1.no <- 5
      shinyInput$fct2.no <- 5

      bmdValues <- lapply(c(2, 3), function(yans) {
            
            shinyInput$yans <- yans
            savedResults <- fitSingleModel(data = metry, shinyInput = shinyInput, 
                selectedModel = 21)
            summaryResult <- summaryModels(savedResults = savedResults)
            
            data.frame(response = summaryResult$extraInfo$responseName,
                summaryResult$summaryTable[grepl("bmd.", names(summaryResult$summaryTable), fixed = TRUE)],
                summaryResult$summaryTable[grepl("bmdl.", names(summaryResult$summaryTable), fixed = TRUE)], 
                summaryResult$summaryTable[grepl("bmdu.", names(summaryResult$summaryTable), fixed = TRUE)])
            
          })
      
      bmdData <- do.call(rbind, bmdValues)
      bmdData
      #    response bmd.sex_1 bmd.sex_2 bmdl.sex_1 bmdl.sex_2 bmdu.sex_1 bmdu.sex_2
      # 1 forest.kk  2.419465  3.453864   1.895329   2.661489   2.987808   4.587193
      # 2  liver.kk 11.535510  4.086500   5.614385   2.362691  19.573449   5.709021

      
      makeForestPlot(bmdData = bmdData, groupNames = covariateNames)
      
    })


      
# Reported bug by Jose
if (FALSE) {
  
  data <- read.table("/home/mvarewyck/Documents/bmd/Background/RatsFullW14.txt",
      header = TRUE)
  proastData <- list(
      info = "testData",
      nvar = ncol(data),
      varnames = names(data),
      data = data,
      dtype = rep(0, ncol(data))) 
  
  shinyInput <- list(
      dtype = 10, 
      cont = TRUE,
      xans = 1, 
      yans = 4, 
      nans = 8,
      sans = 6,
      sd.se = 1,
      CES = 0.1,
      conf.lev = 0.9,
      fct1.no = 2,
      fct2.no = 2
  )
  
  savedResults <- fitAllModels(data = proastData, shinyInput = shinyInput)
  summaryResult <- summaryModels(savedResults = savedResults)
  
  bmdData <- data.frame(response = summaryResult$extraInfo$responseName,
      as.list(summaryResult$extraInfo$minBmdl), 
      as.list(summaryResult$extraInfo$maxBmdu))
  
  covariateNames <- sub(pattern = "bmd.", replacement = "", 
      names(summaryResult$summaryTable)[
          grepl(pattern = "bmd.", names(summaryResult$summaryTable), fixed = TRUE)])
  
  bmdData
  #    response     bmdl     bmdu
  # 1 forest.kk 1.538738 3.512874
  
  makeForestPlot(bmdData = bmdData, groupNames = covariateNames)
        
}
