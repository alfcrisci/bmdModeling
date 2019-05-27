# Fit all covariate combinations when a covariate is included
# Select best model given all covariate combinations
# Model-averaging when covariates are included 

library(bmdModeling)
library(testthat)

dataDir <- system.file("extdata", package = "bmdModeling")

track <- FALSE
doPrint <- FALSE
doPlot <- FALSE


## Quantal data ##

context("Proast - Quantal data (methyleug)")

metry <- f.scan(file.path(dataDir, "methyleug.txt"))

shinyInput <- list(
    dtype = 4, 
    cont = FALSE,
    xans = 1, 
    yans = 2, 
    nans = 4, 
    CES = 0.1,
    ces.ans = 3,
    conf.lev = 0.9
)


test_that("Different covariate combinations for Weibull model", {
      
      allFactors <- expand.grid(fct1.no = c(0, 5), fct2.no = c(0, 5), fct3.no = 0,
          fct4.no = c(0, 5))
      
      
      bestModels <- sapply(seq_len(nrow(allFactors)), function(i) {
            
            if (doPrint)
              print(i)
            
            shinyInput[paste0("fct", 1:ncol(allFactors), ".no")] <- allFactors[i, ]
            nCombinations <- 2^sum(allFactors[i, ] != 0)
            
            savedResults <- fitAllCovariates(data = metry, shinyInput = shinyInput, selectedModel = 19)
            
            if (nCombinations != 1)
              expect_equal(length(savedResults), nCombinations)
            summaryResult <- summaryModels(savedResults = savedResults)
            
            
            summaryTable <- summaryResult$summaryTable
            colnames(summaryTable) <- attr(summaryTable, "columnNames")
            expect_equal(nrow(summaryTable), nCombinations)
            
            if (doPrint)
              print(summaryTable)
            
            summaryResult$extraInfo$bestModel
            
          })
      
    })


test_that("All covariate combinations for all quantal models", {
      
      shinyInput$fct1.no <- 5
      shinyInput$fct2.no <- 5
      shinyInput$fct4.no <- 5
      
      
      savedResults <- fitAllModels(data = metry, shinyInput = shinyInput, 
          fitCovariateCombinations = TRUE)
      expect_equal(length(savedResults),
          sum(2^sapply(getModelNames(dtype = shinyInput$dtype), function(iModel){
                    modelParameters <- getModelParameters(selectedModel = iModel, dtype = shinyInput$dtype)
                    length(modelParameters[!is.na(modelParameters)])
                  })))
      
      summaryTable <- bindModelResults(savedResults = savedResults)
      extraInfo <- findBestModel(summaryTable = summaryTable)
      
      colnames(summaryTable) <- attr(summaryTable, "columnNames")
      if (doPrint)
        print(summaryTable)
      
      expect_equal(as.character(extraInfo$bestModel$model), "Gamma") 
      expect_equal(as.character(extraInfo$bestModel$parameterNames), "b")
      
    })


test_that("Model averaging with one covariate for quantal models", {
      
      shinyInput$fct1.no <- 5
      shinyInput$fct2.no <- 5
      shinyInput$fct4.no <- 5
      
      savedResults <- fitAllModels(data = metry, shinyInput = shinyInput, 
          fitCovariateCombinations = TRUE)
      
      summaryTable <- bindModelResults(savedResults = savedResults)
      bestSummaryTable <- filterBestCovariates(summaryTable = summaryTable)
      aicNull <- summaryTable[1, "aic"]
      
      if (doPrint) 
        print(bestSummaryTable)
      
      
      bestIndices <- matchWithResults(savedResults = savedResults,
          modelNames = bestSummaryTable$model[-c(1,2)], 
          parameterNames = bestSummaryTable$parameterNames[-c(1,2)])
#      bestIndices <- c(1, 2, 6, 7, 22)
      
      estimatedWeights <- calculateWeights(aicValues = bestSummaryTable$aic[-c(1,2)])
      
      # NON-NAIVE Calculation model averaged BMD
      maBmd <- optimizeBmd(weights = estimatedWeights, 
          modelResults = savedResults[bestIndices])
      expect_equal(as.numeric(signif(maBmd, 5)), c(2.1707, 3.0625))
      
      tmp <- sapply(seq_along(maBmd), function(iGroup) {
            
            baseResponse <- averageResponse(weights = estimatedWeights, dose = 0, 
                modelResults = savedResults[bestIndices], groupIndex = iGroup)
            estimatedCes <- bmr3(weights = estimatedWeights, dose = maBmd[iGroup], 
                modelResults = savedResults[bestIndices], 
                baseResponse = baseResponse, groupIndex = iGroup)
            
            expect_equal(round(as.numeric(estimatedCes), 2), 0.1)
            
          })
      
      bootstrapResults <- bootstrapBmd(proastData = metry, 
          weights = estimatedWeights, 
          modelResults = savedResults[bestIndices], 
          shinyInput = shinyInput, aicNull = aicNull, nBootstraps = 20)
      maBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
      expect_equivalent(signif(maBounds, 5), data.frame(bmdl = c(1.7733, 2.4206), 
              bmdu = c(2.8858, 4.1825)))
      
      
      tmp <- sapply(seq_along(maBmd), function(i) {
            expect_lt(maBounds$bmdl[i], maBmd[i]) 
            expect_lt(maBmd[i], maBounds$bmdu[i])      
          })
      
      if (doPlot)
        plotAverageModel(proastData = metry, xans = shinyInput$xans, 
            yans = shinyInput$yans, nans = shinyInput$nans, 
            modelResults = savedResults[bestIndices], bmd = maBmd,
            bootstrapBmd = bootstrapResults$bootstrapBmd,
            bootstrapModelResults = bootstrapResults$modelResults)
      
      # NAIVE Calculation model averaged BMD
      maSimpleBmd <- optimizeBmd(weights = estimatedWeights, 
          modelResults = savedResults[bestIndices], 
          naiveApproach = TRUE)
      expect_equal(as.numeric(signif(maSimpleBmd, 5)), c(2.1668, 3.0372))
      
      expect_equal(sum(bestSummaryTable$bmd.sex.1[-c(1,2)]*estimatedWeights, 
              na.rm = TRUE), as.numeric(maSimpleBmd[1]))
      expect_equal(sum(bestSummaryTable$bmd.sex.2[-c(1,2)]*estimatedWeights, 
              na.rm = TRUE), as.numeric(maSimpleBmd[2]))
      
      bootstrapResults <- bootstrapBmd(proastData = metry, 
          weights = estimatedWeights, 
          modelResults = savedResults[bestIndices], 
          shinyInput = shinyInput, aicNull = aicNull,
          nBootstraps = 20, naiveApproach = TRUE)
      maSimpleBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
      expect_equivalent(signif(maSimpleBounds, 5), data.frame(bmdl = c(1.7648, 2.3941), 
              bmdu = c(2.9596, 4.1384)))
      
      tmp <- sapply(seq_along(maSimpleBmd), function(i) {
            expect_lt(maSimpleBounds$bmdl[i], as.numeric(maSimpleBmd[i]))
            expect_lt(as.numeric(maSimpleBmd[i]), maSimpleBounds$bmdu[i])      
          })
      
    })


test_that("Model averaging with two covariates for quantal models", {
      
      newData <- metry
      
      newData$data <- rbind(metry$data, metry$data)
      newData$data$age <- rep(c(3, 4), each = nrow(newData$data)/2)
      newData$nvar <- 6
      newData$varnames <- c(metry$varnames, "age")
      newData$dtype <- c(metry$dtype, 0)
      
      # Random response values for second age group
      set.seed(1)
      newData$data[-(1:(nrow(newData$data)/2)), shinyInput$yans] <- 
          apply(newData$data[-(1:(nrow(newData$data)/2)),], 1, function(iRow){
                
                # Age group 2 has 10% larger probability of positive response
                nObs <- as.numeric(iRow[shinyInput$nans])
                prob <- min(1, as.numeric(iRow[shinyInput$yans]/nObs + 0.2))
                rbinom(1, nObs, prob)
                
              })
      
      if (doPlot)
        plot(newData$data[, shinyInput$xans], newData$data[, shinyInput$yans],
            col = as.numeric(interaction(newData$data[, 5], newData$data[, 6])))      
      
      shinyInput$fct1.no <- 5
      shinyInput$fct2.no <- 6
      shinyInput$fct4.no <- 5
      
      
      savedResults <- fitAllModels(data = newData, shinyInput = shinyInput, 
          fitCovariateCombinations = TRUE)
      
      summaryTable <- bindModelResults(savedResults = savedResults)
      aicNull <- summaryTable[1, "aic"]
      bestSummaryTable <- filterBestCovariates(summaryTable = summaryTable)
      
      if (doPrint) 
        print(bestSummaryTable)
      
      
      bestIndices <- matchWithResults(savedResults = savedResults,
          modelNames = bestSummaryTable$model[-c(1, 2)], 
          parameterNames = bestSummaryTable$parameterNames[-c(1, 2)])
      estimatedWeights <- calculateWeights(aicValues = bestSummaryTable$aic[-c(1, 2)])
      
#      bestIndices <- matchWithResults(savedResults = savedResults,
#          modelNames = c("Log-probit", "Log-probit", "Weibull"), 
#          parameterNames = c("a, b", "c", "a, b, c"))
#      estimatedWeights <- calculateWeights(aicValues = summaryTable$aic[c(15, 16, 27)])
      
      
      # NON-NAIVE Calculation model averaged BMD
      maBmd <- optimizeBmd(weights = estimatedWeights, 
          modelResults = savedResults[bestIndices])
      expect_equal(as.numeric(signif(maBmd, 5)), c(2.7738, 2.7764, 1.5483, 1.5496))
      
      
      tmp <- sapply(seq_along(maBmd), function(iGroup) {
            
            baseResponse <- averageResponse(weights = estimatedWeights, dose = 0, 
                modelResults = savedResults[bestIndices], groupIndex = iGroup)
            estimatedCes <- bmr3(weights = estimatedWeights, dose = maBmd[iGroup], 
                modelResults = savedResults[bestIndices], 
                baseResponse = baseResponse, groupIndex = iGroup)
            
            expect_equal(round(as.numeric(estimatedCes), 2), 0.1)
            
          })
      
      bootstrapResults <- bootstrapBmd(proastData = newData, 
          weights = estimatedWeights, 
          modelResults = savedResults[bestIndices], 
          shinyInput = shinyInput, aicNull = aicNull, nBootstraps = 20)
      maBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
      expect_equivalent(signif(maBounds, 5), data.frame(bmdl = c(2.3871, 2.3817, 1.4212, 1.4213), 
              bmdu = c(3.1565, 3.1883, 1.7395, 1.7395)))
      
      tmp <- sapply(seq_along(maBmd), function(i) {
            expect_lt(maBounds$bmdl[i], maBmd[i]) 
            expect_lt(maBmd[i], maBounds$bmdu[i])      
          })
      
      if (doPlot)
        plotAverageModel(proastData = newData, xans = shinyInput$xans, 
            yans = shinyInput$yans, nans = shinyInput$nans, 
            modelResults = savedResults[bestIndices], bmd = maBmd,
            bootstrapBmd = bootstrapResults$bootstrapBmd,
            bootstrapModelResults = bootstrapResults$modelResults)
      
      # NAIVE Calculation model averaged BMD
      maSimpleBmd <- optimizeBmd(weights = estimatedWeights, 
          modelResults = savedResults[bestIndices], naiveApproach = TRUE)
      expect_equal(as.numeric(signif(maSimpleBmd, 5)), c(2.7730, 1.5478))
      
      expect_equal(sum(bestSummaryTable$bmd.age.3[-c(1,2)]*estimatedWeights, 
              na.rm = TRUE), as.numeric(maSimpleBmd[1]))
      expect_equal(sum(bestSummaryTable$bmd.age.4[-c(1,2)]*estimatedWeights, 
              na.rm = TRUE), as.numeric(maSimpleBmd[2]))
      
      bootstrapResults <- bootstrapBmd(proastData = newData, 
          weights = estimatedWeights, 
          modelResults = savedResults[bestIndices], 
          shinyInput = shinyInput, aicNull = aicNull,
          nBootstraps = 20, naiveApproach = TRUE)
      maSimpleBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
      expect_equivalent(signif(maSimpleBounds, 5), data.frame(bmdl = c(2.3901, 1.4212),
              bmdu = c(3.1722, 1.7395)))
      
      tmp <- sapply(seq_along(maSimpleBmd), function(i) {
            expect_lt(maSimpleBounds$bmdl[i], as.numeric(maSimpleBmd[i]))
            expect_lt(as.numeric(maSimpleBmd[i]), maSimpleBounds$bmdu[i])      
          })
      
      # For comparison: non-naive BMD average models with naive bootstrap bounds
      if (doPlot)
        plotAverageModel(proastData = newData, xans = shinyInput$xans, 
            yans = shinyInput$yans, nans = shinyInput$nans, 
            modelResults = savedResults[bestIndices], bmd = maSimpleBmd,
            bootstrapBmd = bootstrapResults$bootstrapBmd,
            bootstrapModelResults = bootstrapResults$modelResults) 
      
      
      ## Bind BMD results with different dimensions: naive shorter than non-naive
      bmdTable1 <- cbind(maBmd, maBounds)
      bmdTable2 <- cbind(maSimpleBmd, maSimpleBounds)
      
      if (FALSE) {
        
        write.csv(newData$data, file = file.path(dataDir, "methyleug2cov.csv"),
            row.names = FALSE)
        
      }
      
    })


test_that("Model averaging with double covariates for quantal models", {
      
      rawData <- read.csv(file.path(dataDir, "methyleug2cov.csv"), header = TRUE)
      
      newData <- list()
      newData$info <- metry$info
      newData$nvar <- ncol(rawData)
      newData$varnames <- names(rawData)
      newData$data <- rawData
      newData$dtype <- c(metry$dtype, 0)
      
      fct1.no <- 5:6
      fct2.no <- 5:6
      fct3.no <- 0
      fct4.no <- 0
      fct5.no <- 0
      
      lastIndex <- ncol(rawData)
      
      for (i in 1:5) {
        
        if (length(get(paste0("fct", i, ".no"))) > 1) {
          
          selectedFactors <- rawData[, get(paste0("fct", i, ".no"))]
          names(selectedFactors) <- names(rawData)[get(paste0("fct", i, ".no"))]
          rawData <- cbind(rawData, 
              interaction(selectedFactors, drop = TRUE, sep = "_"))
          names(rawData)[lastIndex + 1] <- paste(names(rawData)[get(paste0("fct", i, ".no"))], collapse = "_")
          assign(paste0("fct", i, ".no"), lastIndex + 1)
          lastIndex <- lastIndex + 1
          
        }
        
      }
      
      newData$data <- rawData
      shinyInput[paste0("fct", 1:5, ".no")] <- mget(paste0("fct", 1:5, ".no"))
      
      savedResults <- fitAllModels(data = newData, shinyInput = shinyInput, 
          fitCovariateCombinations = TRUE)
      
      summaryTable <- bindModelResults(savedResults = savedResults)
      aicNull <- summaryTable[1, "aic"]
      bestSummaryTable <- filterBestCovariates(summaryTable = summaryTable)
      
      if (doPrint) 
        print(bestSummaryTable)
      
      
      bestIndices <- matchWithResults(savedResults = savedResults,
          modelNames = bestSummaryTable$model[-c(1, 2)], 
          parameterNames = bestSummaryTable$parameterNames[-c(1, 2)])
      estimatedWeights <- calculateWeights(aicValues = bestSummaryTable$aic[-c(1, 2)])
      
      
      
      # NON-NAIVE Calculation model averaged BMD
      maBmd <- optimizeBmd(weights = estimatedWeights, 
          modelResults = savedResults[bestIndices])
      expect_equal(as.numeric(signif(maBmd, 5)), c(2.2904, 3.4029, 1.4476, 2.2275))
      
      
      tmp <- sapply(seq_along(maBmd), function(iGroup) {
            
            baseResponse <- averageResponse(weights = estimatedWeights, dose = 0, 
                modelResults = savedResults[bestIndices], groupIndex = iGroup)
            estimatedCes <- bmr3(weights = estimatedWeights, dose = maBmd[iGroup], 
                modelResults = savedResults[bestIndices], 
                baseResponse = baseResponse, groupIndex = iGroup)
            
            expect_equal(round(as.numeric(estimatedCes), 2), 0.1)
            
          })
      
      bootstrapResults <- bootstrapBmd(proastData = newData, 
          weights = estimatedWeights, 
          modelResults = savedResults[bestIndices], 
          shinyInput = shinyInput, aicNull = aicNull, nBootstraps = 20)
      maBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
      expect_equivalent(signif(maBounds, 5), data.frame(bmdl = c(1.6128, 2.5714, 1.1568, 1.5468), 
              bmdu = c(2.9673, 4.1937, 1.8562, 3.2079)))
      
      tmp <- sapply(seq_along(maBmd), function(i) {
            expect_lt(maBounds$bmdl[i], maBmd[i]) 
            expect_lt(maBmd[i], maBounds$bmdu[i])      
          })
      
      if (doPlot)
        plotAverageModel(proastData = newData, xans = shinyInput$xans, 
            yans = shinyInput$yans, nans = shinyInput$nans, 
            modelResults = savedResults[bestIndices], bmd = maBmd,
            bootstrapBmd = bootstrapResults$bootstrapBmd,
            bootstrapModelResults = bootstrapResults$modelResults)
      
      # NAIVE Calculation model averaged BMD
      maSimpleBmd <- optimizeBmd(weights = estimatedWeights, 
          modelResults = savedResults[bestIndices], naiveApproach = TRUE)
      expect_equal(as.numeric(signif(maSimpleBmd, 5)), c(2.2871, 3.3843, 1.4401, 2.2166))
      
      expect_equal(sum(bestSummaryTable$bmd.sex_age.1_3[-c(1,2)]*estimatedWeights, 
              na.rm = TRUE), as.numeric(maSimpleBmd[1]))
      expect_equal(sum(bestSummaryTable$bmd.sex_age.2_3[-c(1,2)]*estimatedWeights, 
              na.rm = TRUE), as.numeric(maSimpleBmd[2]))
      expect_equal(sum(bestSummaryTable$bmd.sex_age.1_4[-c(1,2)]*estimatedWeights, 
              na.rm = TRUE), as.numeric(maSimpleBmd[3]))
      expect_equal(sum(bestSummaryTable$bmd.sex_age.2_4[-c(1,2)]*estimatedWeights, 
              na.rm = TRUE), as.numeric(maSimpleBmd[4]))
      
      # Takes very long time to fit
      if (FALSE) {
        bootstrapResults <- bootstrapBmd(proastData = newData, 
            weights = estimatedWeights, 
            modelResults = savedResults[bestIndices], 
            shinyInput = shinyInput, aicNull = aicNull,
            nBootstraps = 20, naiveApproach = TRUE)
        maSimpleBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
        expect_equivalent(signif(maSimpleBounds, 5), data.frame(bmdl = c(2.3901, 1.4212),
                bmdu = c(3.1722, 1.7395)))
        
        tmp <- sapply(seq_along(maSimpleBmd), function(i) {
              expect_lt(maSimpleBounds$bmdl[i], as.numeric(maSimpleBmd[i]))
              expect_lt(as.numeric(maSimpleBmd[i]), maSimpleBounds$bmdu[i])      
            })
      }
      
    })



context("Proast - Continuous data (das1)")

load(file.path(dataDir, "das1.rda"))

shinyInput <- list(
    dtype = 1, 
    xans = 1, 
    yans = 4, 
    CES = 0.05, 
    cont = TRUE,
    conf.lev = 0.9,
    sf.x = 100)

# Remove missing response values
das1$data <- das1$data[!is.na(das1$data$BW), ]


test_that("All covariate combinations for all continuous models", {
      
      shinyInput$fct1.no <- 5
      shinyInput$fct2.no <- 5
      shinyInput$fct3.no <- 0 
      shinyInput$fct4.no <- 0 
      shinyInput$fct5.no <- 0 
      # Non-zero for the others takes long time to fit all combinations
      
      savedResults <- fitAllModels(data = das1, shinyInput = shinyInput, 
          fitCovariateCombinations = TRUE)
      expect_equal(length(savedResults), sum(2^c(0, 0, 2, 2, 2, 2)))
      
      summaryTable <- bindModelResults(savedResults = savedResults)
      extraInfo <- findBestModel(summaryTable = summaryTable)
      
      colnames(summaryTable) <- attr(summaryTable, "columnNames")
      if (doPrint)
        print(summaryTable)
      
      expect_equal(as.character(extraInfo$bestModel$model), "Hill model 3") 
      expect_equal(as.character(extraInfo$bestModel$parameterNames), "a, b")
      
    })


if (FALSE) {
  
  # Check results for das10 cont summary data: not converged for E3 with b and H3 with b
  load(file.path(dataDir, "das10.rda"))
  
  shinyInput <- list(dtype = 10, 
      xans = 2, 
      # mean response
      yans = 4, 
      # associated sd
      sans = 5, 
      # sample size
      nans = 6, 
      # extra parameters
      sd.se = 1, 
      # critical effect size
      CES = 0.05, 
      # is model continous
      cont = TRUE)
  
  
  test_that("All covariate combinations for all continuous models", {
        
        shinyInput$fct1.no <- 3
        shinyInput$fct2.no <- 3
        shinyInput$fct3.no <- 0
        shinyInput$fct4.no <- 0
        shinyInput$fct5.no <- 0
        # Non-zero for the others takes long time to fit all combinations
        
        savedResults <- fitAllModels(data = das10, shinyInput = shinyInput, 
            fitCovariateCombinations = TRUE)
        expect_equal(length(savedResults),
            sum(2^sapply(getModelNames(dtype = shinyInput$dtype), function(iModel)
                      length(getModelParameters(selectedModel = iModel, dtype = shinyInput$dtype)))))
        
        summaryTable <- bindModelResults(savedResults = savedResults)
        extraInfo <- findBestModel(summaryTable = summaryTable)
        
        colnames(summaryTable) <- attr(summaryTable, "columnNames")
        if (doPrint)
          print(summaryTable)
        
        expect_equal(as.character(extraInfo$bestModel$model), "Hill model 3") 
        expect_equal(as.character(extraInfo$bestModel$parameterNames), "a")
        
      })
  
  
  
# Find error in application for specific data
  rawData <- read.table("/home/mvarewyck/Documents/bmd/Background/data2008.csv",
      sep = ",", header = TRUE)
  
  # Remove problematic record
  rawData <- rawData[-63, ]
  
  proastData <- list(info = "newData", 
      nvar = ncol(rawData),
      varnames = colnames(rawData),
      data = rawData,
      dtype = rep(0, ncol(rawData)))
  
  shinyInput <- list(dtype = 4, 
      xans = 2, yans = 6, nans = 5,
      CES = 0.1, ces.ans = 3, cont = FALSE, conf.lev = 0.9, sf.x = 1)
  
  # Fit all models
  results <- fitAllModels(data = proastData, shinyInput = shinyInput)
  summaryTable <- bindModelResults(savedResults = results)
  aicNull <- summaryTable[1, "aic"]
  extraInfo <- findBestModel(summaryTable = summaryTable)
  summaryTable
  if (doPlot)
    toPlot <- f.plot.all(results[[7]])
  
  tmp <- summaryModels(savedResults = results)
  
  # Model Averaging
  # exclude non-fitted or non-converged models
  estimatedWeights <- calculateWeights(aicValues = summaryTable$aic[-c(1, 2, 8)])
  
  bestIndices <- matchWithResults(savedResults = results,
      modelNames = summaryTable$model[-c(1, 2, 8)], 
      parameterNames = summaryTable$parameterNames[-c(1, 2, 8)])
  
  # Non-naive approach for calculating model averaged BMD
  maBmd <- optimizeBmd(weights = estimatedWeights, 
      modelResults = results[bestIndices])
  
  bootstrapResults <- bootstrapBmd(proastData = proastData, 
      weights = estimatedWeights, 
      modelResults = results[bestIndices], 
      shinyInput = shinyInput, aicNull = aicNull,
      nBootstraps = 100)
  maBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
  
  # Plot results
  if (doPlot)
    plotAverageModel(proastData = proastData, xans = shinyInput$xans, 
        yans = shinyInput$yans, nans = shinyInput$nans, 
        modelResults = results[bestIndices], 
        bmd = maBmd,
        bootstrapBmd = bootstrapResults$bootstrapBmd, 
        bootstrapModelResults = bootstrapResults$modelResults)
  
  # Naive approach for calculating model averaged BMD
  maSimpleBmd <- optimizeBmd(weights = estimatedWeights, 
      modelResults = results[bestIndices], naiveApproach = TRUE)
  
  bootstrapResults <- bootstrapBmd(proastData = proastData, 
      weights = estimatedWeights, 
      modelResults = results[bestIndices], 
      shinyInput = shinyInput, aicNull = aicNull,
      nBootstraps = 10, naiveApproach = TRUE)
  maSimpleBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
  
  
}
