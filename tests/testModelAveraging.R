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
    # quantal
    dtype = 4, cont = FALSE,
    # xans = independent variable = dose
    xans = 1, 
    # response: number of responding animals: forest.kk, liver.kk
    yans = 2, # also 3
    # associated sample size
    nans = 4, # sample.size
    # need for f.mm6.cat
    CES = 0.05,
    # confidence level
    conf.lev = 0.9
)

nullModel <- fitSingleModel(data = metry, shinyInput = shinyInput, 
    selectedModel = 1)
aicNull <- 2 * nullModel$npar - 2 * nullModel$loglik 


allCed <- matrix(NA, nrow = 7 , ncol = 3)
allModelResults <- list()

for (ces.ans in 1:3) {
  
  shinyInput$ces.ans <- ces.ans
  modelIndices <- getModelNames(dtype = shinyInput$dtype)[-c(1,2)]
  allModelResults[[ces.ans]] <- fitAllModels(data = metry, shinyInput = shinyInput, 
      modelIndices = modelIndices, track = track)
  
  allCed[, ces.ans] <- sapply(allModelResults[[ces.ans]], function(tmpResult){
        
        predictedResponses <- f.expect.bin(model.ans = tmpResult$model.ans, 
            x = tmpResult$x, regr.par = tmpResult$regr.par, 
            CES = tmpResult$CES, ces.ans = tmpResult$ces.ans)
        
        if(doPrint)
          cat("predicted responses:", predictedResponses, "\n")
        
        currentCed <- tmpResult$CED.matr * tmpResult$sf.x
        if(doPrint)
          cat("estimated CED:", currentCed, "\n")
        
        
        return(currentCed)        
        
      })
  
}
allCed
#          [,1]      [,2]      [,3]
# [1,] 8.097868 2.1707373 2.1050923
# [2,] 8.473965 1.8813381 1.8086467
# [3,] 6.723154 2.0827486 2.0730446
# [4,] 6.675262 2.1469867 2.1393841
# [5,] 7.414059 1.3193275 1.3112422
# [6,] 7.217400 1.6887146 1.6801213
# [7,] 7.688083 0.9722054 0.9642627



test_that("Estimate weights based on fitted models for quantal response, all ces.ans", {
      
      allWeights <- sapply(allModelResults, function(iModelResults){
            
            aicValues <- sapply(iModelResults, function(iResult)            
                  2 * iResult$npar - 2 * iResult$loglik)
            calculateWeights(aicValues = aicValues)
            
          })
      
      expect_equal(apply(allWeights, 2, sum), rep(1, 3))
      # Note: weights do not depend on ces.ans
      
    })


# Calculate weights 
aicValues <- sapply(allModelResults[[1]], function(iResult)            
      2 * iResult$npar - 2 * iResult$loglik)
estimatedWeights <- calculateWeights(aicValues = aicValues)


test_that("Run f.proast() - model-averaged BMD for quantal response", {
      
      ## For ces.ans = 1
      # Calculate averaged response for range of doses
      averagedResponses <- averageResponse(weights = estimatedWeights, 
          dose = c(0, 3, 10, 30), modelResults = allModelResults[[1]])
      expect_equal(signif(as.numeric(averagedResponses), 5), 
          c(0.010223, 0.143880, 0.696190, 0.988220))
      
      
      ## For ces.ans = 1:3
      # Calculate ma-BMD
      averagedBmd <- sapply(1:3, function(ces.ans){
            
            # Calculate model averaged BMD
            bmd <- optimizeBmd(weights = estimatedWeights, 
                modelResults = allModelResults[[ces.ans]])
            
            # Calculate model averaged BMD
            bmdNaive <- optimizeBmd(weights = estimatedWeights, 
                modelResults = allModelResults[[ces.ans]], naiveApproach = TRUE)
            
            # Plot the model-averaged dose response
#            png(paste0("averagedModel", ces.ans, ".png", width = 800, height = 600)
            if (doPlot)
              plotAverageModel(proastData = metry, xans = shinyInput$xans, 
                  yans = shinyInput$yans, nans = shinyInput$nans, 
                  modelResults = allModelResults[[ces.ans]], bmd = bmd)
#            dev.off()
            
            return(c(bmd = bmd, bmdNaive = bmdNaive))
            
          })
      
      
      # BMD following weighted average-response
      expect_equal(signif(averagedBmd["bmd",], 5), c(6.9860, 1.8264, 1.8172))
      
      # BMD following weighted average-bmd (naive approach)
      expect_equal(signif(averagedBmd["bmdNaive",], 5), c(7.0097, 1.8271, 1.8183))
      
    })



test_that("Run f.proast() - 1-model averaged BMD for quantal response", {
      
      allBmd <- sapply(1:7, function(iModel) {
            
            if (doPrint)
              cat("* model.ans:", c(16, 18, 19, 21, 24, 25, 26)[iModel], "*\n")
            
            # Define weights
            weights <- rep(0, 7)
            weights[iModel] <- 1
            
            # Calculate BMD
            averagedBmd <- sapply(1:3, function(ces.ans)
                  optimizeBmd(weights = weights, 
                      modelResults = allModelResults[[ces.ans]]))
            
            if(doPrint)
              cat("estimated bmd:", averagedBmd, "\n")
            
            
            return(averagedBmd)
            
          })
      
      expect_equal(signif(t(allBmd), 5), signif(allCed, 5))
      
    })



test_that("Estimate BMDL using parametric bootstrap - quantal response", {
      
      averagedBmdl <- sapply(1:3, function(ces.ans){
            
            if(doPrint)
              cat("ces.ans:", ces.ans, "\n")
            
            shinyInput$ces.ans <- ces.ans
            
            # Calculate model averaged BMD
            bmd <- optimizeBmd(weights = estimatedWeights, 
                modelResults = allModelResults[[ces.ans]])
            
            # Calculate model averaged BMDL and BMDU
            bootstrapResults <- bootstrapBmd(proastData = metry, 
                weights = estimatedWeights, 
                modelResults = allModelResults[[ces.ans]], 
                shinyInput = shinyInput, aicNull = aicNull, nBootstraps = 20)
            bounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
            
            expect_lt(bounds$bmdl, bmd)
            expect_lt(bmd, bounds$bmdu)
            
            # Plot results
            if (doPlot)
              plotAverageModel(proastData = metry, xans = shinyInput$xans, 
                  yans = shinyInput$yans, nans = shinyInput$nans, 
                  modelResults = allModelResults[[ces.ans]], 
                  bmd = bmd,
                  bootstrapBmd = bootstrapResults$bootstrapBmd, 
                  bootstrapModelResults = bootstrapResults$modelResults)
            
            
            ## NAIVE approach
            # Calculate model averaged BMD
            bmdNaive <- optimizeBmd(weights = estimatedWeights, 
                modelResults = allModelResults[[ces.ans]], naiveApproach = TRUE)
            
            # Calculate model averaged BMDL
            bootstrapResultsNaive <- bootstrapBmd(proastData = metry, 
                weights = estimatedWeights, 
                modelResults = allModelResults[[ces.ans]], 
                shinyInput = shinyInput, aicNull = aicNull,
                nBootstraps = 20, naiveApproach = TRUE)
            naiveBounds <- findBootstrapBounds(bootstrapBmd = bootstrapResultsNaive$bootstrapBmd)
            
            expect_lt(naiveBounds$bmdl, bmdNaive)
            expect_lt(bmdNaive, naiveBounds$bmdu)
            
            
            return(c(bmdl = as.numeric(bounds$bmdl), 
                    bmdlNaive = as.numeric(naiveBounds$bmdl)))
            
          })
      
      expect_equal(signif(averagedBmdl["bmdl",], 5), c(6.4194, 1.3035, 1.3035))
      expect_equal(signif(averagedBmdl["bmdlNaive",], 5), c(6.4261, 1.3150, 1.3150))
      
    })




if(FALSE){
  # Takes a lot of time
  
  test_that("Time for parametric bootstrap - quantal response", {
        
        bmdl <- sapply(1:3, function(ces.ans){
              
              if(doPrint)
                cat("ces.ans:", ces.ans, "\n")
              
              shinyInput$ces.ans <- ces.ans
              
              # Calculate model averaged BMDL
#              sapply(c(20, 100, 999), function(nBootstraps){
              sapply(c(20), function(nBootstraps){
                    
                    time1 <- Sys.time()
                    bootstrapResults <- bootstrapBmd(proastData = metry, weights = estimatedWeights, 
                        modelResults = allModelResults[[ces.ans]], 
                        shinyInput = shinyInput, 
                        nBootstraps = nBootstraps)
                    bounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
                    time2 <- Sys.time()
                    
                    if(doPrint){
                      cat("nBootstraps:", nBootstraps, "\n")
                      print(time2 - time1)
                    }
                    
                  })
              
              return(bounds["bmdl"])
              
            })
        
      })
}


## Model averaging for other data types ## 

# Binary data

test_that("Model averaging for binary data", {
      
      load(file.path(dataDir, "das2.rda"))
#      write.csv(das2$data, file.path(tempdir(), "das2.csv"))
      
      shinyInput <- list(
          # binary
          dtype = 2, cont = FALSE,
          # xans = independent variable = dose
          xans = 1, 
          # response
          yans = 2, 
          # type of benchmark response (1 = ED50)
          ces.ans = 3,
          # give scaling factor for dose
          sf.x = 1, 
          # need for f.mm6.cat
          CES = 0.1,
          conf.lev = 0.9
      )
      
      modelIndices <- getModelNames(dtype = shinyInput$dtype)[-c(1,2)]
      modelResults <- fitAllModels(data = das2, shinyInput = shinyInput, estimateCed = TRUE,
          modelIndices = modelIndices)
      allCed <- sapply(modelResults, function(iResult) iResult$CED)
      allCed
      
      # Calculate weights 
      aicValues <- unlist(sapply(modelResults, function(iResult)            
                2 * iResult$npar - 2 * iResult$loglik))
      estimatedWeights <- calculateWeights(aicValues = aicValues)
      
      # Calculate model averaged BMD
      bmd <- optimizeBmd(weights = estimatedWeights, 
          modelResults = modelResults)
      
      expect_equal(signif(bmd, 5), 247.78) 
      mean(unlist(allCed), na.rm = TRUE)
      # [1] 242.7994
      
      nullModel <- fitSingleModel(data = das2, shinyInput = shinyInput, 
          selectedModel = 1)
      aicNull <- 2 * nullModel$npar - 2 * nullModel$loglik 
      
      bootstrapResults <- bootstrapBmd(proastData = das2, 
          weights = estimatedWeights, 
          modelResults = modelResults, 
          shinyInput = shinyInput, aicNull = aicNull, nBootstraps = 20)
      bounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
      bounds
      #       bmdl     bmdu
      # 1 114.8877 769.8622
      
    })


# TODO no correct results, waiting for Wout
if (FALSE) {
  
  test_that("Model averaging for quantal clustered data", {
        
        load(file.path(dataDir, "das6.rda"))
#      write.csv(das6$data, file.path(tempdir(), "das6.csv"))
        
        shinyInput <- list(
            dtype = 6, # with litter effect
#          dtype = 4, # without litter effect
            # dose
            xans = 1, 
            # mean response
            yans = 3,
            # associated sample sizes
            nans = 2,    
            # critical effect size
            CES = 0.1, 
            # is model continous
            cont = FALSE,
            # construct <conf.lev>% confidence intervals 
            conf.lev = 0.9,
            #scaling factor
            sf.x = 1,
            ces.ans = 3)
        
        
        modelIndices <- getModelNames(dtype = shinyInput$dtype)[-c(1,2)]
        modelResults <- fitAllModels(data = das6, shinyInput = shinyInput, estimateCed = TRUE,
            modelIndices = modelIndices)
        allCed <- sapply(modelResults, function(iResult) iResult$CED)
        allCed
        bindModelResults(modelResults)
        
        # Calculate weights 
        aicValues <- unlist(sapply(modelResults, function(iResult)            
                  2 * iResult$npar - 2 * iResult$loglik)[-2])
        estimatedWeights <- calculateWeights(aicValues = aicValues)
        
        # Calculate model averaged BMD
        bmd <- optimizeBmd(weights = estimatedWeights, modelResults = modelResults)
        expect_equal(signif(bmd, 5), 1.612)
        
        # Calculate model averaged BMDL and BMDU
        nullShinyInput <- shinyInput
#      nullShinyInput$dtype <- 4  # no litter effect for null model
        nullModel <- fitSingleModel(data = das6, shinyInput = nullShinyInput, 
            selectedModel = 1)
        aicNull <- 2 * nullModel$npar - 2 * nullModel$loglik 
        
        bootstrapResults <- bootstrapBmd(proastData = das6, 
            weights = estimatedWeights, 
            modelResults = modelResults[-2], 
            shinyInput = shinyInput, aicNull = aicNull, nBootstraps = 20)
        bounds <- findBootstrapBounds(bootstrapBmd = bootstrapResults$bootstrapBmd)
        
      })
  
}