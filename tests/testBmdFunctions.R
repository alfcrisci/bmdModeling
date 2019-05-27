
library(bmdModeling)
library(testthat)

dataDir <- system.file("extdata", package = "bmdModeling")

# Compare results with original f.proast()
#source("/home/mvarewyck/git/bmd/exploringProast/fFunctionsProast.R")
#load("/home/mvarewyck/git/bmd/proast61.3/data/das5.rda")
#f.proast(das5)
#track <- TRUE
#track2 <- TRUE

track <- FALSE
# track <- TRUE



## Load raw or proast data ##

# Data proast
fileName <- file.path(dataDir, "methyleug.txt")
separator <- "\t"

# Raw data
fileName <- file.path(dataDir, "methyleug_raw.txt")
separator <- "\t"

tmpData <- read.table(fileName, header = TRUE, 
    sep = separator, quote = '"', stringsAsFactors = FALSE)

if (any(apply(tmpData, 2, class) == "character"))
  dataType <- "proastData" else
  dataType <- "rawData"

dataType

# Fix decimal commas
tmpData$forest.kk <- tmpData$forest.kk/10

newData <- as.data.frame(sapply(tmpData, function(x) sub("\\.", ",", as.character(x))),
    stringsAsFactors = FALSE)
newData
originalData <- as.data.frame(sapply(tmpData, function(x) as.numeric(sub(",", "\\.", as.character(x)))))
summary(originalData)


## Continuous data ##

context("Proast - Continuous data (das1)")

load(file.path(dataDir, "das1.rda"))
#write.csv(das1$data, file = "~/git/bmd/proastUI/inst/extdata/das1.csv", row.names = FALSE)

shinyInput <- list(
    #continuous data
    dtype = 1, quick.ans = 1, 
    #fit model
    main.ans = 4, 
    # dose
    xans = 1, 
    # mean response
    yans = 4, 
    # critical effect size
    CES = 0.05, 
    # is model continous
    cont = TRUE,
    # construct <conf.lev>% confidence intervals 
    conf.lev = 0.9,
    #scaling factor
    sf.x = 100)


# Remove missing response values
das1$data <- das1$data[!is.na(das1$data$BW), ]


test_that("Run f.proast() - fit model for all continuous models", {
      
      shinyInput$main.ans <- 4		
      
      result <- lapply(c(1, 11, 13, 15, 23, 25), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            
            shinyInput$model.ans <- model.ans
            
            tmpResult <- f.proast(odt = das1, shinyInput = shinyInput, track = track)
            
            if(model.ans < 13)
              tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
            
            return(tmpResult)
            
          })
      
      
      expect_equal(sapply(result, function(x) signif(x$loglik, 5)), 
          c(-25.57, -19.51, -19.58, -19.58, -19.58, -19.58))
      
      expect_equal(sapply(result, function(x) {
                values <- signif(x$MLE, 5)
                names(values) <- NULL
                values
              }),
          list(c(0.082333, 341.41),
              c(0.075949, 359.47, 351.97, 358.11, 294.79),
              c(0.076015, 356.42, 42.112, 10),
              c(0.076015, 356.42, 41.578, 0.58753, 10),
              c(0.076015, 356.42, 41.815, 10),
              c(0.076015, 356.42, 41.204, 0.45569, 10)
          ))
      
    })


test_that("Run f.proast() - calculate CED and plot for all models beside null&full", {
      
      shinyInput$main.ans <- c(4, 6)		
      
      result <- lapply(c(13, 15, 23, 25), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            
            shinyInput$model.ans <- model.ans
            
            result <- f.proast(odt = das1, shinyInput = shinyInput, track = track)
            tmp <- f.plot.all(ans.all = result, track = track)
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$CED.origScale, 5)),
          c(4211.2, 4157.8, 4181.5, 4120.4))
      
      expect_equal(sapply(result, function(x) signif(x$conf.int, 5)), 
          matrix(c(0, 438.31, 322.67, 439.85,
                  4485.7, 4485.8, 4474.4, 4473.9), nrow = 2, byrow = TRUE))
      
    })


test_that("Using summary function in performAnalysis.R - continuous response", {
      
      fittedModels <- fitAllModels(data = das1, shinyInput = shinyInput, 
          track = track)
      summaryResult <- summaryModels(savedResults = fittedModels)
      tmpTable <- summaryResult$summaryTable
      print(tmpTable)
      
      tmpPlot <- sapply(summaryResult$extraInfo$bestModelIndex, function(iBest)
            f.plot.all(fittedModels[[iBest]], track = track)) 
          
    })



## Quantal data ##

context("Proast - Quantal data (methyleug)")

metry <- f.scan(file.path(dataDir, "methyleug.txt"))

shinyInput <- list(
    # quantal
    dtype = 4, cont = FALSE,
    # model chosen (1==NULL)
    model.ans = 1,
    quick.ans = 1, #always 1: single model
    # main menu
    main.ans = 4, # fit.model 
    # xans = independent variable = dose
    xans = 1, 
    # response: number of responding animals: forest.kk, liver.kk
    yans = 2, # also 3
    # associated sample size
    nans = 4, # sample.size (or nans?)
    # type of benchmark response (1 = ED50)
    ces.ans = 1,
    # potential covariate (not yet implemented)
#    covar.no = 5, # sex
    # needed?		
    sd.se = 1, 
    # need for f.mm6.cat
    CES = 0.05
)


test_that("Run f.proast() - run model for all quantal models", {
      
      result <- lapply(c(1, 14, 16, 18, 19, 21, 24, 25, 26), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            
            shinyInput$model.ans <- model.ans
            
            # issue with ces.ans == 4
            allResults <- lapply(1:4, function(ces.ans){# 1:4
                  
                  if(!(ces.ans == 4 & model.ans %in% c(16:26))){
                    
                    cat("Ces.ans = ", ces.ans, "\n")
                    
                    shinyInput$ces.ans <- ces.ans		
                    tmpResult <- f.proast(odt = metry, shinyInput = shinyInput, track = track)
                    
                    if(model.ans < 16)
                      tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
                    
                    return(tmpResult)
                    
                  } else	cat("Model not tested, bug in the code")
                  
                })
            
            cat("\n\n")
            
            return(allResults[[1]])
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$loglik, 5)), 
          c(-275.98, -112.78, -115.06, -113.59, -113.81, -113.01, -113.05, -122.32, -117.75))
      
      expect_equal(sapply(result, function(x) {
                values <- signif(x$MLE, 5)
                names(values) <- NULL
                values
              }),
          list(c(0.46),
              c(0.01, 0.14, 0.70, 0.99),
              c(0.0087307, 7.6881, 2.307),
              c(0.010956, 6.7232, 2.4838),
              c(0.008928, 7.4141, 1.4954),
              c(0.010429, 6.6753, 1.4339),
              c(0.0095945, 7.2174, 2.1732),
              c(-1.6546, 8.4740),
              c(-3.1294, 8.0979)))
      
    })



test_that("Run f.proast() - calculate CED for quantal models", {
      
      shinyInput$main.ans <- c(4, 6)		
      
      result <- lapply(c(16, 18, 19, 21, 24, 25, 26), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            
            shinyInput$model.ans <- model.ans
            shinyInput$conf.lev <- 0.9
            
            # issue with ces.ans == 4
            allResults <- lapply(1:3, function(ces.ans){# 1:4
                  
                  cat("Ces.ans = ", ces.ans, "\n")
                  
                  shinyInput$ces.ans <- ces.ans		
                  tmpResult <- f.proast(odt = metry, shinyInput = shinyInput, track = track)
                  tmp <- f.plot.all(ans.all = tmpResult, track = track)
                  
                  return(tmpResult)
                  
                })
            
            cat("\n\n")
            return(allResults[[1]])
            
          })
      
      print(sapply(result, function(x) signif(x$CED.matr, 5)))
      print(sapply(result, function(x) signif(x$conf.int, 5)))
      
      expect_equal(sapply(result, function(x) signif(x$CED.matr, 5)),
          c(7.6881, 6.7232, 7.4141, 6.6753, 7.2174, 8.474, 8.0979))
      
      expect_equal(sapply(result, function(x) signif(x$conf.int, 5)), 
          matrix(c(6.6612, 5.943, 6.498, 5.9392, 6.3977, 0.012249, 7.385,
                  8.7164, 7.5847, 8.3948, 7.4899, 8.111, Inf, 8.8707), nrow = 2, byrow = TRUE))
      
    })



test_that("Using summary function in performAnalysis.R - quantal response", {
      
      shinyInput$conf.lev <- 0.9
      fittedModels <- fitAllModels(data = metry, shinyInput = shinyInput, track = track)
      summaryResult <- summaryModels(savedResults = fittedModels)
      tmpTable <- summaryResult$summaryTable
      print(tmpTable)
      
      tmpPlot <- sapply(summaryResult$extraInfo$bestModelIndex, function(iBest)
            f.plot.all(fittedModels[[iBest]], track = track)) 
      
    })



## Continuous summary data ##

context("Proast - Continuous summary data (atra)")

atra <- f.scan(file.path(dataDir, "atra.txt"))

test_that("Load data", {
      
      expect_match(atra$info, "atrazine")
      expect_equal(names(atra), c("info", "nvar", "varnames", "data", "dtype"))
      
    })

shinyInput <- list(dtype = 10, quick.ans = 1, 
    main.ans = 4, 
    # continous model
    model.ans = 1, 
    # dose
    xans = 1, 
    # mean response
    yans = 2, 
    # associated sd
    sans = 3, 
    # sample size
    nans = 4, 
    # extra parameters
    sd.se = 1, 
    # critical effect size
    CES = 0.05, 
    # is model continous
    cont = TRUE)


test_that("Run f.proast() - fit model for all continuous models", {
      
      shinyInput$main.ans <- 4
      result <- lapply(c(1, 11, 13, 15, 23, 25), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            shinyInput$model.ans <- model.ans
            tmpResult <- f.proast(odt = atra, shinyInput = shinyInput, track = track)
            
            if(model.ans < 13)
              tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
            
            return(tmpResult)
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$loglik, 5)), 
          c(-11.16, 4.73, 4.15, 4.27, 4.16, 4.27))
      
      expect_equal(sapply(result, function(x) signif(x$MLE, 5)),
          list(c(0.068009, 437.790000),
              c(0.05495, 485.11000, 459.57000, 468.54000, 388.52000, 352.22000),
              c(0.055379, 480.020000, 3.545400, 0.688400),
              c(0.055286, 473.730000, 8.348500, 0.732930, 1.535500),
              c(0.055367, 479.100000, 3.984300, 0.771080),
              c(0.055288, 473.450000, 8.251600, 0.687830, 1.763600)))
      
    })


test_that("Run f.proast() - calculate CED with CI for Exp and Hill models", {
      
      shinyInput$main.ans <- c(4, 6)
      shinyInput$conf.lev <- 0.9
      
      result <- lapply(c(13, 15, 23, 25), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            shinyInput$model.ans <- model.ans
            
            tmpResult <- f.proast(odt = atra, shinyInput = shinyInput,
                track = track)
            f.plot.all(ans.all = tmpResult, track = track)
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$CED, 5)),
          c(3.5454, 8.3485, 3.9843, 8.2516))
      
      expect_equal(sapply(result, function(x) signif(x$conf.int, 5)), 
          matrix(c(0.17464, 0.20767, 0.22209, 0.26154,
                  14.665, 26.745, 15.061, 22.382), nrow = 2, byrow = TRUE))
      
    })


test_that("Using summary function in performAnalysis.R - continuous response", {
      
      shinyInput <- list(dtype = 10, xans = 1, yans = 2, nans = 4, sans = 3, 
          sd.se = 1, CES = 0.05, cont = TRUE, conf.lev = 0.9)
      
      fittedModels <- fitAllModels(data = atra, shinyInput = shinyInput)
      summaryResult <- summaryModels(savedResults = fittedModels)
      tmpTable <- summaryResult$summaryTable
      print(tmpTable)
      
      tmpPlot <- sapply(summaryResult$extraInfo$bestModelIndex, function(iBest)
            f.plot.all(fittedModels[[iBest]], track = track)) 
      
    })


## Binary data ##

context("Proast - Binary data (das2)")

load(file.path(dataDir, "das2.rda"))

shinyInput <- list(
    # binary
    dtype = 2, cont = FALSE,
    # model chosen (1==NULL)
    model.ans = 1,
    quick.ans = 1, #always 1: single model
    # main menu
    main.ans = 4, # fit.model 
    # xans = independent variable = dose
    xans = 1, 
    # response: number of responding animals: forest.kk, liver.kk
    yans = 2, # also 3
    # type of benchmark response (1 = ED50)
    ces.ans = 1,
    # needed?		
    sd.se = 1, 
    # give scaling factor for dose >
    sf.x = 100, 
    # need for f.mm6.cat
    CES = 0.05
)

test_that("Run f.proast() - run model for all binary models", {
      
      result <- lapply(c(1, 14, 16, 18, 19, 21, 24, 25, 26), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            
            shinyInput$model.ans <- model.ans
            shinyInput$main.ans <- 4
            
            # issue with ces.ans == 4
            allResults <- lapply(1:4, function(ces.ans){# 1:4
                  
                  if(!(ces.ans == 4 & model.ans %in% c(16:26))){	# object 'CED' not found							
                    
                    cat("Ces.ans = ", ces.ans, "\n")
                    
                    shinyInput$ces.ans <- ces.ans		
                    tmpResult <- f.proast(odt = das2, shinyInput = shinyInput, track = track)
                    
                    # Error for plot if model.ans == 14, also in proast61.3
                    if(model.ans < 16)
                      tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
                    
                    return(tmpResult)
                    
                  } else	cat("Model not tested, bug in the code")
                  
                })
            
            cat("\n\n")
            
            return(allResults[[1]])
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$loglik, 5)), 
          c(-49.42, 0, -46.39, -46.32, -46.39, -46.29, -46.38, -46.62, -46.66))
      
      expect_equal(sapply(result, function(x) {
                values <- signif(x$MLE, 5)
                names(values) <- NULL
                values
              }),
          list(c(0.22826),
              ifelse(das2$data$EHC[order(das2$data$CRD)], 1, 1e-06),
              c(0.084321, 9.9202, 1e-06),
              c(0.11012, 9.3087, 1.4065),
              c(0.092085, 9.752200, 1.059400),
              c(0.12916, 9.1819, 0.91092),
              c(0.10038, 9.60710, 1.16900),
              c(-1.1275, 9.7093),
              c(-1.8610, 9.6704)
          ))
      
    })

test_that("Run f.proast() - calculate CED for binary models", {
      
      shinyInput$main.ans <- c(4, 6)		
      
      result <- lapply(c(16, 18, 19, 21, 24, 25, 26), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            
            shinyInput$model.ans <- model.ans
            shinyInput$conf.lev <- 0.9
            
            # issue with ces.ans == 4
            allResults <- lapply(1:3, function(ces.ans){# 1:4
                  
                  cat("Ces.ans = ", ces.ans, "\n")		
                  
                  # get an issue if sf.x == 1: need at least two non-NA values to interpolate
                  # but not when use scaling factor, e.g. 100
                  if(model.ans == 25 & ces.ans == 2){
                    
                    cat("Model not tested, no convergence f.profile.all\n")
                    
                  } else {				
                    
                    shinyInput$ces.ans <- ces.ans		
                    tmpResult <- f.proast(odt = das2, shinyInput = shinyInput, track = FALSE)
                    f.plot.all(ans.all = tmpResult, track = track)
                    
                  }
                  
                })
            
            cat("\n\n")
            
            return(allResults[[1]])
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$CED.matr*x$sf.x, 5)),
          c(992.02, 930.87, 975.22, 918.19, 960.71, 970.93, 967.04))
      
      expect_equal(sapply(result, function(x) signif(x$conf.int, 5)), 
          matrix(c(591.34, 557.19, 590.75, 558.76, 588.74, 3.1541, 655.34,
                  2846.3, 19904, 18057, 22093, 12902, Inf, 2452.6), nrow = 2, byrow = TRUE))
      
    })


test_that("Using summary function in performAnalysis.R - binary response", {
      
      shinyInput$conf.lev <- 0.9
      shinyInput$ces.ans <- 1
      
      fittedModels <- fitAllModels(data = das2, shinyInput = shinyInput,
          track = track)
      summaryResult <- summaryModels(savedResults = fittedModels)
      tmpTable <- summaryResult$summaryTable
      print(tmpTable)
      
      tmpPlot <- sapply(summaryResult$extraInfo$bestModelIndex, function(iBest)
            f.plot.all(fittedModels[[iBest]], track = track)) 
      
    })


## Ordinal data ##

context("Proast - Ordinal data (das3)")

# load data
load(file.path(dataDir, "das3.rda"))
dataComplete <- f.remove.NAs(variableIndices = 
        # xans, Vyans, sans, covar.no, nans
        as.numeric(c(2,	3, NULL, 1, NULL)), 
    originalData = das3$data, track = track)
das3$data <- dataComplete

shinyInput <- list(
    dtype = 3, quick.ans = 1, 
    # dose
    xans = 2, 
    # mean response
    yans = 3, 
    # extra parameters
    CES.cat = 2, # severity category associated with CED: to add in interface
    sf.x = 1,
    # critical effect size
    CES = 0.05, 
    # is model continous?
    cont = FALSE)


test_that("Run f.proast() - fit model for all ordinal models", {
      
      shinyInput$main.ans <- 4
      result <- lapply(c(1, 14, 13, 15, 23, 25), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            shinyInput$model.ans <- model.ans
            
            tmpResult <- f.proast(odt = das3, shinyInput = shinyInput, track = track)
            if(model.ans %in% c(1, 14))
              tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
            
            return(tmpResult)
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$loglik, 5)), 
          c(-91.26, -38.52, -42.21, -38.48, -37.87, -36.84))
      
      expect_equal(sapply(result, function(x) {
                values <- signif(x$MLE, 5)
                names(values) <- NULL
                values
              }),
          list(c(1.0596, 0, -0.19547, -0.48297, -0.42343, -0.16694, 1),
              c(10, 8.3564, 0.0012483, 0, -1.557, -1.6617, -0.94695, -0.32788, 1),
              c(10, 8.0424, 0.56669, 0, -0.90447, -1.4741, -0.9943, -0.43925, 1),
              c(10, 8.4325, 0.0012194, 1.0758, 0, -1.622, -1.6728, -0.93934, -0.32123, 1),
              c(10, 7.8662, 2.6782, 0, -1.443, -1.6444, -1.0106, -0.41943, 1),
              c(10, 8.1266, 0.00049196, 3.8767, 0, -2.1147, -1.6666, -0.96341, -0.35389, 1)
          ))
      
    })


test_that("Run f.proast() - calculate CED for ordinal models and fit results", {
      
      shinyInput$main.ans <- c(4, 6)		
      
      result <- lapply(c(13, 15, 23, 25), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            
            shinyInput$model.ans <- model.ans
            shinyInput$conf.lev <- 0.9	
            
            tmpResult <- f.proast(odt = das3, shinyInput = shinyInput, track = track)
            
            # Plot for all fitted models
            sapply(c(3, 5), function(plot.type){
                  
                  tmpPlot <- f.plot.all(ans.all = tmpResult, track = track, plot.type = plot.type)
                  
                })
            
            tmpResult
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$CED.matr, 5))[2,],
          c(8.0424, 8.4325, 7.8662, 8.1266))
      
      expect_equal(sapply(result, function(x) signif(x$conf.int, 5)), 
          matrix(c(6.0618, 7.0032, 6.2713, 6.6274,
                  10.208, 9.6048, 9.5717, 9.4547), nrow = 2, byrow = TRUE))
      
    })


test_that("Using summary function in performAnalysis.R - quantal response", {
      
      shinyInput$conf.lev <- 0.9
      fittedModels <- fitAllModels(data = das3, shinyInput = shinyInput, track = track)
      summaryResult <- summaryModels(savedResults = fittedModels)
      tmpTable <- summaryResult$summaryTable
      print(tmpTable)
      
      tmpPlot <- sapply(summaryResult$extraInfo$bestModelIndex, function(iBest)
            f.plot.all(fittedModels[[iBest]], track = track, plot.type = 3)) 
      
    })


## Continuous clustered data ##

context("Proast - Continuous clustered data (das5)")

load(file.path(dataDir, "das5.rda"))
# write.csv(das5$data, file = "~/git/bmd/proastUI/inst/extdata/das5.csv", row.names = FALSE)

shinyInput <- list(
    #continuous clustered data
    dtype = 5, quick.ans = 1, 
    # dose
    xans = 1, 
    # mean response
    yans = 11,
    # nested factor
    nest.no = 10,    
    # critical effect size
    CES = 0.05, 
    # is model continous
    cont = TRUE,
    # construct <conf.lev>% confidence intervals 
    conf.lev = 0.9,
    #scaling factor
    sf.x = 1000)


test_that("Run f.proast() - fit all models", {
      
      shinyInput$main.ans <- 4
      result <- lapply(c(1, 11, 13, 15, 23, 25), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            shinyInput$model.ans <- model.ans
            tmpResult <- f.proast(odt = das5, shinyInput = shinyInput, track = track)
            
            # Plotting takes very long time
#            if(model.ans < 13)
#              tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
            
            return(tmpResult)
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$inter.var, 5)), 
          c(0.00374910, 0.00013444, 0.00025598, 0.00019029, 0.00023903, 0.00019289))
      
      expect_equal(sapply(result, function(x) signif(x$intra.var, 5)), 
          c(rep(0.009204, 6)))
      
    })


test_that("Run f.proast() - calculate CED with CI for Exp and Hill models", {
      
      shinyInput$main.ans <- c(4, 6, 7)
      shinyInput$conf.lev <- 0.9
      shinyInput$nruns <- 30
      
      result <- lapply(c(13, 15, 23, 25), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            shinyInput$model.ans <- model.ans
            tmpResult <- f.proast(odt = das5, shinyInput = shinyInput, track = track)
            # Plotting takes a long time
#            tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
            
            return(tmpResult)
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$sf.x*x$CED, 5)),
          c(478.82, 537.69, 488.70, 542.43))
      
      expect_equal(sapply(result, function(x) signif(x$conf.int, 4)), 
          matrix(c(357.9, 430.5, 373.0, 433.2,
                  557.1, 609.8, 561.6, 612.9), nrow = 2, byrow = TRUE))
      
    })


test_that("Using summary function in performAnalysis.R - continuous response", {
      
      shinyInput$conf.lev <- 0.9
      shinyInput$nruns <- 30
      
      fittedModels <- fitAllModels(data = das5, shinyInput = shinyInput)
      summaryResult <- summaryModels(savedResults = fittedModels)
      tmpTable <- summaryResult$summaryTable
      print(tmpTable)
      
      # Plotting takes a long time
#      tmpPlot <- sapply(summaryResult$extraInfo$bestModelIndex, function(iBest)
#            f.plot.all(fittedModels[[iBest]], track = track)) 
      
    })





## Quantal clustered data ##

context("Proast - Quantal clustered data (das6)")

load(file.path(dataDir, "das6.rda"))
#write.csv(das6$data, file = "~/git/bmd/proastUI/inst/extdata/das6.csv", row.names = FALSE)

shinyInput <- list(
    #continuous clustered data
    dtype = 6, quick.ans = 1, 
    # dose
    xans = 1, 
    # mean response
    yans = 3,
    # associated sample sizes
    nans = 2,    
    # critical effect size
    CES = 0.05, 
    # is model continous
    cont = FALSE,
    # construct <conf.lev>% confidence intervals 
    conf.lev = 0.9,
    #scaling factor
    sf.x = 100)


test_that("Run f.proast() - run model for all quantal models", {
      
      shinyInput$main.ans <- 4
      
      result <- lapply(c(1, 14, 16, 18, 19, 21, 24, 25, 26), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            
            shinyInput$model.ans <- model.ans
            
            # issue with ces.ans == 4: not tested, bug in code
            allResults <- lapply(1:3, function(ces.ans){
                  
                  if(!(ces.ans %in% c(2,3) & model.ans == 25)){
                    
                    cat("Ces.ans = ", ces.ans, "\n")
                    
                    shinyInput$ces.ans <- ces.ans		
                    tmpResult <- f.proast(odt = das6, shinyInput = shinyInput, track = track)
                    
                    if(model.ans < 16)
                      tmpPlot <- f.plot.all(ans.all = tmpResult, track = track)
                    
                    return(tmpResult)
                    
                  } else	cat("Model not tested, need to adapt start values \n")
                  
                })
            
            cat("\n\n")
            
            return(allResults[[1]])
            
          })
      
      
      expect_equal(sapply(result, function(x) signif(x$loglik, 5)), 
          c(-366.28, -310.69, -312.83, -311.57, -310.85, -312.04, -310.74,
              -319.21, -318.8))
      
      expect_equal(sapply(result, function(x) {
                values <- signif(x$MLE, 5)
                names(values) <- NULL
                values
              }),
          list(c(0.16527, 0.31006),
              c(0.36705, 0.037457, 0.070363, 0.1022, 0.25231, 0.76562, 0.91074),
              c(0.34429, 0.06508, 0.26307, 1e-06),
              c(0.35685, 0.048434, 0.16206, 0.96091),
              c(0.3651, 0.040463, 0.19914, 0.69617),
              c(0.35214, 0.054682, 0.15906, 0.57272),
              c(0.36608, 0.038479, 0.21673, 0.60401),
              c(0.29263, -1.25470, 0.47269),
              c(0.29697, -2.1411, 0.44615)))
      
    })



test_that("Run f.proast() - calculate CED for quantal models", {
      
      shinyInput$main.ans <- c(4, 6)		
      
      result <- lapply(c(16, 18, 19, 21, 24, 25, 26), function(model.ans){
            
            cat("Model.ans = ", model.ans, "\n")
            
            shinyInput$model.ans <- model.ans
            shinyInput$conf.lev <- 0.9
            
            # issue with ces.ans == 4: bug in the code
            allResults <- lapply(1:3, function(ces.ans){# 1:4
                  
                  cat("Ces.ans = ", ces.ans, "\n")
                  
                  if(!(ces.ans %in% c(2,3) & model.ans == 25)){
                    
                    shinyInput$ces.ans <- ces.ans		
                    tmpResult <- f.proast(odt = das6, shinyInput = shinyInput, track = track)
                    tmp <- f.plot.all(ans.all = tmpResult, track = track)
                    
                    return(tmpResult)
                    
                  } else	cat("Model not tested, need to adapt start values \n")
                  
                  
                })
            
            cat("\n\n")
            return(allResults[[1]])
            
          })
      
      expect_equal(sapply(result, function(x) signif(x$CED.matr*x$sf.x, 5)),
          c(26.307, 16.206, 19.914, 15.906, 21.673, 47.269, 44.615))
      
      expect_equal(sapply(result, function(x) signif(x$conf.int, 5)), 
          matrix(c(19.808, 10.179, 13.431, 9.6881, 15.254, 19.823, 36.352,
                  35.591, 25.711, 29.627, 25.698, 31.119, 43.459, 54.868), 
              nrow = 2, byrow = TRUE))
      
    })


test_that("Using summary function in performAnalysis.R - quantal response", {
      
      shinyInput$conf.lev <- 0.9
      shinyInput$ces.ans <- 1
      
      fittedModels <- fitAllModels(data = das6, shinyInput = shinyInput, track = track)
      summaryResult <- summaryModels(savedResults = fittedModels)
      tmpTable <- summaryResult$summaryTable
      print(tmpTable)
      
      tmpPlot <- sapply(summaryResult$extraInfo$bestModelIndex, function(iBest)
            f.plot.all(fittedModels[[iBest]], track = track)) 
      
    })



## Continuous summary data, clustered ##

# NOTE: This data type results in errors with proast61.3 => not yet implemented


if(FALSE){
  
  context("Proast - Continuous summary data, clustered (das10)")
  
  load(file.path(dataDir, "das10.rda"))
  #write.csv(das10$data, file = "~/git/bmd/proastUI/inst/extdata/das10.csv", row.names = FALSE)
  
  
  ## Need at least 3 groups for clustered analysis?
#das10$data <- rbind(das10$data, das10$data)
#das10$data[11:20, "sex"] <- c(rep(3, 5), rep(4, 5))
#das10$data[11:20, "mean"] <- round(rnorm(10, mean(das10$data$mean), mean(das10$data$sd)))
#das10$data[11:20, "sd"] <- round(rnorm(10, mean(das10$data$sd), sd(das10$data$sd)), 1)
#das10$data[11:20, "number"] <- round(rnorm(10, mean(das10$data$number), sd(das10$data$number)))
#
#das10$data
  
  
  shinyInput <- list(dtype = 15, quick.ans = 1, 
      # dose
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
      cont = TRUE,
      nest.no = 3)
  
  
  test_that("Run f.proast() - fit model for all continuous models", {
        
        shinyInput$main.ans <- 4
        result <- lapply(c(1, 11, 13, 15, 23, 25), function(model.ans){
              
              cat("Model.ans = ", model.ans, "\n")
              shinyInput$model.ans <- model.ans
              tmpResult <- f.proast(odt = das10, shinyInput = shinyInput, track = track)
              
            })
        
      })
  
}