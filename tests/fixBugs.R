# Project: bmd_git
# 
# Author: mvarewyck
###############################################################################

# Compare results with original f.proast()
#source("/home/mvarewyck/git/bmd/exploringProast/fFunctionsProast.R")
#dataDir <- "~/git/bmd/proastUI/inst/extdata"
#load(file.path(dataDir, "das1.rda"))
#f.proast(das1)


if (FALSE) {
  
  dataDir <- "/home/mvarewyck/Documents/bmd/Background"
  
  # Data 1 #
  rawData <- read.table(file.path(dataDir, "Marco.txt"), header = TRUE)
  proastData <- list(info = "testData", 
      nvar = ncol(rawData),
      varnames = colnames(rawData),
      data = rawData,
      dtype = rep(0, ncol(rawData)))
  
  shinyInput <- list(dtype = 4, 
      xans = 1, 
      yans = 2, 
      nans = 3, 
      CES = 0.1,
      ces.ans = 3,
      cont = FALSE)
  
  savedResults <- fitSingleModel(data = proastData, shinyInput = shinyInput, 
      selectedModel = 25, track = TRUE)
  bindModelResults(savedResults = savedResults)
  #    model parameterNames npar loglik   aic     bmd     bmdl     bmdu converged
  # 1 Probit                   2  -3.14 10.28 74.5686 57.91975 85.49365      TRUE
  
# Expected result (Proast 62.10)
# npar = 2.00, loglik =	-3.14, AIC = 10.28, accepted = yes, BMDL = 57.90, BMDU = 85.50, BMD = 74.60 
  
  # Data 2 #
  rawData <- read.table(file.path(dataDir, "examplequantal.txt"), header = TRUE)
  proastData <- list(info = "testData", 
      nvar = ncol(rawData),
      varnames = colnames(rawData),
      data = rawData,
      dtype = rep(0, ncol(rawData)))
  
  shinyInput <- list(dtype = 4, 
      xans = 1, 
      yans = 3, 
      nans = 2, 
      CES = 0.1,
      ces.ans = 3,
      cont = FALSE)
  
  savedResults <- fitSingleModel(data = proastData, shinyInput = shinyInput, 
      selectedModel = 25, track = TRUE)
  bindModelResults(savedResults = savedResults)
  #              model parameterNames npar loglik    aic      bmd     bmdl    bmdu
  # (Intercept) Probit                   2 -91.14 186.28 4.003365 3.080731 6.32023
  #             converged
  # (Intercept)      TRUE

  
  # loglik value is ok (-97.54) at the end of f.mm4.cat()
  # but value changes when estimating CI (also with proast 62.10)
  
  # Expected result (Proast 62.10) when not calculating CI
  # npar = 2.00, loglik = -97.54, AIC = 199.08
  # Expected result (Proast 62.10) when calculating CI
  # npar = 2.00, loglik = -91.14, AIC = ??, BMDL =	3.0807, BMDU = 6.3202, BMD = 4.0034 
  
}