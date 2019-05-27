
library(bmdModeling)
library(testthat)

dataDir <- system.file("extdata", package = "bmdModeling")

# Compare results with original f.proast()
#source("/home/mvarewyck/git/bmd/exploringProast/fFunctionsProast.R")
#load("/home/mvarewyck/git/bmd/proast61.3/data/das5.rda")
#f.proast(das5)
#metry <- f.scan("/home/mvarewyck/git/bmd/proastUI/inst/extdata/methyleug.txt")
#track <- TRUE
#track2 <- TRUE

track <- FALSE
doPrint <- FALSE 


## Continuous data (dtype == 1) ##

context("Continuous data (das1)")

# load data
load(file.path(dataDir, "das1.rda"))
# Remove missing response values
das1$data <- das1$data[!is.na(das1$data$BW), ]

# shinyInput
shinyInput <- list(
    #continuous data
    dtype = 1, quick.ans = 1, 
    # fit model
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
    sf.x = 100
)


# A. Define constraints
testConstraints(shinyInput = shinyInput, data = das1,
    modelAns = c(1, 11, 13, 15, 23, 25), 
    doPrint = doPrint, track = track)

# B. Define starting values (similar as default)
testStartingValues(shinyInput = shinyInput, data = das1,
    modelAns = c(1, 11, 13, 15, 23, 25), 
    startValues = NULL, results = NULL,
    doPrint = doPrint, track = track)


# C. Include covariates

# wrapper for dtype == 1
testCovariatesDtype1 <- function(...)
  testCovariates(shinyInput = shinyInput, 
      data = das1, modelAns = c(1, 11, 13, 15, 23, 25), 
      doPrint = doPrint, track = track, ...)

context("Covariate in a and b")
results <- list(
    loglik = c(93.31, 155.87, 152.65, 153.61, 152.72, 153.6),
    CED.origScale = matrix(
        c(2447.0, 3870.6, 2468.4, 3835.0,
            1392.1, 1212.2, 1359.4, 1212.1), 
        nrow = 2, byrow = TRUE),
    conf.int = matrix(
        c(1527.90, 2557.80, 1582.90, 1919.90,
            654.14,  812.01,  659.62,  812.81,
            4112.50, 4487.20, 3888.50, 4440.00,
            3524.10, 2089.30, 3810.40, 2101.50), 
        nrow = 4, byrow = TRUE)
)

testCovariatesDtype1(
    shinyInputCovariates = list( 
        fct1.no = 5,
        fct2.no = 5,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 0
    ), results = results
)

# Different covariate in a and b
testCovariatesDtype1(
    shinyInputCovariates = list( 
        fct1.no = 5,
        fct2.no = 6,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 0
    ), results = NULL
)


context("Covariate in var and c")
results <- list(
    loglik = c(6.92, 157.55,  17.79,  18.06,  17.75,  18.01),
    CED.origScale =  c(972.26,  981.91,  996.39, 1007.00),
    conf.int =  matrix(c(
            391.77,  403.75 , 419.57 , 431.11,
            2336.20, 2306.20 ,2328.00, 2309.30), 
        nrow = 2, byrow = TRUE)
)

testCovariatesDtype1(
    shinyInputCovariates = list( 
        fct1.no = 0,
        fct2.no = 0,
        fct3.no = 5,
        fct4.no = 6,
        fct5.no = 0
    ), results = results
)


context("Covariate in a and d")
results <- list(
    loglik =  c(93.31, 155.87, 142.96, 142.96, 142.96, 142.96),
    # why matrix for this model?
    CED.origScale = c(2520.2, 2520.2, 2497.6, 2497.6),
    conf.int = matrix(c(1258.6, 1258.6, 1274.6, 1274.6,
            4107.4, 4107.5, 4181.3, 4181.2), 
        nrow = 2, byrow = TRUE)
)
testCovariatesDtype1(
    shinyInputCovariates = list( 
        fct1.no = 5,
        fct2.no = 0,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 6
    ), results = results
)


if (FALSE) {
  
  # Interactions between factors: takes long time to fit
  testCovariatesDtype1(
      shinyInputCovariates = list( 
          fct1.no = c(5, 6),
          fct2.no = c(5, 6),
          fct3.no = 0,
          fct4.no = c(5, 6),
          fct5.no = c(5, 6)
      ), results = NULL
  )
  
}




## Quantal data (dtype == 4) ##

context("Quantal data (methyleug)")

metry <- f.scan(file.path(dataDir, "methyleug.txt"))

shinyInput <- list(
    # quantal
    dtype = 4, cont = FALSE,
    quick.ans = 1, #always 1: single model
    # xans = independent variable = dose
    xans = 1, 
    # response: number of responding animals: forest.kk, liver.kk
    yans = 2, # also 3
    # associated sample size
    nans = 4, # sample.size (or nans?)
    # type of benchmark response (1 = ED50)
    ces.ans = 1,
    # need for f.mm6.cat
    CES = 0.05,
    # confidence level for CED
    conf.lev = 0.9
)

# A. Define constraints
testConstraints(shinyInput, data = metry,
    modelAns = c(1, 14, 16, 18, 19, 21, 24, 25, 26), 
    doPrint = doPrint, track = track)

# B. Define starting values
startValues <- list(
    "1" = c(0.5),
    "14" = c(0.1, 0.2, 0.8, 0.9),
    "16" = c(0.1, 3, 2),
    "18" = c(0.1, 3, 2),
    "19" = c(0.1, 3, 2),
    "21" = c(0.1, 3, 2),
    "24" = c(0.1, 3, 2),
    "25" = c(-0.3, 5),
    "26" = c(-0.3, 5)
)
results <- list(
    par.start =  list(0.46,
        c(0.01, 0.14, 0.70, 0.99),
        c(0.0087307, 7.6881, 2.307),
        c(0.010956, 6.7232, 2.4838),
        c(0.008928, 7.414100, 1.495400 ),
        c(0.010429,    6.675300,    1.433900),
        c(0.0095945, 7.2174000, 2.1732000),
        c(-1.655, 8.474),
        c(-3.1294,   8.0979)
    ),
    loglik = c(-275.98, -112.78, -115.06, -113.59, -113.81,
        -113.01, -113.05, -122.32, -117.75),
    CED.matr = c(7.6881, 6.7232, 7.4141, 6.6753, 
        7.2174, 8.474, 8.0979),
    conf.int =  matrix(
        c(6.661, 8.716, 5.943, 7.585, 6.498, 8.395, 5.939,
         7.49, 6.398, 8.111, 0.01225, Inf, 7.385, 8.871),
        nrow = 2, byrow = FALSE
    )
)
testStartingValues(shinyInput, data = metry,
    modelAns = c(1, 14, 16, 18, 19, 21, 24, 25, 26), 
    startValues = startValues, results = results,
    doPrint, track)


# C. Include covariates

# wrapper for dtype == 4
testCovariatesDtype4 <- function(...)
  testCovariates(shinyInput = shinyInput, 
      data = metry, modelAns = c(1, 14, 16, 18, 19, 21, 24, 25, 26), 
      doPrint = doPrint, track = track, ...)

context("Covariate in a and b")
results <- list(
    loglik =  c(-275.25, -108.80, -110.10, -110.01, 
        -109.51, -109.52, -109.08, -122.32, -113.21),
    CED.matr =  matrix(c(6.2791, 5.6825, 6.1538, 5.6875, 
            6.0339, 8.474, 7.0469, 9.4084, 8.0824, 8.9837,
            7.9354, 8.5904, 11.818, 9.4152), 
        nrow = 2, byrow = TRUE),
    conf.int = matrix(c(5.3372, 7.9575, 7.3571, 11.067, 
    4.814, 6.7859, 6.684, 9.6055, 5.2101, 7.5578, 7.2396, 
    10.622, 4.8465, 6.7069, 6.6492, 9.472, 5.1319, 7.2969, 
    7.0778, 10.117, 0.012249, 0, Inf, Inf, 6.1703, 8.1689, 
    7.9533, 10.995), nrow = 4, byrow = FALSE)
)
testCovariatesDtype4(
    shinyInputCovariates = list( 
        fct1.no = 5,
        fct2.no = 5,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 0
    ), 
    results = results
)

context("Covariate in var and c")
results <- list(
    loglik = c(-275.98, -112.78, -115.06, -113.59, -113.81, -113.01, -113.05, -122.32, -117.75),
    CED.matr = c(7.6881, 6.7232, 7.4141, 6.6753, 7.2174, 8.474, 8.0979) , #x$CED.matr[1]
    conf.int = matrix(c(6.6612, 8.7164, 5.943, 7.5847, 6.498, 8.3948, 
    5.9392, 7.4899, 6.3977, 8.111, 0.012249, Inf, 7.385, 8.8707), 
        nrow = 2, byrow = FALSE)
)

testCovariatesDtype4(
    shinyInputCovariates = list( 
        fct1.no = 0,
        fct2.no = 0,
        fct3.no = 5,
        fct4.no = 5,
        fct5.no = 0,
        fct3.ref = 1
    ), results = results
)
# In proast61.3, covariate for fct4.no is not taken into account for model 18,
# see f.change.settings() change[17]

# Note: no parameter d in non-continuous models.

context("Covariate in a, b and c")
testCovariatesDtype4(
    shinyInputCovariates = list( 
        fct1.no = 5,
        fct2.no = 5,
        fct3.no = 0,
        fct4.no = 5,
        fct5.no = 0
    ), results = NULL
)


## Continuous summary data (dtype == 10)##

context("Continuous summary data (atra)")

atra <- f.scan(file.path(dataDir, "atra.txt"))

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
    cont = TRUE,
    # confidence level for CED
    conf.lev = 0.9
)

# A. Define constraints
testConstraints(shinyInput, data = atra,
    modelAns = c(1, 11, 13, 15, 23, 25), 
    doPrint = doPrint, track = track)

# B. Define starting values
startValues <- list(
    "1" = c(0.02, 450),
    "11" = c(0.02, rep(450, 5)),
    "13" = c(0.02, 450, 5, 2),
    "15" = c(0.02, 450, 5, 0.5, 2),
    "23" = c(0.02, 450, 5, 2),
    "25" = c(0.02, 450, 5, 0.5, 2))
results <- list(
    par.start = list(c(0.068009, 437.79),
        c(0.05495, 485.11, 459.57, 468.54, 388.52, 352.22),
        c(0.055379, 480.02, 3.5454, 0.6884),
        c(0.055286, 473.73, 8.349, 0.73293, 1.5355),
        c(0.055367, 479.1, 3.9843, 0.77108),
        c(0.055288, 473.5, 8.2516, 0.68783, 1.7636)),
    loglik = c(-11.16, 4.73, 4.15, 4.27, 4.16, 4.27),
    CED = c(3.5454, 8.349, 3.9843, 8.2516),
    conf.int = matrix(
        c(0.1746, 14.67, 0.2077, 26.75, 0.2221, 15.06, 0.2615, 22.38), 
        nrow = 2, byrow = FALSE)
)

testStartingValues(shinyInput, data = atra,
    modelAns = c(1, 11, 13, 15, 23, 25), 
    startValues = startValues, results = results,
    doPrint, track)

# C. Include covariates

# load data with covariates
load(file.path(dataDir, "das10.rda"))

shinyInput <- list(dtype = 10, quick.ans = 1, 
    # continous model
    model.ans = 1, 
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
    # confidence level for CED
    conf.lev = 0.9
)

# wrapper for dtype == 10
testCovariatesDtype10 <- function(...)
  testCovariates(shinyInput = shinyInput, 
      data = das10, modelAns = c(1, 11, 13, 15, 23, 25), 
      doPrint = doPrint, track = track, ...)

context("Covariate in a and b")
testCovariatesDtype10(
    shinyInputCovariates = list( 
        fct1.no = 3,
        fct2.no = 3,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 0
    ), results = NULL
)

context("Covariate in var and c")
testCovariatesDtype10(
    shinyInputCovariates = list( 
        fct1.no = 0,
        fct2.no = 0,
        fct3.no = 3,
        fct4.no = 3,
        fct5.no = 0
    ), results = NULL
)

context("Covariate in a and d")
testCovariatesDtype10(
    shinyInputCovariates = list( 
        fct1.no = 3,
        fct2.no = 0,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 3
    ), results = NULL
)

context("Covariate in a and c")
testCovariatesDtype10(
    shinyInputCovariates = list( 
        fct1.no = 3,
        fct2.no = 0,
        fct3.no = 0,
        fct4.no = 3,
        fct5.no = 0
    ), results = NULL
)


context("Covariate in a, b, c and d")
testCovariatesDtype10(
    shinyInputCovariates = list( 
        fct1.no = 3,
        fct2.no = 3,
        fct3.no = 0,
        fct4.no = 3,
        fct5.no = 3
    ), results = NULL
)

# Note: For the data types below we just checked whether code works, 
# so results are not compared (yet) with those obtained in proast61.3

## Binary data ##

context("Binary data (das2)")

load(file.path(dataDir, "das2.rda"))

shinyInput <- list(
    # binary
    dtype = 2, cont = FALSE,
    quick.ans = 1, #always 1: single model
    # xans = independent variable = dose
    xans = 1, 
    # response: number of responding animals: forest.kk, liver.kk
    yans = 2, 
    # give scaling factor for dose >
    sf.x = 100, 
    # need for f.mm6.cat
    CES = 0.05
)


# A. Define constraints
testConstraints(shinyInput = shinyInput, data = das2,
    modelAns = c(1, 14, 16, 18, 19, 21, 24, 25, 26), 
    doPrint = doPrint, track = track)

# B. Define starting values (similar as default)
testStartingValues(shinyInput = shinyInput, data = das2,
    modelAns = c(1, 14, 16, 18, 19, 21, 24, 25, 26), 
    startValues = NULL, results = NULL,
    doPrint = doPrint, track = track)


# C. Include covariates

# wrapper for dtype == 2
testCovariatesDtype2 <- function(...)
  testCovariates(shinyInput = shinyInput, 
      data = das2, modelAns = c(1, 14, 16, 18, 19, 21, 24, 25, 26), 
      doPrint = doPrint, track = track, ...)

context("Covariate in a and b")
testCovariatesDtype2(
    shinyInputCovariates = list( 
        fct1.no = 3,
        fct2.no = 3,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 0
    ), results = NULL
)

context("Covariate in var and c")
testCovariatesDtype2(
    shinyInputCovariates = list( 
        fct1.no = 0,
        fct2.no = 0,
        fct3.no = 3,
        fct4.no = 3,
        fct5.no = 0, 
        fct3.ref = 1
    ), results = NULL
)

context("Covariate in a and c")
testCovariatesDtype2(
    shinyInputCovariates = list( 
        fct1.no = 3,
        fct2.no = 0,
        fct3.no = 0,
        fct4.no = 3,
        fct5.no = 0
    ), results = NULL
)
# No parameter d in non-continuous models.



## Ordinal data ##

context("Ordinal data (das3)")

# load data
load(file.path(dataDir, "das3.rda"))

das3$data <- f.remove.NAs(variableIndices = 
        # xans, Vyans, sans, covar.no, nans
        as.numeric(c(2,	3, NULL, 1, NULL)), 
    originalData = das3$data, track = track)

shinyInput <- list(
    dtype = 3, quick.ans = 1, 
    # dose
    xans = 2, 
    # mean response
    yans = 3, 
    # extra parameters
    CES.cat = 2, # severity category associated with CED
    # critical effect size
    CES = 0.05, 
    # is model continous?
    cont = FALSE,
    # confidence level for CI
    conf.lev = 0.9)

# A. Define constraints
testConstraints(shinyInput = shinyInput, data = das3,
    modelAns = c(1, 14, 13, 15, 23, 25), 
    doPrint = doPrint, track = track)

# B. Define starting values (similar as default)
testStartingValues(shinyInput = shinyInput, data = das3,
    modelAns = c(1, 14, 13, 15, 23, 25), 
    startValues = NULL, results = NULL,
    doPrint = doPrint, track = track)


# C. Include covariates

# wrapper for dtype == 3
testCovariatesDtype3 <- function(...)
  testCovariates(shinyInput = shinyInput, 
      data = das3, modelAns = c(1, 14, 13, 15, 23, 25), 
      doPrint = doPrint, track = track, ...)

context("Covariate in a and b")
testCovariatesDtype3(
    shinyInputCovariates = list( 
        fct1.no = 1,
        fct2.no = 1,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 0
    ), results = NULL
)

## For ordinal data: in Proast61.3 can only use covariate for parameter a and b


## Continuous clustered data ##

context("Continuous clustered data (das5)")

load(file.path(dataDir, "das5.rda"))

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
    sf.x = 1000,
    # number of bootstrap runs
    nruns <- 10
)

# A. Define constraints
testConstraints(shinyInput = shinyInput, data = das5,
    modelAns = c(1, 11, 13, 15, 23, 25), 
    doPrint = doPrint, track = track)

# B. Define starting values (similar as default)
testStartingValues(shinyInput = shinyInput, data = das5,
    modelAns = c(1, 11, 13, 15, 23, 25), 
    startValues = NULL, results = NULL,
    doPrint = doPrint, track = track)


# C. Include covariates

# wrapper for dtype == 5
testCovariatesDtype5 <- function(...)
  testCovariates(shinyInput = shinyInput, 
      data = das5, modelAns = c(1, 11, 13, 15, 23, 25), 
      doPrint = doPrint, track = track, ...)

context("Covariate in a and b")
testCovariatesDtype5(
    shinyInputCovariates = list( 
        fct1.no = 2,
        fct2.no = 2,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 0
    ), results = NULL
)

context("Covariate in a and c")
# Proast: Litter effects not implemented for different within group variances -> set fct3.no = 0
testCovariatesDtype5(
    shinyInputCovariates = list( 
        fct1.no = 2,
        fct2.no = 0,
        fct3.no = 0,
        fct4.no = 3,
        fct5.no = 0
    ), results = NULL
)

context("Covariate in a and d")
testCovariatesDtype5(
    shinyInputCovariates = list( 
        fct1.no = 2,
        fct2.no = 0,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 3
    ), results = NULL
)


## Quantal clustered data ##

context("Quantal clustered data (das6)")

load(file.path(dataDir, "das6.rda"))

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

# A. Define constraints
testConstraints(shinyInput, data = das6,
    modelAns = c(1, 14, 16, 18, 19, 21, 24, 25, 26), 
    doPrint = doPrint, track = track)

# B. Define starting values
testStartingValues(shinyInput, data = das6,
    modelAns = c(1, 14, 16, 18, 19, 21, 24, 25, 26), 
    startValues = NULL, results = NULL,
    doPrint = doPrint, track = track)


# C. Include covariates

# wrapper for dtype == 6
testCovariatesDtype6 <- function(...)
  testCovariates(shinyInput = shinyInput, 
      data = das6, modelAns = c(1, 14, 16, 18, 19, 21, 24, 25, 26), 
      doPrint = doPrint, track = track, ...)

context("Covariate in a and b")
testCovariatesDtype6(
    shinyInputCovariates = list( 
        fct1.no = 4,
        fct2.no = 4,
        fct3.no = 0,
        fct4.no = 0,
        fct5.no = 0
    ), results = NULL
)

context("Covariate in var and c")
testCovariatesDtype6(
    shinyInputCovariates = list( 
        fct1.no = 0,
        fct2.no = 0,
        fct3.no = 4,
        fct4.no = 4,
        fct5.no = 0,
        fct3.ref = 1
    ), results = NULL
)
# In proast61.3, covariate for fct4.no is not taken into account for model 18,
# see f.change.settings() change[17]

# Note: no parameter d in non-continuous models.