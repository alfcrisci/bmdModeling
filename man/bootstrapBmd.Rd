% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelAveraging.R
\name{bootstrapBmd}
\alias{bootstrapBmd}
\title{Estimate lower and upper bound for model-averaged bmd using parametric bootstrap}
\usage{
bootstrapBmd(proastData, weights, modelResults, shinyInput,
  naiveApproach = FALSE, aicNull = NA, nBootstraps = 200, seed = 1,
  showProgress = FALSE)
}
\arguments{
\item{proastData}{list, data in proast format as returned by f.scan()}

\item{weights}{numeric vector, estimated weights as returned by calculateWeights()}

\item{modelResults}{list, with results for each model, same length as weights. 
For each model a list with at least npar, loglik,  model.ans, regr.par, CES and ces.ans; 
these are by default included in result from f.proast(). Eventually 
contains also fct1 and fct2 if factors are included for the model parameters.}

\item{shinyInput}{list with necessary parameters used for fitted models in 
modelResults}

\item{naiveApproach}{boolean, TRUE if the model-averaged BMD is estimated as the 
weighted average of bmd values, FALSE if the model-averaged BMD is estimated
based on weighted average of response values; default value is FALSE}

\item{aicNull}{numeric, aic value for null model as criterion for accepting 
bootstrap data, if NA all bootstrap data are accepted; default value is NA}

\item{nBootstraps}{integer, the number of bootstrap data sets to generate;
default value is 200}

\item{seed}{integer, allows reproducing results; default value is 1}

\item{showProgress}{boolean, whether progress bar should be shown in shiny
application; important: only use this option when function is called from 
within shiny application; default value is FALSE}
}
\value{
list with modelResults and bootstrapBmd. The modelResults contain 
modelResults for each bootstrap data set; bootstrapBmd is data frame with
all estimated bmd values per group
}
\description{
Estimate lower and upper bound for model-averaged bmd using parametric bootstrap
}
