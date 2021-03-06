% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmdFunctions_visualization.R
\name{f.plot.all}
\alias{f.plot.all}
\title{Wrapper function for plotting the results of proast}
\usage{
f.plot.all(ans.all, track = FALSE, plot.type = NA)
}
\arguments{
\item{ans.all}{list, with all results that were obtained during the analysis}

\item{track}{logical, if TRUE (FALSE by default) print the name of the function which is currently being run}

\item{plot.type, }{integer, defines the type of plot; for cont choose 1 (y vs x), 
for dtype 2 or 3 choose one of 3 (y.up vs. x), 5 (cumulative y vs. x); 
default is 1 for cont and 5 else}
}
\value{
list, updated version of ans.all
}
\description{
Wrapper function for plotting the results of proast
}
