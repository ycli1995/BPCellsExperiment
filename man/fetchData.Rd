% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/SingleCellExperiment.R
\name{fetchData}
\alias{fetchData}
\alias{variableFeatures}
\alias{variableFeatures<-}
\alias{fetchData,SingleCellExperiment-method}
\alias{variableFeatures,SingleCellExperiment-method}
\alias{variableFeatures<-,SingleCellExperiment,character-method}
\alias{variableFeatures<-,SingleCellExperiment,NULL-method}
\title{Fetch cellular data}
\usage{
fetchData(object, ...)

variableFeatures(object, ...)

variableFeatures(object, ...) <- value

\S4method{fetchData}{SingleCellExperiment}(object, vars = NULL, cells = NULL, assay = NULL, use.Exp = NULL)

\S4method{variableFeatures}{SingleCellExperiment}(object, ...)

\S4method{variableFeatures}{SingleCellExperiment,character}(object, ...) <- value

\S4method{variableFeatures}{SingleCellExperiment,`NULL`}(object, ...) <- value
}
\arguments{
\item{object}{A \code{\link[SingleCellExperiment]{SingleCellExperiment}}}

\item{...}{Arguments to be passed to other methods}

\item{value}{An object of a class specified in the S4 method signature.}

\item{vars}{Vector of all variables to fetch}

\item{cells}{Cells to collect data for (default is all cells)}

\item{assay}{Assay to collect data from}

\item{use.Exp}{Which \code{\link[SingleCellExperiment]{altExp}} to collect
data from. If the selected \code{altExp} doesn't exist, will use the main
Experiment.}
}
\value{
A \code{\link[base]{data.frame}} with cells as rows and \code{vars} data as
columns.
}
\description{
Fetch data for a set of observations (columns) in an object
}
