% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_h5.R
\name{readH5SCE}
\alias{readH5SCE}
\alias{readH5SCME}
\title{Read a SingleCellExperiment from HDF5 file}
\usage{
readH5SCE(file, name = "/", verbose = TRUE, ...)

readH5SCME(file, name = "/", verbose = TRUE, ...)
}
\arguments{
\item{file}{An HDF5 file.}

\item{name}{An HDF5 group where the \code{SingleCellExperiment} is stored.}

\item{verbose}{Display progress.}

\item{...}{Arguments passed to other methods}
}
\value{
\itemize{
\item \code{readH5SCE} returns an object inheriting from
\code{\link[SingleCellExperiment]{SingleCellExperiment}}.
\item \code{readH5SCME} returns an object inheriting from
\code{\link{SingleCellMultiExperiment}}
}
The concrete class is automatically deduced.
}
\description{
Read a SingleCellExperiment from HDF5 file
}
