% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_RRBS.R
\name{import_RRBS}
\alias{import_RRBS}
\title{Import RRBS files for formatting}
\usage{
import_RRBS(ctrl, trts, dir = getwd())
}
\arguments{
\item{ctrl}{Character(s) denoting control treatment in file names}

\item{trts}{List of character(s) denoting treatment groups in file names}

\item{dir}{Directory with bismark.cov.gz files}
}
\value{
A list of imported files called 'x', ctrl and trts variables
}
\description{
Import RRBS files for formatting
}
