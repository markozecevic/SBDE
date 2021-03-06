% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotVenn}
\alias{plotVenn}
\title{Plotting a Venn Diagram}
\usage{
plotVenn(analyses, alpha, ...)
}
\arguments{
\item{analyses}{a vector od objects of class \code{DifferentialExpressionAnalysis}}

\item{alpha}{FDR: q values below this threshold are considered significant}

\item{...}{any additional arguments passed to the parameters of a plot function
from the \code{VennDiagram} package}
}
\value{
A Venn Diagram plot.
}
\description{
\code{plotVenn} plots a Venn diagram of the genes identified as differentially
  expressed in a set of analyses and with the FDR threshold passed as arguments.
}
\examples{
plotVenn(c(DEA1, DEA2, DEA3), 0.05, lty = "blank", fill = c("skyblue", "pink1", "mediumorchid"))
}
