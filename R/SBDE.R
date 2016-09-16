#' SBDE: Differential Expression Analysis on Seven Bridges Platforms
#'
#' A set of API scripts written with a goal of simplifying differential expression
#' analysis on Seven Bridges platforms. The DifferentialExpressionAnalysis container
#' is used to keep track of data, chosen alignment and quantification/DE testing
#' workflows and their corresponding tasks. SBDE also extends the functionality of
#' VennDiagram package to allow for comparison of up to five different workflow
#' choices.
#'
#' @section Vignette:
#' LINK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#'
#' @docType package
#' @name SBDE
#'
#' @import sevenbridges

NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Differential Expression Analysis on Seven Bridges Platforms")
}
