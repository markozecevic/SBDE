#' SBDE: Differential Expression Analysis on Seven Bridges Platforms
#'
#' A set of API scripts written with the goal of simplifying comparison of
#' differentially expressed gene lists produced on Seven Bridges platforms.
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
