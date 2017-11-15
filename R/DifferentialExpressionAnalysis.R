#' S4 Class Representing a DE Analysis
#'
#' Class \code{DifferentialExpressionAnalysis} is a container for holding all the relevant
#'   info regarding a differential expression analysis on one of Seven Bridges platforms.
#'
#' @name DifferentialExpressionAnalysis-class
#' @exportClass DifferentialExpressionAnalysis
#' @slot title analysis title
#' @slot task id of the chosen DEA task
#' @slot results the results of DE analysis

DifferentialExpressionAnalysis <- setClass(
  # Set the name for the class
  "DifferentialExpressionAnalysis",

  # Define the slots
  representation(
    title = "character",
    task = "Task",
    results = "data.frame"
  ),

  prototype(
    title = NA_character_
  )
)

# https://www.r-bloggers.com/vectors-of-s4-classes-with-non-trivial-slots/

setMethod("c", signature(x = "DifferentialExpressionAnalysis"), function(x, ...){
  elements = list(x, ...)

  AnalysesList = list()
  for (i in 1:length(elements)){
    AnalysesList[i] = new("DifferentialExpressionAnalysis",
                          title = slot(elements[[i]], "title"),
                          task = slot(elements[[i]], "task"),
                          results = slot(elements[[i]], "results"))
  }

  class(AnalysesList) = "DifferentialExpressionAnalysis"

  AnalysesList

})

as.list.nssItem=function(from) mapply(function(x) slot(from,x),
                                      slotNames("nssItem"),
                                      SIMPLIFY=FALSE)

# # will do someday
# setMethod("show", "DifferentialExpressionAnalysis", function(object){
#   cat("An object of class DifferentialExpressionAnalysis\n")
#   cat(paste0('\nTitle: "',object@title,'"'))
# })

#' New Differential Expression Analysis
#'
#' \code{newDEA} creates a new object of class \code{DifferentialExpressionAnalysis}.
#'
#' @param auth authentication object
#' @param analysis_title the title of your analysis
#' @param task_id task id
#' @return The output is a new \code{DifferentialExpressionAnalysis} object
#' @examples
#' library(sevenbridges)
#' a <- Auth(token = your_token, url = your_platform)
#' x <- newDEA(a, "Hisat_and_deseq2", "210458f1-0569-4f72-b5d2-0a64ca966c88")

#' @export
newDEA <- function(auth, analysis_title, task_id) {
  t <- auth$task(id = task_id)
  arrid <- unlist(strsplit(t$app, "/"))
  appid <- arrid[length(arrid)-1]
  type <- sapply(c("deseq2", "cufflinks"), grepl, appid)
  # until we add "cufflinks" to the public DE forkflow id:
  if (appid == "rna-seq-differential-expression") type[2] = TRUE
  if (sum(type)!=1) stop("No (or multiple) DE tool names recognized in the app id.")
  switch(names(which(type)),
         deseq2={
           rslt <- read_deseq2(t)
         },
         cufflinks={
           rslt <- read_cufflinks(t)
         })
  object <- DifferentialExpressionAnalysis(title = analysis_title, task = t, results = rslt)
  return(object)
}
