#' S4 Class Representing a DE Analysis
#'
#' Class \code{DifferentialExpressionAnalysis} is a container for holding all the relevant
#'   info regarding a differential expression analysis on one of Seven Bridges platforms.
#'
#' @name DifferentialExpressionAnalysis-class
#' @exportClass DifferentialExpressionAnalysis
#' @slot title analysis title
#' @slot project project id
#' @slot reference reference FASTA
#' @slot annotation reference annotation
#' @slot sample_reads a data frame indicating sample FASTQs to be used in the analysis
#' @slot alignment_wf a string representing the chosen alignment workflow
#' @slot alignment_task id of the alignment task
#' @slot aligned_reads a data frame listing the alignment files produced
#' @slot quantification_DE_wf a string representing the chosen quantification/DE workflow
#' @slot quantification_DE_task id of the chosen quantification/DE task
#' @slot analysis_results the results of DE analysis

DifferentialExpressionAnalysis <- setClass(
  # Set the name for the class
  "DifferentialExpressionAnalysis",

  # Define the slots
  representation(
    title = "character",
    project = "Project",
    reference = "character",
    annotation = "character",
    sample_reads = "data.frame",
    alignment_wf = "character",
    alignment_task = "character",
    aligned_reads = "data.frame",
    quantification_DE_wf  = "character",
    quantification_DE_task = "character",
    analysis_results = "data.frame"
  ),

  # NAMESTI DA IDjevi SAMO BUDU NA_character_
  prototype(
    alignment_task = NA_character_,
    quantification_DE_task  = NA_character_
  )
)

# https://www.r-bloggers.com/vectors-of-s4-classes-with-non-trivial-slots/

setMethod("c", signature(x = "DifferentialExpressionAnalysis"), function(x, ...){
  elements = list(x, ...)

  AnalysesList = list()
  for (i in 1:length(elements)){
    AnalysesList[i] = new("DifferentialExpressionAnalysis",
                          title = slot(elements[[i]], "title"),
                          project = slot(elements[[i]], "project"),
                          reference = slot(elements[[i]], "reference"),
                          annotation = slot(elements[[i]], "annotation"),
                          sample_reads = slot(elements[[i]], "sample_reads"),
                          alignment_wf = slot(elements[[i]], "alignment_wf"),
                          alignment_task = slot(elements[[i]], "alignment_task"),
                          aligned_reads = slot(elements[[i]], "aligned_reads"),
                          quantification_DE_wf = slot(elements[[i]], "quantification_DE_wf"),
                          quantification_DE_task = slot(elements[[i]], "quantification_DE_task"),
                          analysis_results = slot(elements[[i]], "analysis_results"))
  }

  class(AnalysesList) = "DifferentialExpressionAnalysis"

  AnalysesList

})

# # will do someday
# setMethod("show", "DifferentialExpressionAnalysis", function(object){
#   cat("An object of class DifferentialExpressionAnalysis\n")
#   cat(paste0('\nTitle: "',object@title,'"'))
# })

#' New Differential Expression Analysis
#'
#' \code{newDEA} creates a new object of class \code{DifferentialExpressionAnalysis}.
#'
#' @param analysis_title the title of your analysis
#' @param token your Seven Bridges Platform authentication token
#' @param project_id project id
#' @param platform chosen Seven Bridges platform
#' @param align_wf a character string indicating the choice of alignment workflow
#'   (choose between \emph{star}, \emph{tophat2} and \emph{hisat2})
#' @param de_wf a character string indicating the choice of quantification/DE workflow
#'   (choose between \emph{cufflinks} and \emph{deseq2})
#' @return The output is a new \code{DifferentialExpressionAnalysis} object,
#'   with \emph{title}, \emph{project}, \emph{alignment_wf} and \emph{quantification_DE_wf}
#'   slots populated. Also, any FASTA or FA files found in the project are going to be listed
#'   in the \emph{reference} slot - change the value manually if needed so it reflects your
#'   choice of reference. Same goes for the \emph{annotation} and \emph{sample_reads} slots.
#'   Ideally, you will create a new project for the analysis and import just the desired
#'   FASTQs, a single FASTA and a single GTF file.
#' @examples
#' testDEA <- newDEA("Hisat_and_deseq2", your_token, your_project, your_platform, "hisat2", "deseq2")

#' @export
newDEA <- function(analysis_title, token, project_id, platform, align_wf, de_wf) {
  if (!(align_wf %in% c("star", "tophat2", "hisat2"))) stop("Please specify a valid aligner!")
  if (!(de_wf %in% c("cufflinks", "deseq2"))) stop("Please specify a valid quantification/differential expression toolkit!")
  a <- Auth(token = token, url = platform)
  p <- a$project(id = project_id)
  object <- DifferentialExpressionAnalysis(title = analysis_title, project = p, reference = find_files(p, "fasta"),
                                           annotation = find_files(p, "gtf"),
                                           sample_reads = fastq_table(p),
                                           alignment_wf = align_wf,
                                           quantification_DE_wf = de_wf)
  return(object)
}
