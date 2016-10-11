#' Reading Analysis Results
#'
#' \code{readResults} is used to read the results once the DE testing is complete.
#'
#' @param analysis an object of class \code{DifferentialExpressionAnalysis},
#'   containing the task id of the finished quantification/DE task
#' @return The output is a \code{DifferentialExpressionAnalysis} object,
#'   holding the same information as the input, but with the analysis results
#'   contained in the \code{analysis_results} slot.
#' @examples
#' testDEA <- readResults(testDEA)

#' @export
readResults <- function(analysis) {
  switch(analysis@quantification_DE_wf,
         deseq2={
           return (read_deseq2(analysis))
         },
         cufflinks={
           return (read_cufflinks(analysis))
         }
  )
}

read_cufflinks <- function(analysis) {
  p <- analysis@project
  resFile <- p$task(id = analysis@quantification_DE_task)$file("cuffdiff.zip")
  resFile$download(tempdir())
  pathZIP <- paste0(tempdir(), "/", resFile$name)
  unzip(pathZIP, files = "gene_exp.diff", exdir = tempdir())
  path <- paste0(tempdir(),"/gene_exp.diff")

  # If there are multiple features with the same name in the gtf, cuffdiff may add
  # " (1 of many)" to their names in the gene_exp.diff. This will unable us to
  # read.table() so we need to remove those.
  x <- readLines(path)
  y <- gsub("(1 of many)", "", x, fixed = TRUE)
  # NEEDS TO BE SAVED IN SOME TEMP FILES FOLDER
  new_path <- paste0(tools::file_path_sans_ext(path), "_mod.", tools::file_ext(path))
  cat(y, file=new_path, sep="\n")

  cdResults <- read.table(file=new_path, header=TRUE, stringsAsFactors=FALSE)
  file.remove(pathZIP, path, new_path)
  analysis@analysis_results <- data.frame(row.names = cdResults$gene_id, "q_value" = cdResults$q_value, stringsAsFactors=FALSE)
  return(analysis)
}

read_deseq2 <- function(analysis) {
  p <- analysis@project
  resFile <- p$task(id = analysis@quantification_DE_task)$file(".csv")
  resFile$download(tempdir())
  path <- paste0(tempdir(), "/", resFile$name)

  ds2Results <- read.csv(file=path, header=TRUE, stringsAsFactors=FALSE)
  file.remove(path)
  analysis@analysis_results <- data.frame(row.names = ds2Results$X, "q_value" = ds2Results$padj, stringsAsFactors=FALSE)
  return(analysis)
}
