#' Sequence Alignment
#'
#' \code{align} creates the alignment task.
#'
#' @param analysis an object of class \code{DifferentialExpressionAnalysis},
#'   containing the choice of workflow and input data
#' @param should_run a logical value indicating whether the draft task is to
#'   be executed automatically
#' @return The output is a \code{DifferentialExpressionAnalysis} object,
#'   holding the same information as the input, but with the task id of the
#'   alignment task in the \code{alignment_task} slot.
#' @examples
#' testDEA <- align(testDEA, should_run = TRUE)

align <- function(analysis, should_run = FALSE) {
  switch(analysis@alignment_wf,
         star={
           return (alignStar(analysis, should_run))
         },
         tophat2={
           return (alignTophat2(analysis, should_run))
         },
         hisat2={
           return (alignHisat2(analysis, should_run))
         }
  )
}

alignStar <- function(analysis, should_run) {
  p <- analysis@project
  files <- p$file(complete = TRUE)

  ## CHECK THAT INPUT FILES ARE VALID - ONE REF, ONE GTF ETC.
  appid <- paste0(p$id,"/star-alignment")
  app <- tryCatch({
    p$app(id = appid)
    },
    error = function(cond) {
      cwl_fl <- system.file("cwl", "rna-seq-alignment-star.json", package = "SBDE")
      p$app_add(short_name = "star-alignment", filename = cwl_fl)
    })

  fastqs_in <- files[which(sapply(files, function (x) x$name) %in% rownames(analysis@sample_reads))]

  fasta_in <- files[[which(sapply(files, function (x) x$name) == analysis@reference)]]
  gtf_in <- files[[which(sapply(files, function (x) x$name) == analysis@annotation)]]

  # # If you don't know your app's input ids:
  # app <- p$app(id = "marko_zecevic/sbde-test/star-alignment/")
  # noInputs <- length(app$raw$inputs)
  # for (i in 1:noInputs) {
  #   cat(app$raw$inputs[[i]]$id)
  #   cat("\n")
  # }

  # shoud inherit the aligner name
  taskName = paste0("Star-alignment ",date())

  p$task_add(name = taskName,
                 description = "Batch alignment task", app = appid,
                 batch = sevenbridges::batch(input = "fastq", criteria = "metadata.sample_id", type = "CRITERIA"),
                 inputs = list(sjdbGTFfile = gtf_in,
                               fastq = fastqs_in,
                               genomeFastaFiles = fasta_in))

  tsk <- p$task(taskName)
  if (should_run) {
    try(tsk$run())
  }
  analysis@alignment_task <- tsk$id
  return(analysis)
}


alignTophat2 <- function(analysis, should_run) {
  p <- analysis@project
  files <- p$file(complete = TRUE)

  ## CHECK THAT INPUT FILES ARE VALID - ONE REF, ONE GTF ETC.
  appid <- paste0(p$id,"/rna-seq-alignment-tophat")
  app <- tryCatch({
    p$app(id = appid)
  },
  error = function(cond) {
    cwl_fl <- system.file("cwl", "rna-seq-alignment-tophat2.json", package = "SBDE")
    p$app_add(short_name = "rna-seq-alignment-tophat", filename = cwl_fl)
  })

  fastqs_in <- files[which(sapply(files, function (x) x$name) %in% rownames(analysis@sample_reads))]

  fasta_in <- files[[which(sapply(files, function (x) x$name) == analysis@reference)]]
  gtf_in <- files[[which(sapply(files, function (x) x$name) == analysis@annotation)]]

  # # If you don't know your app's input ids:
  # app <- p$app(id = "marko_zecevic/sbde1/rna-seq-alignment-tophat/")
  # noInputs <- length(app$raw$inputs)
  # for (i in 1:noInputs) {
  #   cat(app$raw$inputs[[i]]$id)
  #   cat("\n")
  # }

  # shoud inherit the aligner name
  taskName = paste0("Tophat2-alignment ",date())

  p$task_add(name = taskName,
             description = "Batch alignment task", app = appid,
             batch = sevenbridges::batch(input = "fastq", criteria = "metadata.sample_id", type = "CRITERIA"),
             inputs = list(GTF = gtf_in,
                           fastq = fastqs_in,
                           fasta_reference = fasta_in))

  tsk <- p$task(taskName)
  if (should_run) {
    try(tsk$run())
  }
  analysis@alignment_task <- tsk$id
  return(analysis)
}


alignHisat2 <- function(analysis, should_run) {
  p <- analysis@project
  files <- p$file(complete = TRUE)

  ## CHECK THAT INPUT FILES ARE VALID - ONE REF, ONE GTF ETC.
  appid <- paste0(p$id,"/rna-seq-alignment-hisat")
  app <- tryCatch({
    p$app(id = appid)
  },
  error = function(cond) {
    cwl_fl <- system.file("cwl", "rna-seq-alignment-hisat2.json", package = "SBDE")
    p$app_add(short_name = "rna-seq-alignment-hisat", filename = cwl_fl)
  })

  fastqs_in <- files[which(sapply(files, function (x) x$name) %in% rownames(analysis@sample_reads))]

  fasta_in <- files[[which(sapply(files, function (x) x$name) == analysis@reference)]]
  gtf_in <- files[[which(sapply(files, function (x) x$name) == analysis@annotation)]]

  # # If you don't know your app's input ids:
  # app <- p$app(id = "marko_zecevic/sbde1/rna-seq-alignment-hisat/")
  # noInputs <- length(app$raw$inputs)
  # for (i in 1:noInputs) {
  #   cat(app$raw$inputs[[i]]$id)
  #   cat("\n")
  # }

  # shoud inherit the aligner name
  taskName = paste0("Hisat2-alignment ",date())

  p$task_add(name = taskName,
             description = "Batch alignment task", app = appid,
             batch = sevenbridges::batch(input = "fastq", criteria = "metadata.sample_id", type = "CRITERIA"),
             inputs = list(input_GTF = gtf_in,
                           fastq = fastqs_in,
                           reference_files = fasta_in))

  tsk <- p$task(taskName)
  if (should_run) {
    try(tsk$run())
  }
  analysis@alignment_task <- tsk$id
  return(analysis)
}
