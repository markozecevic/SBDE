#' Quantification and Differential Expression Analysis
#'
#' \code{analyzeForDE} creates the quantification/DE analysis task.
#'
#' @param analysis an object of class \code{DifferentialExpressionAnalysis},
#'   containing the choice of workflow and alignment files
#' @param should_run a logical value indicating whether the draft task is to
#'   be executed automatically
#' @return The output is a \code{DifferentialExpressionAnalysis} object,
#'   holding the same information as the input, but with the task id of the
#'   quantification/DE task in the \code{quantification_DE_task} slot.
#' @examples
#' testDEA <- analyzeForDE(testDEA, should_run = TRUE)

analyzeForDE <- function(analysis, should_run = FALSE) {
  switch(analysis@quantification_DE_wf,
         deseq2={
           return (analyzeDESeq2(analysis, should_run))
         },
         cufflinks={
           return (analyzeCufflinks(analysis, should_run))
         }
  )
}

analyzeDESeq2 <- function(analysis, should_run) {
  p <- analysis@project
  files <- p$file(complete = TRUE)

  ## CHECK THAT INPUT FILES ARE VALID - ONE REF, ONE GTF ETC.
  appid <- paste0(p$id,"/htseq-deseq2")
  app <- tryCatch({
    p$app(id = appid)
  },
  error = function(cond) {
    cwl_fl <- system.file("cwl", "rna-seq-diffexp-deseq2.json", package = "SBDE")
    p$app_add(short_name = "htseq-deseq2", filename = cwl_fl)
  })

  conditions <- unique(analysis@aligned_reads$condition)

  ### HAVE A COUPLE OF WORKFLOWS FOR 2/3/4 DIFFERENT CONDITIONS

  group1 <- conditions[1]
  group2 <- conditions[2]

  cond1 <- as.character(analysis@aligned_reads$bam_names[analysis@aligned_reads$condition == group1])
  cond2 <- as.character(analysis@aligned_reads$bam_names[analysis@aligned_reads$condition == group2])

  cond1_in <- files[which(sapply(files, function (x) x$name) %in% cond1)]
  cond2_in <- files[which(sapply(files, function (x) x$name) %in% cond2)]

  gtf_in <- files[[which(sapply(files, function (x) x$name) == analysis@annotation)]]

  # # If you don't know your app's input ids:
  # app <- p$app(id = "marko_zecevic/sbde1/htseq-deseq2/")
  # noInputs <- length(app$raw$inputs)
  # for (i in 1:noInputs) {
  #   cat(app$raw$inputs[[i]]$id)
  #   cat("\n")
  # }

  # shoud inherit the aligner name
  taskName = paste0("DESeq2 ",date())

  p$task_add(name = taskName,
             description = "DE analysis", app = appid,
             inputs = list(features = gtf_in,
                           input_files = cond1_in, input_files_1 = cond2_in,
                           group_name = group1, group_name_1 = group2, ref_condition = group1))

  tsk <- p$task(taskName)
  if (should_run) {
    try(tsk$run())
  }
  analysis@quantification_DE_task <- tsk$id
  return(analysis)
}


analyzeCufflinks <- function(analysis, should_run) {
  p <- analysis@project
  files <- p$file(complete = TRUE)

  ## CHECK THAT INPUT FILES ARE VALID - ONE REF, ONE GTF ETC.
  appid <- paste0(p$id,"/rna-seq-differential-expression")
  app <- tryCatch({
    p$app(id = appid)
  },
  error = function(cond) {
    cwl_fl <- system.file("cwl", "rna-seq-diffexp-cufflinks.json", package = "SBDE")
    p$app_add(short_name = "rna-seq-differential-expression", filename = cwl_fl)
  })

  conditions <- unique(analysis@aligned_reads$condition)

  ### HAVE A COUPLE OF WORKFLOWS FOR 2/3/4 DIFFERENT CONDITIONS

  group1 <- conditions[1]
  group2 <- conditions[2]

  cond1 <- as.character(analysis@aligned_reads$bam_names[analysis@aligned_reads$condition == group1])
  cond2 <- as.character(analysis@aligned_reads$bam_names[analysis@aligned_reads$condition == group2])

  cond1_in <- files[which(sapply(files, function (x) x$name) %in% cond1)]
  cond2_in <- files[which(sapply(files, function (x) x$name) %in% cond2)]

  gtf_in <- files[[which(sapply(files, function (x) x$name) == analysis@annotation)]]
  ref_in <- files[[which(sapply(files, function (x) x$name) == analysis@reference)]]

  # # If you don't know your app's input ids:
  # app <- p$app(id = "marko_zecevic/sbde1/rna-seq-differential-expression/")
  # noInputs <- length(app$raw$inputs)
  # for (i in 1:noInputs) {
  #   cat(app$raw$inputs[[i]]$id)
  #   cat("\n")
  # }

  # shoud inherit the aligner name
  taskName = paste0("Cufflinks ",date())

  ### CHANGE THE INPUT IDs IN THE PUBLIC WORKFLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  p$task_add(name = taskName,
             description = "DE analysis", app = appid,
             inputs = list(Annotations = gtf_in,
                           Group_ERR315421 = cond1_in, Group_ERR315335 = cond2_in,
                           group_name = group1, group_name_1 = group2,
                           Reference = ref_in))

  tsk <- p$task(taskName)
  if (should_run) {
    try(tsk$run())
  }
  analysis@quantification_DE_task <- tsk$id
  return(analysis)
}
