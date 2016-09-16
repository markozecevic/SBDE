
newAnalysisAlreadyAligned <- function(analysis, title, quantwf) {
  NEWanalysis <- analysis
  NEWanalysis@title <- title
  NEWanalysis@quantification_DE_wf <- quantwf
  NEWanalysis@quantification_DE_task  <- NA_character_
  return(NEWanalysis)
}

fastq_table <- function(project) {
  p <- project
  # cat("creating file index")
  # files <- p$file(complete = TRUE)
  # # this used to be necessary before 'complete = TRUE' was introduced
  files_aux <- p$file()
  files <- files_aux
  os <- 100 # this may change if the limit changes, increments too
  while (length(files_aux) == 100) {
    files_aux <- p$file(offset = os)
    files <- FileList(c(as.list(files),as.list(files_aux)))
    os <- os + 100
  }
  fastq_ind <- which(tools::file_ext(sapply(files, function (x) x$name)) == "fastq")
  fastq_name <- sapply(fastq_ind, function (x) files[[x]]$name)
  meta <- data.frame(t(sapply(fastq_ind, function (x) files[[x]]$meta()[c("sample_id", "paired_end")])))
  rownames(meta) <- fastq_name
  return(meta[order(rownames(meta)), ])
}

list_bam <- function(analysis) {
  switch(analysis@alignment_wf,
         star={
           return (list_bamStar(analysis))
         },
         tophat2={
           return (list_bamTophat2(analysis))
         },
         hisat2={
           return (list_bamHisat2(analysis))
         }
  )
}

list_bamStar <- function(analysis) {
  p <- analysis@project
  children <- p$task(parent = analysis@alignment_task, status = "completed", detail = TRUE)
  return(data.frame(bam_names = sapply(children, function(x) x$outputs$sorted_bam$name)))
}

list_bamTophat2 <- function(analysis) {
  p <- analysis@project
  children <- p$task(parent = analysis@alignment_task, status = "completed", detail = TRUE)
  return(data.frame(bam_names = sapply(children, function(x) x$outputs$output_bam_file$name)))
}

list_bamHisat2 <- function(analysis) {
  p <- analysis@project
  children <- p$task(parent = analysis@alignment_task, status = "completed", detail = TRUE)
  return(data.frame(bam_names = sapply(children, function(x) x$outputs$sorted_bam$name)))
}

find_files <- function(project, ext) {
  p <- project
  if (ext %in% c("fa", "fasta")) {
    e <- c("fa", "fasta")
  } else {
    e <- ext
  }
  # cat("creating file index")
  # files <- p$file(complete = TRUE)
  # # this used to be necessary before 'complete = TRUE' was introduced
  files_aux <- p$file()
  files <- files_aux
  os <- 100 # this may change if the limit changes, increments too
  while (length(files_aux) == 100) {
    files_aux <- p$file(offset = os)
    files <- FileList(c(as.list(files),as.list(files_aux)))
    os <- os + 100
  }
  ext_ind <- which(tools::file_ext(sapply(files, function (x) x$name)) %in% e)
  ext_name <- sapply(ext_ind, function (x) files[[x]]$name)
  return(ext_name)
}

status <- function(analysis) {
  p <- analysis@project
  if (is.na(analysis@alignment_task)&&is.na(analysis@quantification_DE_task)) {
    cat("No tasks have been run.\n")
  } else if (is.na(analysis@quantification_DE_task)) {
    tsk <- p$task(id = analysis@alignment_task)
    s <- tsk$execution_status
    if (s$message=="COMPLETED") {
      cat("Alignment is completed. Proceed with the quantification step.\n")
    } else {
      if (is.null(s$completed)) {
        msg <- paste0("Alignment in progress: none of the tasks are completed and ",
                      s$running," still running.\n")
      } else msg <- paste0("Alignment in progress: ", s$completed,
                           " of the tasks completed and ", s$running," still running.\n")
      cat(msg)
    }
  }
}
