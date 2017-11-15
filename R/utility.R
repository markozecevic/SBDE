fastq_table <- function(project) {
  p <- project
  # cat("creating file index")
  aux <- capture.output(files <- p$file(complete = TRUE))
  # # # this used to be necessary before 'complete = TRUE' was introduced
  # files_aux <- p$file()
  # files <- files_aux
  # os <- 100 # this may change if the limit changes, increments too
  # while (length(files_aux) == 100) {
  #   files_aux <- p$file(offset = os)
  #   files <- FileList(c(as.list(files),as.list(files_aux)))
  #   os <- os + 100
  # }
  if (length(files) == 0) stop("Your project is missing a FASTQ reference.")
  fastq_ind <- which(tools::file_ext(sapply(files, function (x) x$name)) == "fastq")
  fastq_name <- sapply(fastq_ind, function (x) files[[x]]$name)
  meta <- data.frame(t(sapply(fastq_ind, function (x) files[[x]]$meta()[c("sample_id", "paired_end")])))
  rownames(meta) <- fastq_name
  return(meta[order(rownames(meta)), ])
}

#' @export
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
  return(data.frame(bam_name = sapply(children, function(x) x$outputs$sorted_bam$name)))
}

list_bamTophat2 <- function(analysis) {
  p <- analysis@project
  children <- p$task(parent = analysis@alignment_task, status = "completed", detail = TRUE)
  return(data.frame(bam_name = sapply(children, function(x) x$outputs$output_bam_file$name)))
}

list_bamHisat2 <- function(analysis) {
  p <- analysis@project
  children <- p$task(parent = analysis@alignment_task, status = "completed", detail = TRUE)
  return(data.frame(bam_name = sapply(children, function(x) x$outputs$sorted_bam$name)))
}

find_files <- function(project, ext) {
  p <- project
  if (ext %in% c("fa", "fasta")) {
    e <- c("fa", "fasta")
  } else {
    e <- ext
  }
  # cat("creating file index")
  aux <- capture.output(files <- p$file(complete = TRUE))
  # # # this used to be necessary before 'complete = TRUE' was introduced
  # files_aux <- p$file()
  # files <- files_aux
  # os <- 100 # this may change if the limit changes, increments too
  # while (length(files_aux) == 100) {
  #   files_aux <- p$file(offset = os)
  #   files <- FileList(c(as.list(files),as.list(files_aux)))
  #   os <- os + 100
  # }
  if (length(files) == 0) stop(paste0("Your project is missing the needed ", toupper(ext), " files."))
  ext_ind <- which(tools::file_ext(sapply(files, function (x) x$name)) %in% e)
  ext_name <- sapply(ext_ind, function (x) files[[x]]$name)
  return(ext_name)
}
