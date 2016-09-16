## ---- echo=FALSE, results='hide', message = FALSE------------------------
devtools::load_all("/Users/marko/R/SBDE")

## ---- echo=TRUE----------------------------------------------------------
your_token <- "29dd0f0d333e473cbb094e8b1014e0d1"
your_project <- "marko_zecevic/sbde-test"
your_platform <- "https://api.sbgenomics.com/v2"

## ---- echo=TRUE----------------------------------------------------------
testDEA <- newDEA("Hisat_and_DESeq2", your_token, your_project, your_platform,
                  align_wf = "hisat2", de_wf = "deseq2")
testDEA

## ---- echo=TRUE----------------------------------------------------------
testDEA@sample_reads

## ---- echo=TRUE, eval = FALSE--------------------------------------------
#  testDEA <- align(testDEA, should_run = TRUE)

## ---- echo=FALSE, results='hide'-----------------------------------------
load("~/R/SBDE/testscr/testDEAaligned.RData")

## ---- echo=TRUE----------------------------------------------------------
testDEA@alignment_task

## ---- echo=TRUE----------------------------------------------------------
testDEA@aligned_reads <- list_bam(testDEA)
testDEA@aligned_reads

## ---- echo=TRUE----------------------------------------------------------
testDEA@aligned_reads$condition <- c("untreated", "untreated", "treated", "treated",
                                     "untreated", "treated", "treated", "untreated")
testDEA@aligned_reads

## ---- echo=TRUE, eval = FALSE--------------------------------------------
#  testDEA <- analyzeForDE(testDEA, should_run = TRUE)

## ---- echo=FALSE, results='hide'-----------------------------------------
load("~/R/SBDE/testscr/testDEAquantified.RData")

## ---- echo=TRUE----------------------------------------------------------
testDEA <- readResults(testDEA)
head(testDEA@analysis_results)

## ---- echo=FALSE, results='hide'-----------------------------------------
load("~/R/SBDE/testscr/DEA1.RData")
load("~/R/SBDE/testscr/DEA2.RData")
load("~/R/SBDE/testscr/DEA3.RData")
load("~/R/SBDE/testscr/DEA4.RData")
load("~/R/SBDE/testscr/DEA5.RData")
load("~/R/SBDE/testscr/DEA6.RData")

## ---- echo=TRUE, warning = FALSE-----------------------------------------
analyses <- c(DEA1, DEA2, DEA3, DEA4, DEA5, DEA6)
sapply(analyses, function (x) x@title)

## ---- echo=TRUE, warning = FALSE-----------------------------------------
barChart(analyses, alpha = 0.05, fill = c("skyblue", "pink1", "mediumorchid", "orange",
                                          "chartreuse3", "firebrick2"))

## ---- echo=TRUE----------------------------------------------------------
clusterDendogram(analyses, alpha = 0.05)

## ---- echo=TRUE----------------------------------------------------------
plotVenn(analyses[-4], 0.05, lty = "blank", fill = c("skyblue", "pink1", "mediumorchid",
                                                     "chartreuse3", "firebrick2"))

## ---- echo=TRUE, message = FALSE-----------------------------------------
plotVenn(analyses[c(-3,-4)], 0.05, lty = "blank", fill = c("skyblue", "pink1",
                                                           "chartreuse3", "firebrick2"))

