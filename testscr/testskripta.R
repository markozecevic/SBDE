devtools::use_package("tools", "Imports")
devtools::use_package("VennDiagram", "Imports")
devtools::use_package("grid", "Imports")
devtools::use_package("sevenbridges", "Imports") # adds to description
devtools::use_package("methods", "Imports") # adds to description

devtools::document() # use roxygen to make documentation and add stuff to namespace
devtools::load_all()

your_token <- "29dd0f0d333e473cbb094e8b1014e0d1"
your_project <- "marko_zecevic/sbde-test"
your_platform <- "https://api.sbgenomics.com/v2"

DEA5 <- newDEA("Hisat_and_cuffdiff", your_token, your_project, your_platform, align_wf = "hisat2", de_wf = "cufflinks")

DEA5 <- align(DEA5, should_run = TRUE)
# wait
DEA5@aligned_reads <- list_bam(DEA5)
DEA5@aligned_reads$condition <- c("untreated", "untreated", "treated", "treated", "untreated", "treated", "treated", "untreated")
DEA5 <- analyzeForDE(DEA5, should_run = TRUE)
DEA6 <- newAnalysisAlreadyAligned(DEA5, "Hisat_and_DESeq2", "deseq2")
DEA6 <- analyzeForDE(DEA6, should_run = TRUE)
# wait
DEA3 <- readResults(DEA3)
DEA5 <- readResults(DEA5)
# res_cuff <- DEA1@analysis_results
# res_deseq2 <- DEA2@analysis_results

# http://www.graphviz.org/doc/info/colors.html
plotVenn(analyses = c(DEA1, DEA2, DEA3, DEA4, DEA6), 0.05, lty = "blank", fill = c("skyblue", "pink1", "mediumorchid", "orange", "chartreuse3"))

plotVenn(analyses = c(DEA1, DEA3, DEA5), 0.05, lty = "blank", fill = c("skyblue", "pink1", "mediumorchid"))

plotVenn(analyses = c(DEA1, DEA2, DEA3, DEA4), 0.05, lty = "blank", fill = c("skyblue", "pink1", "mediumorchid", "orange"))


# setwd("~/R/SBDE/testscr")
# save(DEA1, file = "DEA1.RData")
# save(DEA2, file = "DEA2.RData")
# save(DEA3, file = "DEA3.RData")
# save(DEA4, file = "DEA4.RData")
# save(DEA5, file = "DEA5.RData")
# save(DEA6, file = "DEA6.RData")
