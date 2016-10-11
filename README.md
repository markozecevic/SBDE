Detecting differentially expressed genes across sample groups is among the major goals in statistical analysis of RNA-seq data. There is a number of tools developed for this purpose and some are already wrapped and made available on Seven Bridges Cancer Genomics Cloud (CGC) platform, which also hosts massive genomic datasets, available immediately to authorized researchers.

**SBDE** package is a collection of API functions bundled together to simplify differential expression analysis on Seven Bridges platforms. It abstracts over the differences in output formats and setting up and running different workflows, allowing an inexperienced platform user to perform his analysis while getting familiar with the graphical interface. Some plotting functionality is also included for visual inspection and comparison of results.

A vignette can be previewed [HERE](https://htmlpreview.github.io/?https://github.com/markozecevic/SBDE/blob/master/inst/doc/SBDE.html).

Before installing **SBDE**, make sure to install the latest version of **sevenbridges** package:

```{r, echo=TRUE}
library(devtools)
install_github("sbg/sevenbridges-r")
install_github("markozecevic/SBDE")
```
