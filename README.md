Detecting differentially expressed genes across sample groups is among the major goals in statistical analysis of RNA-seq data. There is a number of tools developed for this purpose and some are already wrapped and made available on [Seven Bridges Cancer Genomics Cloud (CGC)](http://cgc.sbgenomics.com/) platform, which also hosts massive genomic datasets, available immediately to authorized researchers.

**SBDE** package is a collection of API scripts written with the goal of simplifying comparison of differentially expressed gene lists produced on Seven Bridges platforms. It abstracts over the differences in output formats and includes some plotting functionality.

A vignette can be previewed [HERE](https://htmlpreview.github.io/?https://github.com/markozecevic/SBDE/blob/master/inst/doc/SBDE.html).

Before installing **SBDE**, make sure to install the latest version of **sevenbridges** package:

```{r, echo=TRUE}
library(devtools)
install_github("sbg/sevenbridges-r")
install_github("markozecevic/SBDE")
```
