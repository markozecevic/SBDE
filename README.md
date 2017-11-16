Detecting differentially expressed genes across sample groups is among the major goals in statistical analysis of RNA-seq data. There is a number of tools developed for this purpose and some are already wrapped and made available on [Seven Bridges Cancer Genomics Cloud (CGC)](http://cgc.sbgenomics.com/) platform, which also hosts massive genomic datasets, available immediately to authorized researchers.

**SBDE** package is a collection of API scripts written with the goal of simplifying comparison of differentially expressed gene lists produced on Seven Bridges platforms. It abstracts over the differences in output formats and includes some plotting functionality.

Before installing **SBDE**, make sure to install the latest version of **sevenbridges** package:

```{r, echo=TRUE}
library(devtools)
install_github("sbg/sevenbridges-r")
install_github("markozecevic/SBDE")
```

To use it, load **sevenbridges** and create an auth object.

```{r, echo=TRUE}
library(sevenbridges)

a <- Auth(token = "8284fcd1810140ee87474b4b559cd21f", url = "https://api.sbgenomics.com/v2/")
```

Finally, make some plots!

```{r. echo=TRUE}
library(SBDE)

de1 <- newDEA(a, "first test", "c1db11ef-d233-40cc-a176-d96dd6d65b28")
de2 <- newDEA(a, "second test", "1fb7d850-af45-47cd-8370-db97e6c84e3f")
de3 <- newDEA(a, "one more test", "9c4e9390-f611-4230-8e8d-1b7324d854b9")
plotVenn(c(de1, de2, de3), 0.05, lty = "blank", fill = c("skyblue", "pink1", "mediumorchid"))
plotUpset(c(de1,de2,de3), 0.05)
clusterDendogram(c(de1, de2, de3), 0.05)
barChart(c(de1, de2, de3), 0.05, fill = c("skyblue", "pink1", "mediumorchid"))
```

