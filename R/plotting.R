significant <- function(analyses_titles, d) {
  candidates <- d
  for (i in 1:length(analyses_titles)) {
    candidates <- subset(candidates, candidates[analyses_titles[i]] == T)
  }
  nrow(candidates)
}

#' Plotting a Venn Diagram
#'
#' \code{plotVenn} plots a Venn diagram of the genes identified as differentially
#'   expressed in a set of analyses and with the FDR threshold passed as arguments.
#'
#' @param analyses a vector od objects of class \code{DifferentialExpressionAnalysis}
#' @param alpha FDR: q values below this threshold are considered significant
#' @param ... any additional arguments passed to the parameters of a plot function
#'   from the \code{VennDiagram} package
#' @return A Venn Diagram plot.
#' @examples
#' plotVenn(c(DEA1, DEA2, DEA3), 0.05, lty = "blank", fill = c("skyblue", "pink1", "mediumorchid"))

#' @export
plotVenn <- function(analyses, alpha, ...) {
  coln <- sapply(analyses, function(x) x@title)
  sig <- lapply(analyses, function(x) rownames(x@results)[which(x@results < alpha)])
  names(sig) <- coln
  if (!(sum(sapply(sig, length) > 0) > 1)) stop("At least two non-empty sets are needed for an upset plot. There is ",
                                              sum(sapply(sig, length) > 0), ".")
  n <- length(analyses)
  d <- UpSetR::fromList(sig)

  grid::grid.newpage()
  if (n == 1) {
    out <- VennDiagram::draw.single.venn(significant(coln, d), category = coln, ...)
  }
  if (n == 2) {
    out <- VennDiagram::draw.pairwise.venn(significant(coln[1], d), significant(coln[2], d),
                                           significant(coln[1:2], d), category = coln, ...)
  }
  if (n == 3) {
    out <- VennDiagram::draw.triple.venn(significant(coln[1], d), significant(coln[2], d),
                                         significant(coln[3], d), significant(coln[1:2], d),
                                         significant(coln[2:3], d), significant(coln[c(1, 3)], d),
                                         significant(coln, d), category = coln, ...)
  }
  if (n == 4) {
    out <- VennDiagram::draw.quad.venn(significant(coln[1], d), significant(coln[2], d), significant(coln[3], d),
                          significant(coln[4], d), significant(coln[1:2], d), significant(coln[c(1, 3)], d),
                          significant(coln[c(1, 4)], d), significant(coln[2:3], d),
                          significant(coln[c(2, 4)], d), significant(coln[3:4], d), significant(coln[1:3], d),
                          significant(coln[c(1, 2, 4)], d), significant(coln[c(1, 3, 4)], d),
                          significant(coln[2:4], d), significant(coln, d), category = coln, ...)
  }
  if (n == 5) {
    out <- VennDiagram::draw.quintuple.venn(significant(coln[1], d), significant(coln[2], d), significant(coln[3], d),
                               significant(coln[4], d), significant(coln[5], d), significant(coln[1:2], d),
                               significant(coln[c(1, 3)], d), significant(coln[c(1, 4)], d),
                               significant(coln[c(1, 5)], d), significant(coln[2:3], d),
                               significant(coln[c(2, 4)], d), significant(coln[c(2, 5)], d),
                               significant(coln[3:4], d), significant(coln[c(3, 5)], d),
                               significant(coln[4:5], d), significant(coln[c(1, 2, 3)], d),
                               significant(coln[c(1, 2, 4)], d), significant(coln[c(1, 2, 5)], d),
                               significant(coln[c(1, 3, 4)], d), significant(coln[c(1, 3, 5)], d),
                               significant(coln[c(1, 4, 5)], d), significant(coln[c(2, 3, 4)], d),
                               significant(coln[c(2, 3, 5)], d), significant(coln[c(2, 4, 5)], d),
                               significant(coln[c(3, 4, 5)], d), significant(coln[c(1, 2, 3, 4)], d),
                               significant(coln[c(1, 2, 3, 5)], d), significant(coln[c(1, 2, 4, 5)], d),
                               significant(coln[c(1, 3, 4, 5)], d), significant(coln[c(2, 3, 4, 5)], d),
                               significant(coln[c(1, 2, 3, 4, 5)], d), category = coln, ...)
  }
  if (!exists("out")) out <- "Plotting a Venn diagram of more than 5 analyses is not supported."
  # return(out)
}

#' @export
clusterDendogram <- function(analyses, alpha) {
  n <- length(analyses)
  if (n<2) stop("No sense in plotting a cluster dendrogram with just one analysis.")
  coln <- sapply(analyses, function(x) x@title)
  features = character(0)
  for (i in 1:n) {
    features <- c(features, rownames(analyses[[i]]@results))
  }
  features <- sort(unique(features))
  d <- data.frame(row.names = features, sapply(analyses, function(x) 1*(x@results[features,"q_value"]<=alpha)))
  names(d) <- coln
  d[is.na(d)] <- 0
  # analyze only those genes which are found to be DE at least once
  d <- d[rowMeans(d)>0, ]
  distance <- dist(t(d))
  plot(hclust(distance), hang = -1)
}

#' @export
barChart <- function(analyses, alpha, fill) {
  n <- length(analyses)
  # I hope to find an elegant solution to this in class definition
  if (n==1) {
    coln <- analyses@title
    features <- rownames(analyses@results)
    features <- sort(unique(features))
    d <- data.frame(row.names = features, 1*(analyses@results[features,"q_value"]<=alpha))
  } else {
    coln <- sapply(analyses, function(x) x@title)
    features = character(0)
    for (i in 1:n) {
      features <- c(features, rownames(analyses[[i]]@results))
    }
    features <- sort(unique(features))
    d <- data.frame(row.names = features, sapply(analyses, function(x) 1*(x@results[features,"q_value"]<=alpha)))
  }
  names(d) <- coln
  d[is.na(d)] <- 0
  o <- order(colSums(d), decreasing = FALSE)
  d <- colSums(d)[o]
  barplot(d, names.arg = NA, col = fill[o], legend.text = TRUE, main = "Number of differentially expressed genes",
          args.legend = list(x = "bottomright", text.width = max(strwidth(names(d)))))
}

#' @export
plotUpset <- function(analyses, alpha, ...) {
  coln <- sapply(analyses, function(x) x@title)
  sig <- lapply(analyses, function(x) rownames(x@results)[which(x@results < alpha)])
  names(sig) <- coln
  if (sum(sapply(sig, length) > 0) > 1) {
    UpSetR::upset(UpSetR::fromList(sig), order = "freq")
  } else stop("At least two non-empty sets are needed for an upset plot. There is ",
              sum(sapply(sig, length) > 0), ".")
}
