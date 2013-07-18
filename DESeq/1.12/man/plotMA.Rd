\name{plotMA}

\alias{plotMA}
\title{Makes a so-called "MA-plot"}
\description{
A simple helper function that makes a so-called "MA-plot", i.e. a
scatter plot of logarithmic fold changes (on the y-axis) versus the mean of
normalized counts (on the x-axis).}
\usage{
plotMA(x, ylim,
  col = ifelse(x$padj>=0.1, "gray32", "red3"),
  linecol = "#ff000080",
  xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change),
  log = "x", cex=0.45, ...)
}
\arguments{
  \item{x}{a \code{data.frame} with columns \code{baseMean},
    and \code{log2FoldChange}. In addition, if the argument \code{col} is left at
    its default, this \code{data.frame} also needs to have a column named \code{padj}.}
  \item{linecol}{colour used for the horizontal line at y=0.}
  \item{ylim, col, xlab, ylab, log, cex, ...}{arguments that are passed on to \code{plot.default}.}
}

\details{This is a trivial helper function. Do not be afraid to edit and
    modify it to your needs.}

\value{
   The function is called for its side effect.
}
\author{
   Wolfgang Huber
}
\examples{
## see vignette
}
