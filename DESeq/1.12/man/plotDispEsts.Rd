\name{plotDispEsts}

\alias{plotDispEsts}
\title{Plot dispersion estimates and fitted values}
\description{
A simple helper function that plots the per-gene dispersion
estimates together with the fitted mean-dispersion relationship.
}
\usage{
plotDispEsts(cds, name=NULL, ymin, linecol="#ff000080",
  xlab = "mean of normalized counts", ylab = "dispersion",
  log = "xy", cex = 0.45, ... )
}
\arguments{
  \item{cds}{a \code{CountDataSet}.}
  \item{name}{this argument, together with \code{cds}, is passed on to \code{\link{fitInfo}}.}
  \item{ymin}{a scalar numeric, indicating the lower limit of the y-axis. The y-axis is plotted on the
    logarithmic scale. For the purpose of this plot, per-gene dispersion estimates that are below this
    value (in particular, those that happen to be zero) are shifted up
    to this value. If missing, the function attempts to guess a
    reasonable default.}
  \item{linecol}{colour used for the regression line}
  \item{xlab, ylab, log, cex, ...}{arguments that are passed on to \code{plot.default}.}
}

\details{This is a trivial helper function. Do not be afraid to edit and
    modify it to your needs.}

\value{
   The function is called for its side effect.
}
\author{
   Simon Anders, sanders@fs.tum.de
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
plotDispEsts(cds)
}
