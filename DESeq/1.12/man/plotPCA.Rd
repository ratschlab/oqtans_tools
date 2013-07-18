\name{plotPCA}
\alias{plotPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Sample PCA plot from variance-stabilized data
}
\description{This plot helps to check for batch effects and the like.}
\usage{
plotPCA(x, intgroup = "condition", ntop = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an ExpressionSet, as obtained from \code{\link{varianceStabilizingTransformation}}
}
  \item{intgroup}{

}
  \item{ntop}{ how many of the most variable genes should be used in calculating the PCA
}
}
\value{
a plot is produced
}
\author{
Wolfgang Huber}

\seealso{
\code{\link{varianceStabilizingTransformation}}
}
\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="blind" )
vsd <- varianceStabilizingTransformation( cds )
plotPCA( vsd )
}
