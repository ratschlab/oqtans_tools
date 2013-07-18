\name{getVarianceStabilizedData}
\Rdversion{1.1}
\alias{getVarianceStabilizedData}
\alias{varianceStabilizingTransformation}

\title{
   Apply a variance stabilizing transformation (VST) to the count data
}
\description{
   This function calculates a variance stabilizing transformation (VST) from the
   fitted dispersion-mean relation(s) and then transforms the count data (normalized
   by division by the size factor), yielding a matrix
   of values which are now approximately homoskedastic. This is useful as input 
   to statistical analyses requiring homoskedasticity. 
}
\usage{
varianceStabilizingTransformation(cds)
getVarianceStabilizedData(cds)
}
\arguments{
  \item{cds}{a \code{CountDataSet} which also contains the fitted dispersion-mean relation}
}
\details{
   For each sample (i.e., column of \code{counts(cds)}), the full variance function
   is calculated from the raw variance (by scaling according to the size factor and adding 
   the shot noise). The function requires a blind estimate of the variance function, i.e.,
   one ignoring conditions. Usually, this is achieved by calling \code{\link{estimateDispersions}}
   with \code{method="blind"} before calling it. 
   A typical workflow is shown in Section \emph{Variance stabilizing transformation} in the package vignette.

   If \code{\link{estimateDispersions}} was called with \code{fitType="parametric"},
   a closed-form expression for the variance stabilizing transformation is used on the normalized
   count data. The expression can be found in the file \file{vst.pdf} which is distributed with the vignette.

   If \code{\link{estimateDispersions}} was called with \code{fitType="locfit"},
   the reciprocal of the square root of the variance of the normalized counts, as derived
   from the dispersion fit, is then numerically
   integrated, and the integral (approximated by a spline function) is evaluated for each
   count value in the column, yielding a transformed value. 
   
   In both cases, the transformation is scaled such that for large
   counts, it becomes asymptotically (for large values) equal to the
   logarithm to base 2.

   Limitations: In order to preserve normalization, the same
   transformation has to be used for all samples. This results in the
   variance stabilizition to be only approximate. The more the size
   factors differ, the more residual dependence of the variance on the
   mean you will find in the transformed data. As shown in the vignette, you can use the function
   \code{meanSdPlot} from the package \pkg{vsn} to see whether this
   is a problem for your data.
}

\value{
   For \code{varianceStabilizingTransformation}, an \code{ExpressionSet}.

   For \code{getVarianceStabilizedData}, a \code{matrix} of
   the same dimension as the count data, containing the transformed
   values.  
} 


\author{Simon Anders <sanders@fs.tum.de>}

\examples{
cds <- makeExampleCountDataSet()
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cds )
colsA <- conditions(cds) == "A"
plot( rank( rowMeans( vsd[,colsA] ) ), genefilter::rowVars( vsd[,colsA] ) )
}
