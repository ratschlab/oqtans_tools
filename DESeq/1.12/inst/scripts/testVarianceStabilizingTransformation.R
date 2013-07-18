## W. Huber, 1.4.2012
## This script tests the function 'varianceStabilizingTransformation' on the
## pasillaGenes dataset. It also calls 'arrayQualityMetrics' to doublecheck 
## that the produced ExpressionSet can be used by that function.


library("pasilla")
library("vsn")
library("arrayQualityMetrics")

data("pasillaGenes")

cdsBlind = estimateDispersions(estimateSizeFactors(pasillaGenes), method="blind" )
pasillaVst = varianceStabilizingTransformation(cdsBlind)
 
meanSdPlot(pasillaVst)

arrayQualityMetrics(pasillaVst, intgroup=c("condition", "type"), force=TRUE)
