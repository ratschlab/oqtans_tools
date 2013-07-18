plotMA = function(x, ylim,
  col = ifelse(x$padj>=0.1, "gray32", "red3"),
  linecol = "#ff000080",
  xlab = "mean of normalized counts", ylab = expression(log[2]~fold~change),
  log = "x", cex=0.45, ...)
{
  if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x))))
    stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")

  x = subset(x, baseMean!=0)
  py = x$log2FoldChange
  if(missing(ylim))
      ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  plot(x$baseMean, pmax(ylim[1], pmin(ylim[2], py)),
       log=log, pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  abline(h=0, lwd=4, col=linecol)
}

plotDispEsts = function( cds, name=NULL, ymin, linecol="#ff000080",
  xlab = "mean of normalized counts", ylab = "dispersion",
  log = "xy", cex = 0.45, ... )
{
  px = rowMeans( counts( cds, normalized=TRUE ) )
  sel = (px>0)
  px = px[sel]

  py = fitInfo(cds, name=name)$perGeneDispEsts[sel]
  if(missing(ymin))
      ymin = 10^floor(log10(min(py[py>0], na.rm=TRUE))-0.1)

  plot(px, pmax(py, ymin), xlab=xlab, ylab=ylab,
    log=log, pch=ifelse(py<ymin, 6, 16), cex=cex, ... )
  xg = 10^seq( -.5, 5, length.out=100 )
  lines( xg, fitInfo(cds, name=name)$dispFun( xg ), col=linecol, lwd=4)
}

plotPCA = function(x, intgroup="condition", ntop=500)
{
  rv = rowVars(exprs(x))
  select = order(rv, decreasing=TRUE)[seq_len(ntop)]
  pca = prcomp(t(exprs(x)[select,]))

  fac = factor(apply( pData(x)[, intgroup, drop=FALSE], 1, paste, collapse=" : "))
  if( length(fac) >= 3 )
     colours = brewer.pal(nlevels(fac), "Paired")
  else  
   colours = c( "green", "blue" )

  xyplot(PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), pch=16, cex=2,
    aspect = "iso", col=colours,
    main = draw.key(key = list(
      rect = list(col = colours),
      text = list(levels(fac)),
      rep = FALSE)))
}
