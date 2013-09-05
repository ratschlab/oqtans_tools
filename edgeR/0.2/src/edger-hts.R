### get arguments 1: COUNT_MATRIX, 2: OUTFILE 3:SIZE
args <- commandArgs()
INFILE<-args[4]
OUTFILE<-args[5]

INFILE_COUNTS=c(paste(INFILE, "_COUNTS.tab", sep=""))

suppressMessages(library(edgeR))

toc <- read.table(INFILE_COUNTS, sep="\t", comment="", as.is=T)
groups <- sapply(toc[1, -1], strsplit, ":")
for(i in 1:length(groups)) { g <- make.names(groups[[i]][2]); names(groups)[i] <- g; groups[[i]] <- groups[[i]][-2] }
colnames(toc) <- make.names(toc[2,])
toc[,1] <- gsub(",", ".", toc[,1])
tagnames <- toc[-(1:2), 1]
rownames(toc) <- toc[,1]
toc <- toc[-(1:2), -1]
for(i in colnames(toc)) toc[, i] <- as.numeric(toc[,i])
norm_factors <- calcNormFactors(as.matrix(toc))

pw_tests <- list()
uniq_groups <- unique(names(groups))
for(i in 1:(length(uniq_groups)-1)) for(j in (i+1):length(uniq_groups)) pw_tests[[length(pw_tests)+1]] <- c(uniq_groups[i], uniq_groups[j])
DGE <- DGEList(toc, lib.size=norm_factors*colSums(toc), group=names(groups))

group_fact <- factor(names(groups))
design <- model.matrix(~ -1 + group_fact)
colnames(design) <- sub("group_fact", "", colnames(design))
disp <- estimateGLMCommonDisp(DGE, design)
disp <- estimateGLMTrendedDisp(disp, design)
disp <- estimateGLMTagwiseDisp(disp, design)
fit <- glmFit(disp, design)

tab <- data.frame(ID=rownames(fit$fitted.values), fit$fitted.values, stringsAsFactors=F)
write.table(tab, OUTFILE, quote=F, sep="\t", row.names=F)

## Under functional testing and active development of the code, not tested on a production server yet, 
