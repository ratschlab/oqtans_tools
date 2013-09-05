### load DESeq package
suppressMessages(require("DESeq2"))

### get arguments 1: INFILE, 2: OUTFILE 3:SIZE
args <- commandArgs()
INFILE<-args[4]
OUTFILE<-args[5]

INFILE_COUNTS=c(paste(INFILE, "_COUNTS.tab", sep=""))
INFILE_CONDS=c(paste(INFILE, "_CONDITIONS.tab", sep=""))

### read count data from file
countsTable <- read.delim( INFILE_COUNTS, header=TRUE, stringsAsFactors=TRUE, row.names=1)
condsTable <- read.delim( INFILE_CONDS, header=TRUE, stringsAsFactors=TRUE, row.names=1)

colData <- data.frame(row.names=colnames(countsTable), condition=condsTable[[1]], libType=condsTable[[2]])

condition= factor(condsTable[[1]])

dds <- DESeqDataSetFromMatrix(countData=countsTable, colData=colData, design=~condition)

colData(dds)$condition <- factor(colData(dds)$condition, levels=unique(condsTable[[1]]))

dds <- DESeq(dds)

res <- results(dds)

write.csv(as.data.frame(res), file = OUTFILE)
