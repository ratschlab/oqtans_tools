### load DESeq package
suppressMessages(require("DESeq"))

### get arguments 1: INFILE, 2: OUTFILE 3:SIZE
args <- commandArgs()
INFILE<-args[4]
OUTFILE<-args[5]

INFILE_COUNTS=c(paste(INFILE, "_COUNTS.tab", sep=""))
INFILE_CONDS=c(paste(INFILE, "_CONDITIONS.tab", sep=""))

### read count data from file
countsTable <- read.delim( INFILE_COUNTS, header=TRUE, stringsAsFactors=TRUE )
condsTable <- read.delim( INFILE_CONDS, header=TRUE, stringsAsFactors=TRUE )

### use gene IDs as row names
rownames( countsTable ) <- countsTable$gene
countsTable <- countsTable[ , -1 ]
head( countsTable )

conds <- factor( condsTable[ , 2] )
#head( countsTable )

cds <- newCountDataSet( round(countsTable), conds )
#head( counts(cds) )

cds <- estimateSizeFactors( cds )
#sizeFactors( cds )

### estimate variance function, use blind only, if no replicates are provided
if (length(levels(conds)) < length(conds))
{
    cds <- estimateDispersions( cds )
} else {
    writeLines("\n***You did not enter any replicates! - The results may be less valuable without replicates!***\n")
    cds <- estimateDispersions( cds, method='blind', sharingMode='fit-only')
}
experiments <- levels(conds)

res_1<-c()
res_2<-c()
res_3<-c()
res_4<-c()
res_5<-c()
res_6<-c()
res_7<-c()
res_8<-c()
table_col_names<-c()

for (i in 1:(length(experiments)-1))
{
   for( j in (i+1):(length(experiments)))
   {
       print(c(i,j))
       tempres <- nbinomTest(cds,experiments[i],experiments[j])
       res_1 = cbind(res_1,tempres[,1])
       res_2 = cbind(res_2,tempres[,2])
       res_3 = cbind(res_3,tempres[,3])
       res_4 = cbind(res_4,tempres[,4])
       res_5 = cbind(res_5,tempres[,5])
       res_6 = cbind(res_6,tempres[,6])
       res_7 = cbind(res_7,tempres[,7])
       res_8 = cbind(res_8,tempres[,8])
       table_col_names = cbind(table_col_names,paste('cond_', experiments[i], '_vs._cond_', experiments[j], sep='', 'test')) 
   }
}

DiffTable<-cbind(res_1,res_2,res_3,res_4,res_5,res_6,res_7,res_8)
colnames(DiffTable)<-c('feature ID', 'base  mean', 'base mean A', 'base mean B', 'fold change', 'log2 fold change','p value', 'adjusted p value')
write.table(DiffTable, file = OUTFILE, quote = FALSE, sep ="\t", eol ="\n", na = "1.000", dec = ".", row.names = TRUE,col.names =TRUE)
