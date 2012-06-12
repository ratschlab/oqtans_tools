library( DESeq )

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
    writeLines("\nYou did not enter any replicates! - The results may be less valuable without replicates!\n")
    cds <- estimateDispersions( cds, method='blind', sharingMode='fit-only')
}
experiments <- levels(conds)

res<-c()
table_col_names<-c()
for (i in 1:(length(experiments)-1))
{
   for( j in (i+1):(length(experiments)))
   {
       print(c(i,j))
       tempres <- nbinomTest(cds,experiments[i],experiments[j])
       res = cbind(res,tempres[,7])
       #res = cbind(res,tempres[,8])
       table_col_names = cbind(table_col_names,paste('cond_', experiments[i], '_vs._cond_', experiments[j], sep='')) 
   }
}

DiffTable<-res
rownames(DiffTable)<-rownames(countsTable)
colnames(DiffTable)<-table_col_names
write.table(DiffTable, file = OUTFILE, quote = FALSE, sep ="\t", eol ="\n", na = "1.000", dec = ".", row.names = TRUE,col.names =TRUE)
