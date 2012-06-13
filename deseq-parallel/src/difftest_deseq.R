library( DESeq )

### get arguments 1: INFILE, 2: OUTFILE 3:SIZE
args <- commandArgs()
INFILE<-args[4]
OUTFILE<-args[5]
SIZE<-args[6]

INFILE=c(paste(INFILE, "_COUNTS.tab", sep=""))

countsTable <- read.delim( INFILE,header=TRUE, stringsAsFactors=TRUE )

cname=countsTable[ , 1]
rownames( countsTable ) <- countsTable$gene
countsTable <- countsTable[ , -1 ]
countsTable=countsTable[,1:SIZE]

conds <- names(countsTable)

head( countsTable )

cds <- newCountDataSet( round(countsTable), conds )
head( counts(cds) )

cds <- estimateSizeFactors( cds )
sizeFactors( cds )
cds <- estimateVarianceFunctions( cds ,method='blind')

experiments <- names(countsTable)
names <- names(countsTable)

res<-c()
table_col_names<-c()
for (i in 1:(length(experiments)-1))
{
   for( j in (i+1):(length(experiments)))
   {
       print(c(i,j))
       tempres<- nbinomTest(cds,experiments[i],experiments[j])
       res=cbind(res,tempres[,7])
       #res=cbind(res,tempres[,8])
       table_col_names=cbind(table_col_names,paste(names[i],'_vs._',names[j],sep='')) 
   }
}

DiffTable<-res
rownames(DiffTable)<-rownames(countsTable)
colnames(DiffTable)<-table_col_names
write.table(DiffTable, file = OUTFILE, quote = FALSE, sep ="\t", eol ="\n", na = "1.000", dec = ".", row.names = TRUE,col.names =TRUE)
