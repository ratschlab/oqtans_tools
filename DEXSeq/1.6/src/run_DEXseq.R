### load DEXSeq package
suppressMessages(require("DEXSeq"))

### get arguments 1: INFILE, 2: OUTFILE 3:SIZE
args <- commandArgs()
INFILE<-args[4]
EXTRAPATH<-args[5]
annodb<-args[6]
OUTFILE<-args[7]

INFILE_CONDS=c(paste(INFILE, "_CONDITIONS.tab", sep=""))

### read count data from file
condsTable <- read.delim( INFILE_CONDS, header=TRUE, stringsAsFactors=TRUE, row.names=1 )
#head(condsTable)

conditions<-factor( condsTable[ , 1] )
#print(conditions)

## unique condition to define the pair of tests 
uniq_conds <- unique(conditions)
#print(uniq_conds)

## all possible pairs of conditions
pw_tests <- list()
for(i in 1:(length(uniq_conds)-1)) {
    for(j in (i+1):length(uniq_conds)) {
        pw_tests[[length(pw_tests)+1]] <- c(uniq_conds[i], uniq_conds[j])
    }
}
#print(pw_tests)

tab <- NULL
## testing all possible pairs of conditions
for(i in 1:length(pw_tests)) {
    ## header name 
    test_pair_name <- c(paste(pw_tests[[i]][1], "__vs__", pw_tests[[i]][2], sep=""))
    print(test_pair_name)

    sub.data <- subset(condsTable, (conditions %in% c(pw_tests[[i]][1],pw_tests[[i]][2])))
    sub.data[[1]]<-as.factor(sub.data[[1]])

    ecs = read.HTSeqCounts(countfiles=file.path(EXTRAPATH, row.names(sub.data)), design=sub.data, flattenedfile=annodb)
    
    ## Normalisation
    ecs <- estimateSizeFactors(ecs)
    print(sizeFactors(ecs))

    suppressMessages(require("parallel"))
    ## Dispersion estimation 
    ecs <- estimateDispersions(ecs, nCores=2)
    ecs <- fitDispersionFunction(ecs)

    ## Testing for differential exon usage
    ecs <- testForDEU(ecs, nCores=2)

    ## estimate the fold change 
    ecs <- estimatelog2FoldChanges(ecs, nCores=2)

    ## Result table
    resultTable <- DEUresultTable(ecs)
    print(head(resultTable))

    colnames(resultTable) <- paste(test_pair_name, colnames(resultTable), sep=":")
    #print(colnames(resultTable))

    if(is.null(tab)) {
        tab<- resultTable
    }
    else tab<- cbind(tab, resultTable)
}
## printing the result to out file 
write.table(tab, OUTFILE, quote=F, sep="\t", row.names=F)
