### load DESeq package
suppressMessages(require("DESeq2"))

### get arguments 1: INFILE, 2: OUTFILE 3:SIZE
args <- commandArgs()
FITTYP<-args[4]
INFILE<-args[5]
OUTFILE<-args[6]

INFILE_COUNTS=c(paste(INFILE, "_COUNTS.tab", sep=""))
INFILE_CONDS=c(paste(INFILE, "_CONDITIONS.tab", sep=""))

### read count data from file
countsTable <- read.delim( INFILE_COUNTS, header=TRUE, stringsAsFactors=TRUE)
condsTable <- read.delim( INFILE_CONDS, header=TRUE, stringsAsFactors=TRUE)

tagnames <- countsTable[-(1:2), 1]
#print(tagnames)

## use gene IDs as row names
rownames( countsTable ) <- countsTable$gene
countsTable <- countsTable[ , -1 ]
#print(countsTable)

conditions<-factor( condsTable[ , 2] )
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
    #print(test_pair_name)
    ## colnames respective to the test pair 
    sub.data <- subset(condsTable, (conditions %in% c(pw_tests[[i]][1],pw_tests[[i]][2])))
    #print(sub.data)
    #print(sub.data[[1]]) # sample file name 
    #print(sub.data[[2]]) # condition 
    #print(sub.data[[3]]) # replicates 
    colData <- data.frame(row.names=sub.data[[1]], condition=sub.data[[2]], libType=sub.data[[3]])
    #print(colData)
    #print(countsTable[(sub.data[[1]])])
    dds <- DESeqDataSetFromMatrix(countData=countsTable[(sub.data[[1]])], colData=colData, design=~condition)
    colData(dds)$condition <- factor(colData(dds)$condition, levels=unique(sub.data[[2]]))
    dds <- DESeq(dds, fitType=FITTYP)
    ## concatenate the results
    tested_pairs <- results(dds) 
    #print(typeof(tested_pairs))
    
    colnames(tested_pairs) <- paste(test_pair_name, colnames(tested_pairs), sep=":") 
    #print(colnames(tested_pairs))
    #print(tested_pairs)
    
    tab_tmp <- tested_pairs[tagnames,]
    if(is.null(tab)) {
        tab<- as.data.frame(tab_tmp)
    }
    else tab<- cbind(tab, as.data.frame(tab_tmp))
}
## TODO cbind creates a X character to the string start place. 
colnames(tab) <- gsub("^X", "", colnames(tab))
## adding gene names to the row
tab <- cbind(Feature=row.names(tab_tmp), tab)
## priting the result
write.table(tab, OUTFILE, quote=F, sep="\t", row.names=F)
