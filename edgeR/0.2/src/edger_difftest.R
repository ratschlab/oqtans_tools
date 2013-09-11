### load edgeR package
suppressMessages(require("edgeR"))

### get arguments 1: INFILE, 2: OUTFILE 3:SIZE
args <- commandArgs()
ANLTYP<-args[4]
FDRMET<-args[5]
INFILE<-args[6]
OUTFILE<-args[7]

INFILE_COUNTS=c(paste(INFILE, "_COUNTS.tab", sep=""))
INFILE_CONDS=c(paste(INFILE, "_CONDITIONS.tab", sep=""))

### read count data from file
countsTable <- read.delim( INFILE_COUNTS, header=TRUE, stringsAsFactors=TRUE )
condsTable <- read.delim( INFILE_CONDS, header=TRUE, stringsAsFactors=TRUE )

tagnames <- countsTable[-(1:2), 1]
### use gene IDs as row names
rownames( countsTable ) <- countsTable$gene
countsTable <- countsTable[ , -1 ]
#head( countsTable )

conds <- factor( condsTable[ , 2] )

pw_tests <- list()
uniq_conds <- unique(conds)
for(i in 1:(length(uniq_conds)-1)) {
    for(j in (i+1):length(uniq_conds)) {
        pw_tests[[length(pw_tests)+1]] <- c(uniq_conds[i], uniq_conds[j])
    }
}

DGE <- DGEList(countsTable, group=conds)
tested <- list()

##Pairwise comparisons between two or more groups 
##(classic)
if (ANLTYP== "pw") {

    norm_factors <- calcNormFactors(as.matrix(DGE))
    disp <- estimateCommonDisp(DGE, rowsum.filter=5)
    disp <- estimateTagwiseDisp(disp, trend="movingave")

    for(i in 1:length(pw_tests)) {
        tested[[i]] <- exactTest(disp, pair=pw_tests[[i]])
        names(tested)[i] <- paste(pw_tests[[i]][2], "-", pw_tests[[i]][1], sep="")
    }
}

##Generalized Linear Models 
##
if (ANLTYP=="glm" | ANLTYP=="limma") {
    suppressMessages(require("splines"))
    cont <- NULL
    for(i in 1:length(pw_tests)) {
        cont <- c(cont, paste(pw_tests[[i]][2], "-", pw_tests[[i]][1], sep=""))
    }

    group_fact <- factor(conds)
    design <- model.matrix(~ -1 + group_fact)
    colnames(design) <- sub("conds_fact", "", colnames(design))
}
if (ANLTYP=="glm") {

    disp <- estimateGLMCommonDisp(DGE, design)
    disp <- estimateGLMTrendedDisp(disp, design)
    disp <- estimateGLMTagwiseDisp(disp, design)
    fit <- glmFit(disp, design)

    cont <- makeContrasts(contrasts=cont, levels=design)
    for(i in colnames(cont)) {
        tested[[i]] <- glmLRT(fit, contrast=cont[,i])
    }
}

##LIMMA Liner Models for RNA-seq
## 
if (ANLTYP=="limma") { 

    norm_factors <- calcNormFactors(as.matrix(DGE))
    y <- voom(DGE, design, plot=FALSE)

    fit <- lmFit(y, design)
    tab <- data.frame(ID=rownames(y$E), y$E, stringsAsFactors=F)

    cont <- NULL
    for(i in 1:length(pw_tests)) {
        cont <- c(cont, paste(pw_tests[[i]][2], "-", pw_tests[[i]][1], sep=""))
    }
    cont <- makeContrasts(contrasts=cont, levels=design)

    fit2 <- contrasts.fit(fit, cont)
    fit2 <- eBayes(fit2)

    tab <- NULL
    for(i in colnames(fit2)) {
        tab_tmp <- topTable(fit2, coef=i, n=Inf, sort.by="none", adjust.method=FDRMET)
        colnames(tab_tmp)[-1] <- paste(i, colnames(tab_tmp)[-1], sep=":")
        if(is.null(tab)) {
            tab <- tab_tmp
        } else tab <- cbind(tab, tab_tmp[,-1])
    }
    write.table(tab, OUTFILE, quote=F, sep="\t", row.names=F)
}
## 
## General result file 
if (ANLTYP=="glm" | ANLTYP=="pw") {
    tab <- NULL
    for(i in names(tested)) {
        tab_tmp <- topTags(tested[[i]], n=Inf, adjust.method=FDRMET)[[1]]
        colnames(tab_tmp) <- paste(i, colnames(tab_tmp), sep=":")
        tab_tmp <- tab_tmp[tagnames,]
        if(is.null(tab)) {
            tab <- tab_tmp
        } 
        else tab <- cbind(tab, tab_tmp)
    }
    tab <- cbind(Feature=rownames(tab), tab)
    write.table(tab, OUTFILE, quote=F, sep="\t", row.names=F)
}
##
