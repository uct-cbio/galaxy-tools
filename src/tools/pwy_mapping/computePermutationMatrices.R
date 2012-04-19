# TODO: Add comment
# 
# Author: ludwig
###############################################################################

### LIBRARIES
suppressWarnings(suppressPackageStartupMessages(library(simpleaffy)))

### FUNCTIONS
recalc <- function(){
	eset.per <- ESET
	pData(eset.per)[[GROUP.COL]] <- gr[sample(1:le)]
	fc.tt <- get.fold.change.and.t.test(eset.per, GROUP.COL, labels)
	l <- list()
	l$fc <- fc.tt@fc
	#l$p <- fc.tt@tt
	l$p <- p.adjust(fc.tt@tt, method=ADJ.METHOD)
	return(l) 
}

### MAIN
# get command line arguments
ESET.FILE <- commandArgs()[6]
REPS <- commandArgs()[7]
GROUP.COL <- commandArgs()[8]
ADJ.METHOD <- commandArgs()[9]
FC.MAT.FILE <- commandArgs()[10]
P.MAT.FILE <- commandArgs()[11]

## load eset
ESET <- get(load(ESET.FILE))

## set params for permutation computation
reps <- as.numeric(REPS)
le <- nrow(pData(ESET))
gr <- as.vector(pData(ESET)[[GROUP.COL]])
labels <- levels(factor(gr))

## compute permutation matrices
fc.p <- replicate(reps, recalc())

## extract fold change permutations
FC.MAT <- sapply(1:reps, 
		function(x) unlist(fc.p[1,x]))
rownames(FC.MAT) <- featureNames(ESET)
save(FC.MAT, file=FC.MAT.FILE)
rm(FC.MAT)

## extract p value permutations
P.MAT <- sapply(1:reps, 
		function(x) unlist(fc.p[2,x]))
rownames(P.MAT) <- featureNames(ESET)
save(P.MAT, file=P.MAT.FILE)

## P.MAT.ADJ <- apply(P.MAT, 
##         MARGIN=2, FUN=p.adjust, method=ADJ.METHOD)
## P.MAT <- P.MAT.ADJ
## rownames(P.MAT) <- featureNames(ESET)
## save(P.MAT, file=paste(OUT.DIR, "/p_", 
##                 ADJ.METHOD, "_perm_mat.RData", sep=""))

