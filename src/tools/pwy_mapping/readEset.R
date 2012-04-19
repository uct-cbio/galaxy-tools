# TODO: Add comment
# 
# Author: ludwig geistlinger
# Date:	May 27th, 2010
#
# reads in data from plain text files and stores it in an expression
# set
###############################################################################

### LIBRARIES
suppressWarnings(suppressPackageStartupMessages(library(simpleaffy)))
suppressWarnings(suppressPackageStartupMessages(library(MASS)))

### MAIN
# get command line arguments
EXPRS.FILE <- commandArgs()[6]
PDAT.FILE <- commandArgs()[7]
FDAT.FILE <- commandArgs()[8]
GROUP.COL <- commandArgs()[9]
ADJ.METHOD <- commandArgs()[10]
OUT.FILE <- commandArgs()[11]
DISTR.FILE <- commandArgs()[12]
VOLC.FILE <- commandArgs()[13]

# read the eset
eset <- readExpressionSet(
		exprsFile=EXPRS.FILE, phenoDataFile=PDAT.FILE, 
		quote="\"'", colClasses="character", check.names=F)
fData(eset) <- read.table(FDAT.FILE, header=T, row.names=1L) 

# perform de analysis
fc.tt <- get.fold.change.and.t.test(
			eset, GROUP.COL, levels(factor(eset[[GROUP.COL]])))
fData(eset)[["FC"]] <- fc.tt@fc
fData(eset)[["RAW.P"]] <- fc.tt@tt
fData(eset)[[paste("ADJ.P", ADJ.METHOD, sep=".")]] <-
		p.adjust(fData(eset)[["RAW.P"]], method=ADJ.METHOD)

# output the eset
save(eset, file=OUT.FILE)

# plot de
# (a) pval distribution
jpeg(file=DISTR.FILE, quality=100)
truehist(fData(eset)$RAW.P, nbins=100, prob=T,
		main="P Value Distribution", xlab="P Value", ylab="Frequency")
dev.off()

# (b) volcano plot
jpeg(file=VOLC.FILE, quality=100)
plot(x=fData(eset)$FC, y=-log(fData(eset)$RAW.P,base=10), col="red",
		main="Volcano Plot", xlab="Fold Change", ylab="-log(p)")
dev.off()