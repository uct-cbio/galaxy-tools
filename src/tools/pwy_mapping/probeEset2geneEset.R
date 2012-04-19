# 
# Author: ludwig geistlinger
# Date:	May 27th, 2010
#
# transforms a probe eset to a gene eset
###############################################################################

### LIBRARIES
suppressWarnings(suppressPackageStartupMessages(library(simpleaffy)))
suppressWarnings(suppressPackageStartupMessages(library(MASS)))

### FUNCTIONS

## fetch all probes that are annotated to a particular gene
get.probes.4.gene <- function(gene)
	featureNames(eset)[which(fData(eset)$ENTREZ == gene)]

## select the best, i.e. the one with lowest pvalue, probe of a gene
get.best.probe.4.gene <- function(gene){
	probes <- get.probes.4.gene(gene)
	ps <- fData(eset)[probes, "ADJ.P.BH"]
	return(probes[which.min(ps)])
}

## calculates expression level of one gene in one sample
probe.exprs.2.gene.exprs <- function(gene, sample, method="mean"){
	if(method=="mean") 
		gene.expr <- mean(as.numeric(exprs(eset)[gene2probes[[as.character(gene)]],sample]))
	if(method=="best") 
		gene.expr <- exprs(eset)[get.best.probe.4.gene(gene),sample]
	return(gene.expr)
}

## performs calculations that are necessary to summarize probe
## expressions to gene expression
probe.eset.2.gene.eset <- function(eset, method="mean")
{
	#(a) get the genes by ENTREZ ID
	genes <- unique(fData(eset)$ENTREZ)
	
	# (b) build a map gene -> belonging probes
	gene2probes <- sapply(genes, get.probes.4.gene)
	names(gene2probes) <- genes
	assign("gene2probes", gene2probes, envir=.GlobalEnv)
	
	
	# (c) compute gene expression matrix by averaging over spots
	e.genes <- sapply(sampleNames(eset), function(sample)
				sapply(genes, probe.exprs.2.gene.exprs, sample=sample, method=method))
	
	rownames(e.genes) <- genes
	
	# (d) delete NA 
	e.genes <- e.genes[-which(is.na(rownames(e.genes))), ]
	
	# (e) build a new Expression Set
	gene.eset <- new("ExpressionSet", exprs=e.genes)
	pData(gene.eset) <- pData(eset)
	return(gene.eset)
}

## annotates pData and fData of gene eset using the
## information that is provided by the probe eset
annotate <- function(probe.eset, gene.eset)
{
	pData(gene.eset) <- pData(probe.eset)
	curr <- unique(fData(probe.eset)$GENE.NAME)
	fData(gene.eset)$GENE.NAME <- curr[-which(is.na(curr))]
	curr <- unique(fData(probe.eset)$SYMBOL)
	fData(gene.eset)$SYMBOL <- curr[-which(is.na(curr))]
	return(gene.eset)
}

## performs differential expression analysis on gene eset
de.ana <- function(gene.eset, group.col, adj.meth)
{
	fc.tt <- get.fold.change.and.t.test(gene.eset, 
			group.col, levels(factor(gene.eset[[group.col]])))
	fData(gene.eset)[["FC"]] <- fc.tt@fc
	fData(gene.eset)[["RAW.P"]] <- fc.tt@tt
	fData(gene.eset)[[paste("ADJ.P", adj.meth, sep=".")]] <-
			p.adjust(fData(gene.eset)[["RAW.P"]], method=adj.meth)
	return(gene.eset)
}


### MAIN
# get command line arguments
PROBE.ESET.FILE <- commandArgs()[6]
GENE.ESET.FILE <- commandArgs()[7]
GROUP.COL <- commandArgs()[8]
ADJ.METH <- commandArgs()[9]
DISTR.FILE <- commandArgs()[10]
VOLC.FILE <- commandArgs()[11]

# load the probe eset ... 
eset <- get(load(PROBE.ESET.FILE))

# ... transform it to a gene eset ...
gene.eset <- probe.eset.2.gene.eset(eset)

# ... perform de analysis ...
gene.eset <- de.ana(gene.eset, GROUP.COL, ADJ.METH)

# ... and put out the gene eset
save(gene.eset, file=GENE.ESET.FILE)

# plot de
# (a) pval distribution
jpeg(file=DISTR.FILE, quality=100)
truehist(fData(gene.eset)$RAW.P, nbins=100, prob=T,
		main="P Value Distribution", xlab="P Value", ylab="Frequency")
dev.off()

# (b) volcano plot
jpeg(file=VOLC.FILE, quality=100)
plot(x=fData(gene.eset)$FC, y=-log(fData(gene.eset)$RAW.P,base=10), col="red",
		main="Volcano Plot", xlab="Fold Change", ylab="-log(p)")
dev.off()