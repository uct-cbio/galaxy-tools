# 
# Author: ludqwig geistlinger
# Date: 14 June 2010 
#
# An implementation of overrepresentation analysis (ORA) based on
# subject sampling and a self-contained null hypothesis.
# The pseudocode is taken from Goeman and Buhlmann 
# (Bioinformatics 2007).
#
# CALL: Rscript ora.R <eset> <pmat> <pwys> <alpha> <res.out>
#
###############################################################################

## LIBRARIES
suppressWarnings(suppressPackageStartupMessages(library(genefilter)))
suppressWarnings(suppressPackageStartupMessages(library(multtest)))
suppressWarnings(suppressPackageStartupMessages(library(KEGGgraph)))
## libraries
suppressWarnings(suppressPackageStartupMessages(library(MASS)))

if(!("GGEA.DIR" %in% objects())) GGEA.DIR <- "/home/cbio/ludwig/progs/ggea" 
setwd(GGEA.DIR)

## CONSTANTS
#  GSEA R source program
GSEA.PROG.LOC <- "gsea/GSEA.1.0.R"

## Result HTML conversion
RES2HTML.PROG.LOC <- 
		"visualization/results2html.R"

## sources
source(RES2HTML.PROG.LOC)
source(GSEA.PROG.LOC, verbose=F, max.deparse.length=9999)

## FUNCTIONS

## overepresentation function: implements tuckeys higher criticism for
## computation of p-values on gene set level
tuckey.gene.set.p <- function(mgd, alpha){
	# count the number of p-values that are lower than alpha
	# in each permutation
	counts <- apply(p.mat, function(x) sum(x < alpha), MARGIN=2)
	# gene set set p-value: proportion of counts that are larger than
	# the number of de gene set genes
	p <- sum(counts > mgd) / ncol(p.mat)
	return(p)
} 

## parse geneset database
parse.genesets.from.GMT <- function(gs.gmt){
	content <- readLines(gs.gmt, warn=FALSE)
	le <- length(content)
	genesets <- list()
	for(i in 1:le){
		line <- content[i]
		spl <- unlist(strsplit(line, "\t"))
		name <- spl[1]
		spl <- spl[-c(1,2)]
		genesets[[name]] <- spl 
	}
	return(genesets)
}

# count the number of genes of the pathway that are among the
# de genes
count.mgd <- function(set, de.names){
	mgd <- sum(set %in% de.names)
	return(mgd)		
}

## plotting of course over all samples for significant genes
plot.gene.courses <- function(fDat, signif.genes){
	if(nrow(signif.genes) > 6) signif.genes <- signif.genes[1:6,]
	signif.names <- rownames(signif.genes)
	## order sample names according to affection status
	sample.names.ord <- sampleNames(eset)[
			order(pData(eset)[,group.col])]
	length.cond1 <- table(pData(eset)[,group.col])[1]
	conds <- levels(as.factor(pData(eset)[,group.col]))
	## extract the expression values of the significant genes					
	expr <- exprs(eset)[signif.names, sample.names.ord]
	min.e <- min(expr, na.rm=T)
	m <- matrix(rep(1:length(sampleNames(eset)),length(signif.names)),
			nrow=length(sampleNames(eset)),
			ncol=length(signif.names), byrow=F)
	## plot the expression courses of the significant genes
	if(ncol(m) > 1) matplot(x=m, y=t(expr), type="l", xlab="sample", 
				ylab="expression", main="Course of significant genes")
	else plot(x=as.vector(m), y=expr, type="l", xlab="sample",
				ylab="expression", main="Course of significant genes")
	lines(x=1:length.cond1, y=rep(min.e, length.cond1), lwd=2, col="blue")
	text(x=length.cond1 / 2, y=min.e-0.2, labels=conds[1])
	lines(x=(length.cond1+1):length(sampleNames(eset)),
			y=rep(min.e,length(sampleNames(eset))- 
							length.cond1), lwd=2)
	text(x=length.cond1 + (length(sampleNames(eset))- length.cond1)/2,
			y=min.e-0.2, labels=conds[2])
	legend("topright", legend=signif.genes$GENE, 
			col=1:length(signif.names), lty=1:length(signif.names))
}

## plot 2x2 ORA summary
create.plot <- function(fDat, plot.file){
	## open a pdf device to plot into
	pdf(file=plot.file)
	## (1) pval distribution
	truehist(as.numeric(fDat$ADJ.P.BH), nbins=nrow(fDat), 
				xlab="p value", main="P-Value distribution")
	## (2) volcano plot
	plot(x=as.numeric(fDat$FC), y=-log(as.numeric(fDat$ADJ.P.BH),
			base=10), col="red", xlab="fold change",
				ylab="-log p", main="Volcano Plot")
	## (3) heatmap
	expr <- exprs(eset)[rownames(fDat),]
	GSEA.HeatMapPlot(V = expr, row.names = as.vector(fDat$GENE), 
			col.labels = pData(eset)[,group.col], 
			col.classes = levels(as.factor(pData(eset)[,group.col])), 
			col.names = sampleNames(eset), 
			main =" Heat Map for Genes in Gene Set", xlab=" ", ylab=" ")
	## (4) gene curves
	## determine signficant genes
	signif.genes <- subset(fDat, IS.SIGNIFICANT)
	if(nrow(signif.genes) > 0) plot.gene.courses(fDat, signif.genes)
	## close the pdf device
	dev.off()
}

## write text report and plot ORA summary
write.report <- function(set, out.dir, org){
	set.id <- unlist(strsplit(names(set), "_"))[1]
	set <- set[[1]]
	fDat <- fData(eset)
	is.feature <- set %in% rownames(fDat)
	gene.ids <- set[is.feature]
	kegg.ids <- paste(org, gene.ids, sep=":")
	fDat <- cbind(kegg.ids, fDat[gene.ids,])
	fDat <- fDat[order(as.numeric(fDat$ADJ.P.BH)),]
	is.signif <- as.numeric(fDat$ADJ.P.BH) < alpha
	fDat <- cbind(fDat, is.signif)
	colnames(fDat)[1] <- "GENE"
	colnames(fDat)[length(colnames(fDat))] <- "IS.SIGNIFICANT"
	rep.file <- paste(out.dir, sub("path:","",set.id),
												"_report.txt", sep="")
	plot.file <- sub("report.txt", "plot.pdf", rep.file)
	create.plot(fDat=fDat, plot.file=plot.file)
	write.table(fDat, file=rep.file, row.names=F)
}

## ORA EXECUTION
ora <- function(eset, p.mat, genesets, alpha){
	assign("eset", eset, envir=.GlobalEnv)
	assign("p.mat", p.mat, envir=.GlobalEnv)
	
	# determine diff exp genes 
	de.genes <- subset(fData(eset), ADJ.P.BH < alpha)
	de.names <- rownames(de.genes)
	
	# determine mgd (nr of genes that are in gene set AND in diff exp genes) for each set
	mgds <- sapply(genesets, count.mgd, de.names=de.names)
	
	# determine pathway pvalues
	p.vals <- sapply(mgds, tuckey.gene.set.p, alpha=alpha)
	p.vals <- round(p.vals, digits=3)
	return(p.vals)
}

## MAIN
main <- function(){
	
	## get command line arguments
	ESET.FILE <- commandArgs()[6]
	p.mat.FILE <- commandArgs()[7]
	GS.GMT <- commandArgs()[8]
	assign("group.col", commandArgs()[9], env=.GlobalEnv)
	ALPHA <- commandArgs()[10]
	HTML.OUT <- commandArgs()[11]
	
	eset <- get(load(ESET.FILE))
	assign("eset", eset, env=.GlobalEnv)
	p.mat <- get(load(p.mat.FILE))
	assign("alpha", as.numeric(ALPHA), env=.GlobalEnv)

	## parse genesets
	genesets <- parse.genesets.from.GMT(GS.GMT)
	
	## determine organism
	org <- substr(names(genesets)[1], 1, 3)
	
	## execute ORA
	p.vals <- ora(eset, p.mat, genesets, alpha)
	
	# get pathway annotation
	set.names <- names(genesets)
	
	# design a data frame storing the pathway pval ranking
	set.ids <- substr(set.names, 1, 8)
	set.names <- sub("[a-z]{3}[0-9]{5}_", "", set.names)
	res <- cbind(set.names, p.vals)
	res <- cbind(set.ids, res)
	res <- as.data.frame(res)
	colnames(res) <- c("ID", "NAME", "P.VALUE")
	rank <- order(as.numeric(res$P.VALUE))
	le <- length(genesets)
	top.sets <- genesets[rank]
	if(length(top.sets) > 20) top.sets <- top.sets[1:20]
	res <- res[rank , ]
	rownames(res) <- 1:length(genesets)
	
	## create output files
	out.dir <- sub("\\.dat","_files/", HTML.OUT)
	dir.create(out.dir)
	out.dir <- paste(out.dir, "ORA_", sep="")
	## (a) create result file
	res.file <- paste(out.dir, "result.txt" ,sep="")
	write.table(res, file=res.file, row.names=F)
	## (b) create report and plot
	sapply(1:length(top.sets), function(i)
	         write.report(top.sets[i], out.dir=out.dir, org=org))
	## (c) create argument log file
	tags <- c("eset", "pmat", "gsdb", "grcol", "alpha", "out")
	args <- commandArgs()[6:11]
	names(args) <- tags
	save(args, file=paste(out.dir, "args.RData", sep=""))
 
	# display output
	ora2html(ora.tbl=res, ora.file=res.file, html.out=HTML.OUT)
}

if(length(commandArgs()) - 1 && 
		length(grep("ora.R", commandArgs()[4]))) main()
