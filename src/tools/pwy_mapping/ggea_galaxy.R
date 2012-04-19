# 
# The runnable GGEA version for integration in Galaxy
# 
###############################################################################

## LIBRARIES
suppressWarnings(suppressPackageStartupMessages(library(MASS)))

if(!("GGEA.DIR" %in% objects())) GGEA.DIR <- "/home/cbio/ludwig/progs/ggea" 
setwd(GGEA.DIR)


## SOURCES
## SOURCE 1: Result HTML conversion
RES2HTML.PROG.LOC <- 
		"visualization/results2html.R"
source(RES2HTML.PROG.LOC)

## SOURCE 2: GGEA visualization
GGEAVIZ.PROG.LOC <-
		"visualization/ggeaviz.R"
source(GGEAVIZ.PROG.LOC)

### BUILD UP GGEA ENVIRONMENT
# SOURCE 3: load general ggea functions
source("ggea_func.R")
source("v01/ggea_v01_func.R")
source("v02/ggea_v02_func.R")
source("v03/ggea_v03_func.R")
source("v04/ggea_v04_func.R")

## FUNCTIONS
make.ggea.plots <- function(pwy, out.dir){
	## prepare ggea visualization by computation of render infos
	ggea.prep <- prepare.ggea.viz(pwy, ESET)
	## create a pdf device to plot into
	plot.file <- paste(	out.dir, 
						sub("path:", "", getName(pwy)),
						"_plot.pdf", sep="")
	pdf(plot.file)
	## (1) plot the whole pathway with all edges
	suppressWarnings(ggea.viz(ggea.prep = ggea.prep, what="all"))
	## (1b) plot the legend
	GGEA.graph.legend()
	## (2) plot the best edges of the pwy
	suppressWarnings(ggea.viz(ggea.prep = ggea.prep, what="best"))
	## (3) edge score & node degree statistics
	## get edge scores
	edge.scores <- as.numeric(ggea.prep$et["label",])
	## compute graph measures
	u.graph <- ugraph(ggea.prep$gr)
	d <- degree(u.graph)
	cc <- clusteringCoefficient(u.graph, selfLoops=T)
	av.cc <- round(mean(cc, na.rm=T), digits=2)
	## av.num.e <- round(aveNumEdges(ggea.prep$gr), digits=2)
	## num.no.e <- round(numNoEdges(ggea.prep$gr), digits=2)
	par(mfrow=c(2,2))
	## (3a) histogramms
	truehist(edge.scores, nbins=length(edge.scores)/5, ylim=c(0,1), 
			ylab="frequency",xlab="edge score", 
			main="Edge score distribution")
	truehist(d, nbins=length(d)/5, ylab="frequency",
			xlab="node degree", main="Node degree distribution")
	## (3b) boxplots
	boxplot(edge.scores, ylab="edge score")
	boxplot(d, ylab="node degree", 
				sub=paste("average cluster coef. =", av.cc))	
	## (4) hub analysis
	par(mfrow=c(2,1))
	## (4a) node degree vs p-value plot
	ps <- ggea.prep$nt$ADJ.P.BH
	ds <- d
	nas <- which(is.na(ggea.prep$nt$ADJ.P.BH))
	if(length(nas) > 0){
		ps <- ggea.prep$nt$ADJ.P.BH[-nas]
		ds <- d[-nas]
	}	
	plot(x=ds, y=ps, col="red", xlab="node degree", ylab="p value",
			main="Correlation between node degree & DE", 
			sub=paste("correlation coefficient =", 
					round(suppressWarnings(cor(x=ds,y=ps)), digits=2)))
	## (4b) node degree vs edge score plot
	aespn <- sapply(rownames(ggea.prep$nt), 
			function(x){ 
				ind <- grep(x, colnames(ggea.prep$et))
				scores <- as.numeric(ggea.prep$et["label", ind])
				return(mean(scores, na.rm=T))
			})
	ds <- d
	nas <- which(is.na(aespn))
	if(length(nas) >0){
		ds <- d[-nas]
		aespn <- aespn[-nas]
	}	
	plot(x=ds, y=aespn, col="red", 
			xlab="node degree", ylab="average edge score",
			main="Correlation between node degree & edge scores", 
			sub=paste("correlation coefficient =", 
					round(suppressWarnings(cor(x=ds,y=aespn)),
							digits=2)))
	
	## (5) ggea plot of the main hub, i.e. highest degree
	par(mfrow=c(1,1))
	main.hub <- mostEdges(ggea.prep$gr)
	suppressWarnings(ggea.viz(ggea.prep = ggea.prep, what=main.hub))
	## close the plotting device
	dev.off()
}

is.valid.pwy <- function(pwy){
	valid <- TRUE
	## A valid pathway ...
	## (a) ... is non metabolic, ... 
	valid <- valid & length(getReactions(pwy)) == 0
	## (b) ... has edges, ...
	valid <- valid & length(edges(pwy)) > 0
	## (c) ... has a valid node, ...
	valid <- valid & any(sapply(nodes(pwy), getType) %in% c("gene", "group"))
	## (d) ... and has a valid edge
	has.valid.edge <- valid & any(sapply(edges(pwy), 
							function(e){
								subty <- getSubtype(e)$subtype
								val <- ifelse(is.null(subty), NA,  getValue(subty))
								return(val)
							})	
							%in% c("-->", "--|", "---", "-+-"))
	
	return(valid)
}

### MAIN
main <- function(){
	# get command line arguments
	ESET.FILE <- commandArgs()[6]
	PWY.ZIP <- commandArgs()[7]
	FC.MAT.FILE <- commandArgs()[8]
	P.MAT.FILE <- commandArgs()[9]
	HTML.OUT <- commandArgs()[10]
	
	### GLOBAL SETTINGS
	assign("PMETHOD", "best", envir=.GlobalEnv)
	assign("GMETHOD", "mean", envir=.GlobalEnv)
	assign("MMETHOD", "mean", envir=.GlobalEnv)
	assign("SMETHOD", "bidirect", envir=.GlobalEnv)
	assign("FMETHOD", "unweigthed", envir=.GlobalEnv)
	
	# parse & load data
	assign("ESET", get(load(ESET.FILE)), envir=.GlobalEnv)
	## FC.MAT <- get(load(FC.MAT.FILE))
	## P.MAT <- get(load(P.MAT.FILE))
	
	# read in & parse pathways
	pwy.dir <- dirname(PWY.ZIP)
	unzip(PWY.ZIP, exdir=pwy.dir)
	pwy.files <- paste(pwy.dir, 
			as.vector(unzip(PWY.ZIP, list=TRUE)$Name), sep="/")
	pwys <- sapply(pwy.files, parseKGML)
	pwys <- pwys[sapply(pwys, is.valid.pwy)]
	
	# execute ggea
	res <- sapply(pwys, ggea, pvalue=FALSE)
	res <- t(res)
	res <- round(res, digits=3)
	
	# annotate and order result
	nn <- sapply(pwys, getName)
	tt <- sapply(pwys, getTitle)
	ann <- cbind(nn, tt)
	res <- cbind(ann, res)
	colnames(res) <- c("ID", "NAME", "RAW", "NORM")
	rownames(res) <- 1:length(pwys)
	rank <- order(as.numeric(res[,3]), decreasing=T)
	top.pwys <- pwys[rank][1:20]
	res <- res[rank,]
	## dummy p value
	dummy.ps <- function(raw.scores){
		while(raw.scores[1] > 10) raw.scores <- raw.scores / 10
		while(raw.scores[1] < 1) raw.scores <- raw.scores * 10
		ps <- 10^(-raw.scores)
		for(i in 1:length(ps)) if(ps[i] > 1) ps[i] <- 1
		ps <- signif(ps, digits=3)
		return(ps)
	}
	res <- cbind(res, dummy.ps(as.numeric(res[,"RAW"])))
	colnames(res)[5] <- "P.VALUE"
	
	out.dir <- sub("\\.dat","_files/", HTML.OUT)
	dir.create(out.dir)
	out.dir <- paste(out.dir, "GGEA_", sep="")
	
	## (a) create result file
	res.file <- paste(out.dir, "result.txt" ,sep="")
	write.table(res, file=res.file, row.names=F)
	## (b) create ggea plots
	sapply(top.pwys, make.ggea.plots, out.dir=out.dir)
	## (c) create argument log file
	tags <- c("eset", "pwy.zip", "fcmat", "pmat", "out")
	args <- commandArgs()[6:10]
	names(args) <- tags
	save(args, file=paste(out.dir, "args.RData", sep=""))
	
	## cleanup
	sapply(pwy.files, file.remove)
	
	# display output
	ggea2html(ggea.tbl=res, ggea.file=res.file, html.out=HTML.OUT)
	#return(list(gt=res, gf=res.file))
}

if(length(commandArgs()) - 1 && 
		length(grep("ggea_galaxy.R", commandArgs()[4]))) main()