# 
# 
# Author: ludwig geistlinger
# Date: 14 June 2010
#
# a wrapper script for executing PathWave
#
# Call: Rscript pathwave.R <eset> <pwys> <group.col> <res.out>
#
###############################################################################

## LIBRARIES
suppressWarnings(suppressPackageStartupMessages(library(PathWave)))

if(!("GGEA.DIR" %in% objects())) GGEA.DIR <- "/home/cbio/ludwig/progs/ggea" 
setwd(GGEA.DIR)

## CONSTANTS
OPT.GRID <- "gsea/optimalGridHSA.RData"

## SOURCES
## Result HTML conversion
RES2HTML.PROG.LOC <- 
		"visualization/results2html.R"

## PathWave visualization
PWVIZ.PROG.LOC <-
		"visualization/pathwaveviz.R"

## FUNCTIONS
make.pathwave.plots <- function(pwy, pwy.info, out.dir){
	## extract the pwy id
	pwy.id <- sub("path:", "", getName(pwy))
	## create a pdf device to plot into
	plot.file <- paste(out.dir, pwy.id, "_plot.pdf", sep="")
	pdf(plot.file)
	##(a) plot the reaction regulation volcano like
	par(xaxt='n')
	plot(x=pwy.info$reaction.regulation, 
			y=-log(pwy.info$reaction.p.value, base=10), col="red",
			main="Reaction regulation", 
			xlab="Direction of regulation",
			ylab="-log(reaction p.value)")
	par(xaxt='s')
	axis(side=1, at=c(-1,0,1), labels=c("down", "neutral", "up"), tick=T)
	##(b) feature significance
	truehist(pwy.info$feature.p.value, 
			main="Feature significance", xlab="feature p-value")
	##(c) pathwave plot 1 (all reactions)
	pathwave.viz(pwy=pwy,pw.result=pwy.info,what="all")
	## (d) general legend
	make.general.legend()
	##(e) pathwave plot 2 (best reactions)
	cmpds <- pathwave.viz(pwy=pwy,pw.result=pwy.info,what="best")
	## (f) particular legend
	make.particular.legend(cmpds)
	## close device
	dev.off()
}

#(a) Map entrez IDs on reactions
map.entrez.on.reactions <- function(keggXml){
	reac.entrez=NULL
	for(maps in names(keggXml$genes)){
		toDel=NULL   
		
		#Check if more than one gene maps to a reaction
		duplicate.reactions=
				unique(names(keggXml$genes[[maps]])[
								duplicated(names(keggXml$genes[[maps]]))])
		for(i in duplicate.reactions){
			
			reac.entrez=rbind(
					reac.entrez,c(i,paste((keggXml$genes[[maps]])[
											grep(paste("^",i,"$",sep=""),
													names(keggXml$genes[[maps]]))],collapse="_")))    
			toDel=c(toDel,grep(i,names(keggXml$genes[[maps]])))
		} 
		if(is.null(toDel)){
			rest.reactions=names(keggXml$genes[[maps]])
		} else{
			rest.reactions=names(keggXml$genes[[maps]])[-toDel]
		}
		for(i in rest.reactions){
			reac.entrez=rbind(
					reac.entrez,c(i,(keggXml$genes[[maps]])[
									grep(paste("^",i,"$",sep=""),
											names(keggXml$genes[[maps]]))]))
		}
	}
	index=order(reac.entrez[,1])
	reac.entrez=reac.entrez[index,]

	return(reac.entrez)
}

# (b) map expression values on entrez ids of reactions
map.expression.on.entrez <- function(reac.entrez, eset){
	x=NULL
	x.names=NULL
	for(i in 1:nrow(reac.entrez)){
		reac.entrezIDs=unlist(strsplit(reac.entrez[i,2],"_"))
		#Bind expression data together
		# select rows of expression matrix by entrez IDs that are
		# involved in the reaction
		rows=featureNames(eset) %in% reac.entrezIDs
		
		if(any(rows)){
			reac.exprs=exprs(eset)[rows,]
			
			#If more than one gene maps on a reaction: 
			# calculate mean value
			if(!is.null(dim(reac.exprs))) 
				reac.exprs <- apply(reac.exprs, 2, mean)
			
			x=rbind(x,reac.exprs)
			x.names=c(x.names,reac.entrez[i,1])
		}
	}
	rownames(x)=x.names
	return(x)
}
execute.pathwave <- 
		function(eset, pwy.dir, group.col, alpha, perm, diff.reac=5)
{
	## parse pathways
	keggXml <- pwKEGGxml(pwy.dir)
	# Build class factor
	y <- as.factor(as.vector(pData(eset)[,group.col]))
	#Map gene expression values onto reactions - START
	#(a) Map entrez IDs on reactions
	reac.entrez <- map.entrez.on.reactions(keggXml)
	# (b) map expression values on entrez ids of reactions
	x <- map.expression.on.entrez(reac.entrez, eset)
	#Map gene expression values onto reactions - END
	
	#Build adjacency matrices
	adj=pwAdjMatrices(keggXml$reactions,keggXml$compounds)
	
	#Get optimised grids
	#Download the RData file optimalGridHSA.RData and load it
	load(OPT.GRID)
	
	#Some reactions might be in x that are not in optimalGridHSA
	all.reac=unlist(lapply(optimalGridHSA,
					function(entry){unique(as.vector(entry$M))[
								unique(as.vector(entry$M))!="0"]}))
	all.reac=unique(all.reac)
	
	toDel=which(is.na(match(rownames(x),all.reac)))
	if(length(toDel)>0){
		x=x[-toDel,]
	}
	
	#Convert the expression values into Z-scores
	x=t(apply(x,1,function(entry){(entry-mean(entry))/sd(entry)}))
	
	print(paste("KEGG reactions:",length(all.reac),sep=" "))
	print(paste("reactions with expression values:",nrow(x),sep=" "))
	
	#Perform analysis
	result=pathWave(x, y, optimalM=optimalGridHSA, pvalCutoff=alpha, 
					genes=keggXml$genes, nperm=perm, diffReac=diff.reac)
	return(result)
}

## MAIN
main <- function(){
	## libraries
	suppressWarnings(suppressPackageStartupMessages(library(KEGGgraph)))
	suppressWarnings(suppressPackageStartupMessages(library(MASS)))
	
	## sources
	source(RES2HTML.PROG.LOC)
	source(PWVIZ.PROG.LOC)
	
	## get command line arguments
	ESET.FILE <- commandArgs()[6]
	PWY.ZIP <- commandArgs()[7]
	GROUP.COL <- commandArgs()[8]
	ALPHA <- as.numeric(commandArgs()[9])
	PERM <- as.numeric(commandArgs()[10])
	HTML.OUT <- commandArgs()[11]
	
	# load eset
	eset <- get(load(ESET.FILE))

	# read in & parse pathways
	pwy.dir <- dirname(PWY.ZIP)
	unzip(PWY.ZIP, exdir=pwy.dir)
	pwy.files <- paste(pwy.dir, 
			as.vector(unzip(PWY.ZIP, list=TRUE)$Name), sep="/")
	
	result <- execute.pathwave(eset, pwy.dir, GROUP.COL, ALPHA, PERM, 0)

	# parse pathway names and annotate them
	ind <- unlist(sapply(names(result), grep, x=pwy.files))
	pwys <- sapply(pwy.files[ind], parseKGML)
	tt <- sapply(pwys, getTitle)
	ann <- cbind(names(result), tt)
	
	colnames(ann) <- c("ID", "NAME")
	
	#get table 1.
	pval=NULL
	score=NULL
	all=NULL
	down=NULL
	up=NULL
	
	for(i in names(result)){
		pval=c(pval,result[[i]]$p.value)
		score=c(score,result[[i]]$score)
		all=c(all,length(getReactions(pwys[[grep(i,names(result))]])))
		up=c(up,length(result[[i]]$reaction.regulation[result[[i]]$reaction.regulation<0]))
		down=c(down,length(result[[i]]$reaction.regulation[result[[i]]$reaction.regulation>0]))
	}
	
	res=cbind(up+down,all,up,down,signif(pval, digits=3),score)
	rownames(res)=names(result)
	colnames(res) <- c("UP+DOWN", "ALL", "UP", "DOWN", "P.VALUE", "SCORE")

	res <- cbind(ann, res)
	
	out.dir <- sub("\\.dat","_files/", HTML.OUT)
	dir.create(out.dir)
	out.dir <- paste(out.dir, "PathWave_", sep="")
	
	res.file <- paste(out.dir, "result.txt" ,sep="")
	
	## OUTPUT
	## (0) save pathwave complete result
	save(result, file=paste(out.dir, "result.RData"))
	
	## (a) write result table
	write.table(res, file=res.file, row.names=F)
	
	## (b) create the plots
	if(length(ind) > 20) pwys <- pwys[1:20]
	sapply(1:length(pwys), function(i) 
				make.pathwave.plots(pwys[[i]], result[[i]], out.dir))
	sapply(pwy.files, file.remove)
	
	## (c) create argument log file
	tags <- c("eset", "pwy.zip", "grcol", "alpha", "perm", "out")
	args <- commandArgs()[6:11]
	names(args) <- tags
	save(args, file=paste(out.dir, "args.RData", sep=""))
	
	
	## (d) display html output
	pathwave2html(pw.tbl=res, pw.file=res.file, html.out=HTML.OUT)
}

if(length(commandArgs()) - 1 && 
		length(grep("pathwave.R", commandArgs()[4]))) main()
