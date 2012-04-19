# 
# 
# Author: ludwig geistlinger
# Date: 12 July 2010
#
# ---ENRICHMENT BROWSER---
#
# seamless navigation through combined results 
# of set and graph enrichment
#
# CALL: Rscript enricher.R <gsea> <ora> <ggea> <pw> <result>
#
###############################################################################

if(!("GGEA.DIR" %in% objects())) GGEA.DIR <- "/home/cbio/ludwig/progs/ggea" 
setwd(GGEA.DIR)


## CONSTANTS
## (a) html titles
EBRO.TITLE 	<- "Seamless Navigation through Enrichment"
EBRO.NAV 	<- "Browse individual results"
EBRO.TABLE 	<- "Combined result table"

## (b) tool names
TOOLS 		<- c("ORA", "GSEA", "GGEA", "PathWave")
TOOLS.SHORT <- c("ora", "gsea", "ggea", "pw")

## (c) source code for tools
ORA.PROG.LOC 		<- "gsea/ora.R"
GSEA.PROG.LOC 		<- "gsea/gsea.R"
GGEA.PROG.LOC 		<- "ggea/ggea_galaxy.R"
PathWave.PROG.LOC 	<- "gsea/pathwave.R"
PROG.LOCS 			<- c(ORA.PROG.LOC, GSEA.PROG.LOC, 
							GGEA.PROG.LOC, PathWave.PROG.LOC)
names(PROG.LOCS)	<- TOOLS					

## (d) support columns					
SUPPORT.COLS 		<- c("P.VALUE", "RANK", "GRAPHICS")
names(SUPPORT.COLS) <- c("pval", "rank", "graph")

## Result HTML conversion
RES2HTML.PROG.LOC <- 
		"/home/cbio/ludwig/progs/ggea/visualization/results2html.R"
source(RES2HTML.PROG.LOC)

## trim gsea top table to relevant columns for combined view 
trim.gsea.ttable <- function(gsea.ttable){
	## collapse the add on graphics together
	addon.cols <- paste("GSEA", c("REPORT", "PLOT", "BROWSE"), sep=".")
	collapsed.addon <- apply(gsea.ttable[,addon.cols], FUN=paste, 
			MARGIN=1, collapse="")
	## create the trimmed version 
	trimmed <- cbind(gsea.ttable[,1:2], 
			signif(as.numeric(gsea.ttable[,"FDR.q.val"]), digits=2))
	trimmed <- cbind(trimmed, 1:nrow(gsea.ttable))
	trimmed <- cbind(trimmed, collapsed.addon)
	colnames(trimmed) <- 
			c("ID", "NAME", paste("GSEA", SUPPORT.COLS, sep="."))
	return(trimmed)
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

stouffer.merge <- function(pvals, two.sided=TRUE){
	## remove NA
	pvals <- pvals[!is.na(pvals)]
	## correct for two sided tests
	if(two.sided) pvals <- pvals/2
	## z transformation
	z <- sum(-qnorm(pvals)) / sqrt(length(pvals))
	## backtransformation
	p <- 1 - pnorm(abs(z))
	## correct for two sided tests
	if(two.sided) p <- p * 2
	return(p)
}

## create ora graphics (report, plot & browse) 
## for eBrowse ids that are not in the ora top table
create.missed.graphics <- function(tool, id, ids, out.dir, add.args)
{
	ind <- which(ids == id)
	pwy <- pwys[ind]
	
	## (1) create tool specific graphics
	report <- browse <- FALSE
	plot <- TRUE
	
	if(tool == "GGEA")
		make.ggea.plots(pwy=pwy[[1]], out.dir=out.dir)
	else if(tool == "PathWave")
		make.pathwave.plots(pwy=pwy[[1]], out.dir=out.dir,
									pwy.info=add.args$pwy.info[[id]])
	else if(tool == "ORA"){
		write.report(set=pwy, out.dir=out.dir, org=add.args$org)
		report <- TRUE
	}
	
	## (c) create html links to the graphics
	report.col <- plot.col <- browse.col <- ""

	if(report){
		report.file <- 
				paste(out.dir, tool, "_", id, "_report.txt" ,sep="")
		file.rename(from=paste(out.dir, id, "_report.txt" ,sep=""),
					to=report.file)
		report.col <- make.href.col(ref=basename(report.file), 
									tag=REPORT.TAG[2])	
	}
	if(plot){
		plot.file <- paste(out.dir, tool, "_", id, "_plot.pdf" ,sep="")
		file.rename(from=paste(out.dir, id, "_plot.pdf" ,sep=""),
					to=plot.file)
		plot.col <- make.href.col(	ref=basename(plot.file), 
									tag=PLOT.TAG[2])	
	}	
	if(browse){
		url <- browse.ora(report.file)
		browse.col <- make.href.col(ref=url, tag=BROWSE.TAG)
	}	
	
	return(paste(report.col, plot.col, browse.col, sep=""))
}

## 	extract already computed information for shared IDs
extract.shared.support <- function(tool, shared.ids, ttbl){
	## (a) determine ranks
	ranks.shared <- sapply(shared.ids, 
			function(id) which(ttbl[,"ID"] == id))
	## (b) extract pvalues and graphics
	pval.and.graphics.cols <- 
			c("P.VALUE", grep(tool, colnames(ttbl), value=TRUE))
	pvals.and.graphics.shared <- 
			as.matrix(ttbl[ranks.shared, pval.and.graphics.cols])
	## (c) bind pvalues and ranks as first two cols of the table
	shared.tbl <- cbind(as.vector(
					pvals.and.graphics.shared[,1]), 
			ranks.shared)
	## (d) collapse graphics together as a single col
	le <- length(pval.and.graphics.cols)
	if(le > 2) shared.graphics.col <- 
						apply(pvals.and.graphics.shared[,2:le], 
									FUN=paste, MARGIN=1, collapse=" ")
	else shared.graphics.col <- pvals.and.graphics.shared[,2]			
	## (e) bind graphics to table (third col)						
	shared.tbl <- cbind(shared.tbl, shared.graphics.col)
	## (f) give names for rows & cols
	rownames(shared.tbl) <- shared.ids
	colnames(shared.tbl) <- paste(tool, SUPPORT.COLS, sep=".")
	
	return(shared.tbl)
}

## compute missing information for IDs that are not in the top table
compute.missing.support <- 
		function(tool, miss.ids, gtbl, out.dir, pwy.ids, add.args)
{
	miss.tbl <- NULL
	
	## determine ids which are available in the global result
	available <- miss.ids %in% gtbl$ID
	avail.ids <- miss.ids[available]
	
	if(any(available)){	
		## (a) grep p.value and ranking from the global result
		ranks.miss <- sapply(avail.ids, grep, x=gtbl[,"ID"])
		pvals.miss <- gtbl[ranks.miss, "P.VALUE"] 
					
		## (b) create add on (report, plot & browse)
		if(tool=="PathWave") 
			add.args <- list(pwy.info=add.args$pwres[avail.ids])
		
		miss.graphics.refs <- sapply(avail.ids, create.missed.graphics, 
						tool=tool, ids=pwy.ids,  out.dir=out.dir, 
						add.args=add.args)
		
		## (c) pvals and ranks as first and second col
		miss.tbl <- cbind(pvals.miss, ranks.miss) 
		## (d) bind graphics as third col
		miss.tbl <- cbind(miss.tbl, miss.graphics.refs)
	}
	## set NA for not available IDs
	le <- length(miss.ids) - sum(available)
	if(le > 0) for(i in 1:le) miss.tbl <- rbind(miss.tbl, rep(NA, 3))
	
	## (e) give names for rows & cols
	rownames(miss.tbl) <- c(avail.ids, miss.ids[!available])
	miss.tbl <- miss.tbl[miss.ids,]
	colnames(miss.tbl) <- paste(tool, SUPPORT.COLS, sep=".")
	
	return(miss.tbl)
}

add.tool.support <- 
		function(tool, eBro.ids, ttbl, gtbl, out.dir, pwy.ids, add.args)
{
	## get tool specific functions
	source(PROG.LOCS[tool])
	## determine which eBrowse IDs are 
	## (a) in the top table of the tool (shared)
	shared <- eBro.ids %in% ttbl[,"ID"]
	shared.ids <- eBro.ids[shared]
	shared.tbl <- extract.shared.support(tool, shared.ids, ttbl)
	
	## (b) not in the top table of the tool (missing)
	miss.ids <- eBro.ids[!shared]
	miss.tbl <- compute.missing.support(tool, miss.ids, gtbl, 
										out.dir, pwy.ids, add.args)
	
	## bring information for shared & missing together
	support.tbl <- rbind(shared.tbl, miss.tbl)
	support.tbl <- support.tbl[eBro.ids,]
	
	return(support.tbl)	
}

## create navigation bar for browsing through results
## of each tool
create.navigation <- function(results){
	## init navigation table
	elems <- c(
			paste("<h4 class=\"title30\">", EBRO.NAV, "</h4>", sep=""), 
			"<table width=\"800\">", "<tr>")
	images <- list.files(IMAGE.LOC)
	## add a navigation cell for each tool 
	for(tool in TOOLS){
		## get the icon of current tool
		tool.icon <- grep(tool, images, value=TRUE)
		tool.icon <- paste(IMAGE.REM, tool.icon, sep="")
		## link the result
		elems <- c(elems, "<td width=\"200\" align=\"center\">")
		elems <- c(elems, paste(
					"<a href=\"", results[tool], 
					"\"><img src=\"", tool.icon, 
					"\" alt=\"", tool, 
					"\" height=\"60\" width=\"120\" border=\"0\"/>",
					"</a>", 
					sep=""))
		elems <- c(elems, "</td>")
	}						
	elems <- c(elems, "</tr>")
	elems <- c(elems, "</table>")
	elems <- c(elems, "<br />")
	return(paste(elems, collapse="\n"))
}

## create the eBrowse html interface
eBro2html <- function(out.dir, html.out, results){
	header <- create.html.header(tool="eBro", title=EBRO.TITLE)
	cat(header, file=html.out)
	navigation <- create.navigation(results)
	cat(navigation, file=html.out, append=TRUE)
	cat(paste("<h4 class=\"title30\">", EBRO.TABLE, "</h4>", sep=""),
			file=html.out, append=TRUE)
	combined <- combine.results(out.dir)
	cat(combined, file=html.out, append=TRUE)
	cat("</body></html>", file=html.out, append=TRUE)
}

## load additional RData objects for a set of <tools>
## located in <out.dir> and having extension <ext>
load.tools.data <- function(out.dir, tools, ext){
	data.files <- sapply(tools, function(tool) 
				paste(out.dir, tool, ext, sep=""))
	data <- sapply(data.files, 
			function(x) get(load(x)), simplify=F)
	names(data) <- tools
	return(data)
}

## read the global result tables for specified <tools>
## located in <out.dir> and having extension <ext>
read.gtables <- function(out.dir, tools, ext){
	gtable.files <- sapply(tools, function(tool) 
				paste(out.dir, tool, ext, sep=""))
	gtables <- sapply(gtable.files, 
						function(x) read.table(x, header=T))
	names(gtables) <- tools
	return(gtables)
}

## load arguments from the GSEA log file
load.gsea.args <- function(gsea.log){
	assign("eset", get(load(gsea.log["eset"])), envir=.GlobalEnv)
	assign("pwys", 
			parse.genesets.from.GMT(gsea.log["gsdb"]), envir=.GlobalEnv)
	assign("group.col", gsea.log["grcol"], envir=.GlobalEnv)
}

## read pathways for GGEA & PathWave
read.pathways <- function(pwy.zip){
	pwy.dir <- dirname(pwy.zip)
	unzip(pwy.zip, exdir=pwy.dir)
	pwy.files <- paste(pwy.dir, 
			as.vector(unzip(pwy.zip, list=TRUE)$Name), sep="/")
	assign("pwys", sapply(pwy.files, parseKGML), envir=.GlobalEnv)
	sapply(pwy.files, file.remove)
}

## load a particular argument having value <value> and the tag <key>
load.arg <- function(key, value){
	if(key == "eset")
		assign("eset", get(load(value)), envir=.GlobalEnv)
	if(key == "gsdb")
		assign("pwys", parse.genesets.from.GMT(value), envir=.GlobalEnv)
	if(key == "pwy.zip") 
		read.pathways(value)
	if(key == "grcol")
		assign("group.col", value, envir=.GlobalEnv)
	if(key == "alpha")
		assign("alpha", as.numeric(value), envir=.GlobalEnv)
}

## load all non-existing arguments of the specified tool
load.args <- function(tool, logs){
	log <- logs[[tool]]
	nam <- names(log)
	for(n in nam) if(!exists(n)) load.arg(n, log[n])
}

compute.support <- 
		function(tool, ttables, logs, out.dir, gtables, eBrowse.IDs)
{
	## load the arguments
	load.args(tool, logs)
	
	## tool specific settings
	if(tool=="ORA"){
		pwy.ids <- substr(names(pwys), 1, 8)
		add.args <- list(org = substr(names(pwys)[1], 1, 3))
	}
	else if(tool=="GGEA"){
		pwy.ids <- sub("path:", "", sapply(pwys, getName)) 
		ttables$GGEA[,"ID"] <- sub("path:", "", ttables$GGEA[,"ID"])
		gtables$GGEA$ID <- sub("path:", "", gtables$GGEA$ID)
		assign("ESET", eset, envir=.GlobalEnv)	
		add.args <- list()
	}
	else if(tool=="PathWave"){
		pwy.ids <- sub("path:", "", sapply(pwys, getName))
		pw.result <- get(load(
						paste(out.dir, "PathWave_result.RData", sep="")))
		add.args <- list(pwres = pw.result)
	}
	
	## create support columns for tool
	scols <- paste(tool, SUPPORT.COLS, sep=".")
	names(scols) <- names(SUPPORT.COLS)
	
	## compute the support table
	## (a) for ORA, GGEA & PathWave
	if(tool != "GSEA") supp.tbl <- add.tool.support(
							tool=tool, eBro.ids=eBrowse.IDs, 
							ttbl=ttables[[tool]], gtbl=gtables[[tool]], 
							out.dir=out.dir, pwy.ids=pwy.ids,
							add.args=add.args)
	## (b) for GSEA			
	else supp.tbl <- ttables$GSEA[,scols]			
	
	## store pvals & ranks for averaged/combined measure				
	pvals <- as.numeric(supp.tbl[,scols["pval"]])
	ranks <- as.numeric(supp.tbl[,scols["rank"]])
	
	## transform NA's to empty strings
	for(i in 1:nrow(supp.tbl)) 
		if(all(is.na(supp.tbl[i,]))) supp.tbl[i,] <- rep("", 3)
	
	## add a line break before the graphics
	supp.tbl[,scols["graph"]] <- 
			paste("<br/>", supp.tbl[,scols["graph"]], sep="")
	
	## write ranks in brackets for collapsed version
	supp.tbl[,scols["rank"]] <- 
			paste("(", supp.tbl[,scols["rank"]], ")", sep=".")
	
	## collapse ranks, pvals & graphics refs to a single column
	supp.col <- apply(supp.tbl, FUN=paste, MARGIN=1, collapse=" ")
	
	return(list(scol=supp.col, p=pvals, r=ranks))							
}

combine.results <- function(out.dir){
	## load the top tables of all tools
	ttables <- load.tools.data(out.dir, TOOLS, "_rtable.RData")
	## ttables <- sapply(ttables, as.data.frame)
	
	## load the command args log files of all tools
	logs <- load.tools.data(out.dir, TOOLS, "_args.RData")

	## read global result tables for all non-gsea tools
	gtables <- read.gtables(out.dir, TOOLS[c(1,3,4)], "_result.txt")
	
	## init combined table with top IDs (= eBrowse IDs) of GSEA run
	gsea.table <- trim.gsea.ttable(ttables$GSEA)
	comb.table <- subset(gsea.table, select=c("ID", "NAME"))
	eBrowse.IDs <- comb.table[,"ID"]
	ttables$GSEA <- gsea.table
	
	## compute & add support for each tool
	supports <- sapply(TOOLS, compute.support, ttables=ttables,
						logs=logs, out.dir=out.dir, gtables=gtables, 
								eBrowse.IDs=eBrowse.IDs, simplify=F)
	ranks <- NULL
	pvals <- NULL
	for(t in TOOLS[c(2,1,3,4)]){ 
		comb.table <- cbind(comb.table, supports[[t]][["scol"]])
		ranks <- cbind(ranks, supports[[t]][["r"]])
		pvals <- cbind(pvals, supports[[t]][["p"]])
	} 
	colnames(comb.table)[3:6] <- 
			paste(TOOLS[c(2,1,3,4)], "SUPPORT", sep=".")
	## ------------------
	## DUMMY TABLE
	## for(i in 1:4) comb.table <- cbind(comb.table, gsea.support.col)
	## colnames(comb.table)[3:6] <- 
	##         paste(TOOLS[c(2,1,3,4)], "SUPPORT", sep=".")
	## comb.table <- cbind(comb.table, 1:20)
	## comb.table <- cbind(comb.table, pvals)
	## colnames(comb.table)[7:8] <- c("AV.RANK", "COMB.PVALUE")
	## DUMMY TABLE
	## -------------------
	
	## compute average ranks & pvals
	av.ranks <- apply(ranks, MARGIN=1, FUN=mean, na.rm=T)
	comb.pvals <- apply(pvals, MARGIN=1, FUN=stouffer.merge)
	
	## bind average measures to table
	comb.table <- cbind(comb.table, 
						signif(cbind(av.ranks, comb.pvals), digits=3))
	colnames(comb.table)[ncol(comb.table) - c(1,0)] <- 
			c("AV.RANK", "COMB.PVALUE")
	
	## transform R table to HTML table
	comb.html.table <- r.table.2.html.table(comb.table)
	
	return(comb.html.table)
}

## MAIN
main <- function(){
	## get command line arguments
	gsea.result <- commandArgs()[6]
	ora.result <- commandArgs()[7]
	ggea.result <- commandArgs()[8]
	pw.result <- commandArgs()[9]
	enrichment.window <- commandArgs()[10]
	
	results <- c(ora.result, gsea.result, ggea.result, pw.result)
	
	## create extra files directory for eBro
	out.dir <- sub("\\.dat","_files/", enrichment.window)
	dir.create(out.dir)
	
	## copy html result files of all tools
	html.results <- sub("\\.dat", ".html", basename(results))
	file.copy(from=results, to=paste(out.dir, html.results, sep=""))
	
	## copy all result files of all tools to eBro extra dir
	sapply(TOOLS, 
			function(tool){
				tool.dir <- sub("\\.dat","_files/", results)
				tfiles <- list.files(tool.dir, full.names=TRUE)
				efiles <- paste(out.dir, basename(tfiles), sep="")
				file.copy(from=tfiles, to=efiles)
			})
	
	## eBro 2 html
	names(html.results) <- TOOLS
	eBro2html(out.dir, enrichment.window, html.results)		
}

if(length(commandArgs()) - 1 && 
		length(grep("eBro.R", commandArgs()[4]))) main()






