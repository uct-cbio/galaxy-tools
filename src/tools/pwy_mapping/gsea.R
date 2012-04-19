# 
# 
# Author: ludwig geistlinger
# Date: 14 June 2010
#
# wrapper script for gene set enrichment analysis
#
# CALL: Rscript gsea.R <eset> <pwys> <groupcol> <html.out>
#
###############################################################################

## LIBRARIES
suppressWarnings(suppressPackageStartupMessages(library(KEGGgraph)))
suppressWarnings(suppressPackageStartupMessages(library(Biobase)))

if(!("GGEA.DIR" %in% objects())) GGEA.DIR <- "/home/cbio/ludwig/progs/ggea" 
setwd(GGEA.DIR)

## SOURCES
#  GSEA R source program
GSEA.PROG.LOC <- "gsea/GSEA.1.0.R"

## Result HTML conversion
RES2HTML.PROG.LOC <- 
		"visualization/results2html.R"
#  GSEA R source program 
source(GSEA.PROG.LOC, verbose=F, max.deparse.length=9999)
## Result HTML conversion
source(RES2HTML.PROG.LOC)

## FUNCTIONS
build.doc.string <- function(eset.name, db.name){
	basenames <- sapply(c(eset.name, db.name), basename)
	basenames <- sapply(basenames, 
			function(x) unlist(strsplit(x, "\\."))[1])
	return(paste(basenames, collapse="_"))
}

## MAIN
## get command line arguments
ESET.FILE <- commandArgs()[6]
GS.GMT <- commandArgs()[7]
GROUP.COL <- commandArgs()[8]
AFF.COND <-commandArgs()[9]
ALPHA <- as.numeric(commandArgs()[10])
PERM <- as.numeric(commandArgs()[11])
HTML.OUT <- commandArgs()[12]

# load eset
eset <- get(load(ESET.FILE))

# build class list
cls <- list()
cls$phen <- levels(as.factor(as.vector(pData(eset)[,GROUP.COL])))
cls$class.v <- ifelse(
		as.vector(pData(eset)[,GROUP.COL]) == cls$phen[1], 0, 1)

out.dir <- sub("\\.dat","_files/", HTML.OUT)
dir.create(out.dir)

# run GSEA
GSEA(input.ds = as.data.frame(exprs(eset)), input.cls = cls,
		gs.db = GS.GMT, output.directory = out.dir,  
		doc.string            = build.doc.string(ESET.FILE, GS.GMT),
		non.interactive.run   = T,
		nperm                 = PERM,
		fdr.q.val.threshold   = ALPHA,
		topgs                 = 20,
		gs.size.threshold.min = 10
)

# Overlap and leading gene subset assignment analysis of the GSEA results

GSEA.Analyze.Sets(
		directory = out.dir, topgs = 20, height = 16, width = 16)

## create argument log file
tags <- c("eset", "gsdb", "grcol", "aff.cond", "alpha", "perm", "out")
args <- commandArgs()[6:12]
names(args) <- tags
save(args, file=paste(out.dir, "GSEA_args.RData", sep=""))

## display output
gsea2html(gsea.dir=out.dir, html.out=HTML.OUT, aff.cond=AFF.COND)


