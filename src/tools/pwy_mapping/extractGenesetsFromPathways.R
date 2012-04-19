# 
# 
# Author: ludwig geistlinger
# Date: 16 June 2010
#
# script for extracting genesets from kgml files (KEGG XML)
# of pathways. The genesets are subsequently written in gmt format
# that is suitable for input to GSEA.
#
# CALL: Rscript extractGenesetsFromPathways.R <pwys.zip> <gmt.out>
#
###############################################################################

if(!("GGEA.DIR" %in% objects())) GGEA.DIR <- "/home/cbio/ludwig/progs/ggea" 
setwd(GGEA.DIR)

source("preprocessing/utils.R")

## LIBRARIES
load.pkgs.silently(c("KEGGgraph", "KEGGSOAP", "limma"))

## MAIN
## get command line arguments
PWY.ZIP <- commandArgs()[6]
GS.GMT <- commandArgs()[7]

# read in & parse pathways
pwys <- extract.pwys(PWY.ZIP)

# get pathway annotations
nn <- sapply(pwys, getName)
tt <- sapply(pwys, getTitle)

# fetch genes of pathway by SOAP functionality
genesets <- sapply(nn, function(x) 
			sub("hsa:", "", get.genes.by.pathway(x)), simplify=F)

# build first gmt column: the ID (format: <pwy.nr>_<pwy.title>)
nn <- sub("path:", "", nn)
n <- paste(nn, tt)
nw <- trimWhiteSpace(n)
n <- gsub(" ", "_", nw)

## collapse geneset members to one tab separated string 
gs.strings <- sapply(genesets, function(x) paste(x,collapse="\t"))

## paste an not annotated second column (artifact of gmt format)
ann <- paste(n, rep(NA,length(genesets)), sep="\t")

## paste all together
all <- paste(ann, gs.strings, sep="\t")

## collapse all together to a single newline separated string
all.str <- paste(all, collapse="\n")
all.str <- paste(all, "\n", sep="")

## write the genesets in gmt format
cat(all.str, file=GS.GMT, sep="")
