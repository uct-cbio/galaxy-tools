############################################################
#
# author: Ludwig Geistlinger
# date: 03 March 2010
# 
# CALL: Rscript <script> <cluster.table>
#
#
############################################################

### LIBRARY
suppressPackageStartupMessages(library(KEGG.db))

### INPUTS 
cl.table <- read.delim(commandArgs()[6])
cluster.col <- commandArgs()[7]
feat.col <- commandArgs()[8]
pathwy.col <- commandArgs()[9]
gene.col <- commandArgs()[10]

### MAIN
# reconstruct clusters
extr.cl <-  function(cl.index) as.vector(subset(cl.table, get(cluster.col)==cl.index, select=get(feat.col), drop=T))
x <- as.list(sapply(levels(cl.table[, cluster.col]), extr.cl))

# map probe 2 pathways
pathways <- as.vector(cl.table[ , pathwy.col])
names(pathways) <- as.vector(cl.table[ , feat.col])
probe2pathway <- sapply(as.list(pathways), function(i) unlist(strsplit(i, ",")))

# get pathway connections in clusters
is.in.pathway <- function(probe, pathway) pathway %in% probe2pathway[[probe]]
get.pathway.con <- function(cluster){
        # select probes of corresponding cluster
        pr <-  x[[cluster]]
        # check for pathways with more than one occurence in the cluster
        tbl <- table(unlist(sapply(pr, function(p) probe2pathway[[p]])))
        p.con <- names(tbl)[tbl > 1]
	if(any(p.con == "NA")) p.con <- p.con[-which(p.con == "NA")]

        # determine for these pathways the involved probes
        ret <- sapply(p.con, function(pwy) pr[sapply(pr, is.in.pathway, pathway=pwy)], simplify=F)
        # add pathID -> pathway name
        names(ret) <- paste(names(ret), sapply(names(ret), get, KEGGPATHID2NAME))
        # convert probeID -> geneID
        ret <-  sapply(ret, function(prs) paste(prs, cl.table[cl.table[,feat.col] %in% prs, gene.col]), simplify=F)
        return(ret)                                                                          
}



