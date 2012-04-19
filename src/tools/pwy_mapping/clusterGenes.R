############################################################
#
# author: Ludwig Geistlinger
# date: 24 Feb 2010  
#
# CALL: Rscript	<eset>
#		<distance method>
#		<cluster method>
# 		<nr of res cluster> 
#		<outdir>
#
#
#
############################################################

### LIBRARY

suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(bioDist))

### INPUT
esetfile <- commandArgs()[6]
dist.meth <- commandArgs()[7]
clust.meth <- commandArgs()[8] 
nr.of.cluster <- as.integer(commandArgs()[9])
outdir <- commandArgs()[10]

### MAIN

eset <- get(load(esetfile))
x <- exprs(eset)

d <- switch(dist.meth,
		euclidean = euc(x),
		pearson = cor.dist(x),
		spearman = spearman.dist(x),
		kendall = tau.dist(x),
		manhattan = man(x),
		kullback.leibler.discrete = KLdist.matrix(x),
		kullback.leibler.continuous =  KLD.matrix(x),
		mutual.information.independence = mutualInfo(x),
		mutual.information.general = MIdist(x))

hc <- hclust(d, clust.meth)
plot(hc)
clusters <- rect.hclust(hc, k=nr.of.cluster, border="red")

is.gr <- function(i) ifelse(length(clusters[[i]]) > 2, TRUE, FALSE)
m2 <- which(sapply(seq(1,nr.of.cluster),is.gr))

clusters.m2 <- clusters[m2]


### OUTPUT

# (1) save clusters as R object
#save(clusters.m2, file=paste(outdir, "clusters.RData", sep="")

# (2) write clusters in a table file
cl.col <- rep(1:length(m2), times=sapply(clusters.m2, length))
feat.col <- names(unlist(clusters.m2))
fdat.mat <- fData(eset)[feat.col, ]
cl.table <- cbind(cl.col, feat.col)
cl.table <- cbind(cl.table, fdat.mat)

colnames(cl.table) <- c("CLUSTER", "FEATURE", colnames(fData(eset)))

write.table(as.matrix(cl.table), file=outdir, sep="\t", quote=F, row.names=F)


#write.clusters <- function(i) 
#	write.table(
#		as.matrix(fData(eset)[as.vector(clusters.m2[[i]]),]), 
#		file=paste(outdir, "cluster",  i, ".txt",sep=""), 
#		file=outdir, sep="\t", quote=F)
#sapply(seq(1, length(m2)), write.clusters)
#write.clusters(1)

