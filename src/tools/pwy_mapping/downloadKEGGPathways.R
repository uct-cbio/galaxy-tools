# TODO: Add comment
# 
# Author: Ludwig Geistlinger
# Date: 5th August 2010
#
# download all pathways of a particular organism from the KEGG FTP
#
# CALL: Rscript downloadKEGGPathways.R <organism> <out.zip.file>  
#
###############################################################################

suppressWarnings(suppressPackageStartupMessages(library(RCurl)))

KEGG.URL <- "ftp://ftp.genome.jp/pub/kegg/xml/kgml/"
META <- "metabolic/organisms/"

download.ftp.dir <- function(url, out.dir){
	filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
	filenames = paste(url, strsplit(filenames, "\r*\n")[[1]], sep = "")
	con = getCurlHandle( ftp.use.epsv = FALSE)
	contents = sapply(filenames, getURL, curl = con)
	out.names <- paste(out.dir, basename(filenames), sep="/")
	sapply(1:length(contents), function(i) cat(contents[i], file=out.names[i]))
}

download.pathways.of.organism <- function(org, out.dir){
	## download metabolic pathways
	meta.url <- paste(KEGG.URL, META, org, "/", sep="")
	download.ftp.dir(meta.url, out.dir)
	## download non-metabolic pathways
	non.meta.url <- paste(KEGG.URL, "non-", META, org, "/", sep="")
	download.ftp.dir(non.meta.url, out.dir)
}

retrieve.pathway.archive.of.organism <- function(org, out.file){
	out.dir <- dirname(out.file)
	out.dir <- sub("\\.dat","_files/", out.file)
	dir.create(out.dir)
	## (1) download pathways
	download.pathways.of.organism(org, out.dir)
	setwd(out.dir)
	pattern <- paste(out.dir, "*.xml", sep="/")
	## (2) zip pathways
	zip.file <- sub("dat$", "zip", out.file)
	system(paste("zip", zip.file, pattern))
	file.rename(from=zip.file, to=out.file)
	## (3) clean up 
	sapply(list.files(out.dir, full.names=TRUE, pattern="*.xml"), file.remove)
}

main <- function(){
	org <- commandArgs()[6]
	out.file <- commandArgs()[7]
	retrieve.pathway.archive.of.organism(org, out.file)
}

if(length(commandArgs()) - 1 && 
		length(grep("downloadKEGGPathways.R", commandArgs()[4]))) main()