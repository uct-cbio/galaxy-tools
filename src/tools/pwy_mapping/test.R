HEADER <- "<html><head><title>Files for Composite Dataset </title></head><p/>This composite dataset is composed of the following files:<p/><ul>\n"
BODY <- "<li><a href=\"fname\">fname</a>\n"
TAIL <- "</ul></html>\n"

inp <- commandArgs()[6]
html <- commandArgs()[7]

comp.dir <- sub("\\.dat","_files/", html)
dir.create(comp.dir)

f1 <- paste(comp.dir,"res1.txt", sep="")
f2 <- paste(comp.dir, "res2.txt", sep="")

cat("hallo", file=f1)
cat("welt", file=f2)

cat(HEADER, file=html)
sapply(c(f1,f2), function(x) cat(gsub("fname",basename(x),BODY), file=html, append=TRUE))
cat(TAIL, file=html, append=TRUE)
