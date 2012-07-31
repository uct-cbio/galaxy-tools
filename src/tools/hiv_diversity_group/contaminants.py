#!/usr/bin/python
#check for matching seqs in the g-drive folder using pairwise comparisons
#to run example:: python contaminants.py contamtest.fas gag contaminantsOUT.txt
import os, sys, re

###function to read fasta files##################
def fastarecords(file):         #file should be opened
        def record(l):          #get title, seq record ..excludes last record in file, read below
                title = lines[l].strip()[1:]
                seq = ""
                while 1:
                        l = l+1
                        if not lines[l].startswith(">"):
                                seq = seq + lines[l].strip()
                        else:
                                break
                else:
                                pass
                return(title,seq)

        lines = file.readlines()
	seqd = []
        for l in range(len(lines)):
               	if lines[l].startswith(">"):
                       	seqd.append(l)
	records = {}
        for sd in seqd[:-1]:                    #reads all but the last record
        	drecord = record(sd)
		records[drecord[0]] = drecord[-1]
# get data for the last record in the fasta file
        lasttitle = lines[seqd[-1]].strip()[1:]
	lseq = ""
        for ls in range(seqd[-1]+1,len(lines)):
       	       	lseq = lseq+lines[ls].strip()
	records[lasttitle] = lseq
        return records
####################################################

###function to calculate Hamming distance###########
def hamdist(str1, str2):		#use after aligning the seqs
	#counts the number of differences btwn equal length str1 and str2
	diffs = 0
	for ch1, ch2 in zip(str1, str2):
		if ch1 != ch2:
			diffs +=1
	return float(diffs)/len(str1)	
####################################################

infile, genename, outf = sys.argv[-3],sys.argv[-2], sys.argv[-1]		#infile is fasta file of new sequences

###get input data 
ff = open(infile, "r")  
newrecords = fastarecords(ff)
titlelist = []
for newrecord in newrecords:
	if newrecord not in titlelist:
		titlelist.append(newrecord)

#outfile = open(infile[:-4]+"_contaminantsOUT.txt","w")		#with details of possible contaminants
outfile = open(outf,"w")
hammings = []			#list all hamming distances <=0.3
outfile.write("Sequences with =>98% similarity i.e., Hamming distance <= 0.02, are assumed to be possible contaminants\n\n")

##OPEN Gdrive and search for the appropriate gene files and calculate hamming distance.. output possible contaminats in above outfile

gdrive = os.environ.get("ALL_SEQUENCES_PATH") + "/"
wholegenome = os.environ.get("WHOLEGENOMES_PATH") + "/"

###first get data from whole genome
wlist = os.walk(wholegenome)
wgfiles = []			#paths to all fasta files in whole_genomes folder
for wl in wlist:
	wfolder = wl[0]+"/"
	wfiles = wl[-1]
	for wfile in wfiles:
		if wfile.endswith(".fasta") or wfile.endswith(".fas"):
			wpath = wfolder+wfile
			w = open(wpath,"r")
			wlines = w.readlines()
			if not wlines[0].strip().startswith(">"):
				w.close()
			elif wlines[0].strip().startswith(">"):
				wgfiles.append(wpath)

#if folder.endswith(genename)... then open
if os.path.isdir(gdrive+genename):
	rrunfiles = []                  #list of fastafiles to analyse
	gfolderlist = os.listdir(gdrive+genename)
	for gfolder in gfolderlist:	#all subfolders for the specified gene (person names)
		gfpath = gdrive+genename+"/"+gfolder
#		gflist = os.listdir(gfpath)
		gflist = os.walk(gfpath)
		for gfile in gflist:
			fflder = gfile[0]+"/"		#subfolder name
			ffiles = gfile[-1]		#list of files/subfolders in this subfolder
		############get fasta file only.. include a check of infole formating
			for ffile in ffiles: 
				if ffile.endswith(".fasta") or ffile.endswith(".fas"):	
#					print ffile, fflder+ffile

					fpath = fflder+ffile
					f = open(fpath,"r")
					#first check if actually in fasta 
					chlines = f.readlines()
					if not chlines[0].strip().startswith(">"):
						f.close()
					elif chlines[0].strip().startswith(">"):
						if fpath not in rrunfiles:
							rrunfiles.append(fpath)

	####also include the whole_genome files########
	runfiles = rrunfiles+ wgfiles
		
	#######compare data and calculate hamming distance##########
	for rfile in runfiles:
#		print rfile
		fff = open(rfile,"r")
		records = fastarecords(fff)
		oldtitles = []				#lsit of seqnames
		for title in records:
			oldtitles.append(title)
	#		print title
			seq = records[title]
			seq2 = seq.strip("-")	#remove trailing gaps both ends
			seq2 = seq2.strip("~")
			for newt in newrecords:
				newseq = newrecords[newt].strip("-")
				newseq = newseq.strip("~")
				clustaltempfile = str(oldtitles.index(title)+1)+"_"+str(titlelist.index(newt)+1)+"clustaltemp.fas"
				clustf = open(clustaltempfile,"w")
				clustf.write(">"+newt+str(titlelist.index(newt)+1)+"\n"+newseq+"\n")	#added index to avoid same names with old data
				clustf.write(">"+title+"\n"+seq2+"\n")
				clustf.close()
		#run clustalw on created tempfile
				clustoutname = clustaltempfile[:-3]+"aln"
				clustl = "clustalw2 "+clustaltempfile+" -outfile="+clustoutname
				os.system(clustl)
		#open clustaw output file and get the seqs
		
				pairs = open(clustoutname,"r")
				prs = pairs.readlines()
				seqn1 = prs[3].strip().split(" ")[0].strip()	#get newseq name
				seqn2 = prs[4].strip().split(" ")[0].strip()	#get oldseq name
				seqd1raw = ""
				for i in range(3,len(prs)-2,4):
					seqdata1 = prs[i].strip().split(" ")[-1].strip()
					seqd1raw = seqd1raw+seqdata1
				#strip trailing gaps from aligning process, i.e. in the new seq if shorter than existing, genome cases
				seqd1 = seqd1raw.strip("-")		
				sspan = re.search(seqd1,seqd1raw).span()

				seqd2raw = ""
				for j in range(4,len(prs)-1,4):
					seqdata2 = prs[j].strip().split(" ")[-1].strip() 
					seqd2raw = seqd2raw+seqdata2
				seqd2 = seqd2raw[sspan[0]:sspan[-1]]
				
				
#				#get hamming distance and print those with distance <= 30% i.e with >70 similarity	
				diff = hamdist(seqd1, seqd2)
				if diff <= 0.02:
					hammings.append(diff)
#					print "Possible contamination: Hamming distance (difference) = ", hamdist(seqd1, seqd2)
					outfile.write("\nThe new sequence "+newt+" has hamming distance = "+str(diff)+" from sequence "+title+" in "+rfile+"\n")
#					outfile.write(seqd1+"\n"+seqd2+"\n\n")
#				#delete the clustal temp input and output files
				rmcommand = "rm "+clustaltempfile+" "+clustoutname+" "+clustaltempfile[:-3]+"dnd"
				os.system(rmcommand)
				
##if no possible contminats
if len(hammings)==0:
	outfile.write("\nNO POSSIBLE CONTAMINANTS FOUND FOR YOUR NEW SEQUENCE DATA")
else:
	pass
outfile.close()

