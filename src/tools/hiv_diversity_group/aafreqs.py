#!/usr/bin/python

#get and compare AA frequencies bwtn 2 alignments using fisher's exact test
#to run:  /opt/rocks/bin/python aafreqs.py aafreqsSensitive.txt aafreqsResistant.txt

import sys
import rpy
seqfile1, seqfile2 = sys.argv[-2], sys.argv[-1]			#tab_del aa alignment files

def aafrequency(alignment1, alignment2):			#for each alignment
	f1 = open(alignment1,"r")
	lines1 = f1.readlines()
	sub1 = []				#list of file1 seqs
	allresults = {}
	for line1 in lines1:
		lin = line1.strip().split("\t")
		ID = lin[0].strip()
		seq = lin[-1].strip()
		sub1.append(seq)
	f2 = open(alignment2,"r")
        lines2 = f2.readlines()
        sub2 = []                               #list of file2 seqs
        for line2 in lines2:
                lin2 = line2.strip().split("\t")
                ID2 = lin2[0].strip()
                seq2 = lin2[-1].strip()
                sub2.append(seq2)

	for i in range(len(sub1[0])):			#for each pstn in alignment
		s = str(i+1)			#actual site
		allaa = []			#list of all aa at this pstn from both C and B
		for sq in sub1:
			if sq[i] not in allaa:	#add all res to allaa list
				allaa.append(sq[i])
		for sq in sub2:
			if sq[i] not in allaa:
				allaa.append(sq[i])
	#NOW GET THE aa  in each alignment at this site
		afreqs = {}	#aa:[count1,count2]	
		aa1 = []	#list of actual freqs in 1 for stats test
		aa2 = []	#list of actual freqs in 2 for stats test
		for a in allaa:	#get freq of each res 
			afreqs[a] = [0,0]
			for ssq in sub1:
				if ssq[i] == a:
					afreqs[a][0] = afreqs[a][0]+1
			for ssq in sub2:
				if ssq[i] == a:
					afreqs[a][-1] = afreqs[a][-1]+1
	#	print i+1, afreqs		
		
		for r in afreqs:
			aa1.append(afreqs[r][0])
			aa2.append(afreqs[r][-1])
	#	print aab, aac
		allresults[str(i+1)] = [afreqs,aa1,aa2]
	f1.close()
	f2.close()
	return allresults

data = aafrequency(seqfile1, seqfile2)
outf = open("aafreqs_FISHERresults.txt","w")
outf.write("Site\tPvalue\tPvalue(Bonferroni_Corrected)\n")
outf2 = open("aafreqs_SITEresidueCounts.txt","w")
outf2.write("Site\tResidue\tCountINfile1\tCountINfile2\n")
r_fisher_file = os.environ.get("R_FISHER_FILE")

from rpy import *
#for s in range(1,5):
#site = str(s+1)
for site in data:
#	print site, data[site]
	for aa in data[site][0]:		#get res count per site into outf2
		outf2.write(site+"\t"+aa+"\t"+str(data[site][0][aa][0])+"\t"+str(data[site][0][aa][-1])+"\n")

	if len(data[site][1])>1:		#do fisher's exact test on aa counts
		r.source(r_fisher_file)
	#		print data[site][1],data[site][-1]
		fisher = r.my_fisher(data[site][1],data[site][-1])
#		if fisher['p.value'] < 0.045:
#			print site, data[site][1], data[site][-1],"fisher pvalue = : ", fisher['p.value'], fisher['p.value']*len(data)
		outf.write(site+"\t"+str(fisher['p.value'])+"\t"+str(fisher['p.value']*len(data))+"\n")
	else:
		outf.write(site+"\t"+"data_small"+"\n")
outf.close()
outf2.close()