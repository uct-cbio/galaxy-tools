#!/usr/bin/python
#uzip an archive and combine all file contents into one file
#to run eg: python unzip_combine.py jinnynew.zip jinnynew_filelist
import sys,os

infile, outf =  sys.argv[-2], sys.argv[-1]
cmd = "unzip "+infile +" -d "+infile[:-4]+"_filelist" 
os.system(cmd)

ndir = infile[:-4]+"_filelist"
st = os.listdir(ndir)
for sst in st:
	flist = os.walk(ndir+"/"+sst)
	for fl in flist:
		for i in range(len(fl[-1])):
			fpath = fl[0]+"/"+fl[-1][i]
			cmd = "cat "+fpath+" >> "+infile[:-4]+"_ALLsequences.fas"	
			os.system(cmd)

rmc = "rm -R "+infile[:-4]+"_filelist"
os.system(rmc)
