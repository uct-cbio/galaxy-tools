#!/usr/bin/python

# Select a FASTA sequences from a input list. The name in the list should match the header exactly. Except if the -g flag is specified, then the program
# will look for a GI number match in the header.
# Input:
#    1) file with fasta sequences
#    2) list with names of fasta seqeuences
# Output
#    1) file with the selected sequences 

import sys
import re
import string
from optparse import OptionParser
import os
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    usage = "usage: %prog -i INPUT_FILE -l LIST_FILE -o OUTPUT_FILE -g"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--input", dest="input_file", help="Input FASTA file.")
    parser.add_option("-l", "--list", dest="list_file", help="List of sequence names or GI numbers.")
    parser.add_option("-o", "--output", dest="output_file", help="FASTA output file, containing the sequences with GI numbers in the input list.")
    parser.add_option("-g", "--gi_match", action="store_true", dest="gi_match", default=False, help="Look for GI number matches instead of exact header matches.")
 
    (options, args) = parser.parse_args()
    
    if not options.input_file:
        print "Please specify the FASTA input file (-i INPUT_FILE)"
        return - 1
    if not options.list_file:
        print "Please specify the input file with the list of GI numbers (-l LIST_FILE)"
        return - 2
    if not options.output_file:
        print "Please specify the FASTA output file (-o OUTPUT_FILE)"
        return - 3
    if (len(args) > 0):
        print "Too many arguments"
        return - 4

    input_file = options.input_file
    list_file_csv = csv.reader(open(options.list_file), delimiter='\t') 
    output_file = options.output_file
    fd_out = open(output_file, "w")
    
    if(not options.gi_match):
        fd_in = open(input_file , "r");
        seqs = SeqIO.to_dict(SeqIO.parse(fd_in, "fasta")) # later conider using SeqIO.index in biopython 1.52> . 
        seqs_selected = {}
        for row in list_file_csv:
            seq_name = row[0]
            if not seq_name in seqs_selected: # do not write duplicates
                if(seq_name in seqs):
                    SeqIO.write([seqs[seq_name]], fd_out, "fasta")
                    seqs_selected[seq_name] = seq_name
            else:
                print "Duplicate: ", seq_name
        fd_in.close()
    else:
        gis = set()
        for row in list_file_csv:
            gi = row[0]
            if gi not in gis:
                gis.add(gi)
        for record in SeqIO.parse(input_file, "fasta"):
            m = re.search(r'gi\|(\d+)\|', record.id)
            if(m): # check if their was a match
                if(m.lastindex == 1): # just checking if the gi number was found
                    gi = m.group(1)
                if gi in gis:
                    SeqIO.write([record], fd_out, "fasta")
                    gis.remove(gi)
            
    fd_out.close()

if __name__ == "__main__":
    sys.exit(main())


