#!/usr/bin/python

# Convert a GFF with a corresponding FASTA to GenBank Format. No Validation of input files are done. 
# E.g. checking if there is a fasta sequence for a genbank entry.
# gff to genbank converter biopython

# BCBio library downloaded from https://github.com/chapmanb/bcbb/tree/master/gff

import sys
import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from BCBio import GFF
from optparse import OptionParser

def main():
    usage = "usage: %prog -g gff -f fasta -o genbank"
    parser = OptionParser(usage=usage)
    parser.add_option("-g", "--gff", dest="gff_file", help="GFF file.")
    parser.add_option("-f", "--fasta", dest="fasta_file", help="FASTA file")
    parser.add_option("-o", "--genbank", dest="genbank_file", help="GenBank output file")
    (options, args) = parser.parse_args()
    
    if not options.gff_file:
        print "Please specify the GFF file (-g gff)"
        return - 1
    if not options.fasta_file:
        print "Please specify the FASTA file (-f fasta)"
        return - 2
    if not options.genbank_file:
        print "Please specify the GenBank output file (-o genbank)"
        return - 3
    if (len(args) > 0):
        print "Too many input arguments"
        return - 4
    
    gff_file = os.path.abspath(options.gff_file)
    fasta_file = os.path.abspath(options.fasta_file)
    genbank_file = os.path.abspath(options.genbank_file)
    if os.path.isfile(genbank_file):
        os.system("rm -rf " + genbank_file)
    
    fasta_input = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta", generic_dna))
    gff_iter = GFF.parse(gff_file, fasta_input)
    SeqIO.write(gff_iter, genbank_file, "genbank")

if __name__ == "__main__":
    sys.exit(main())
