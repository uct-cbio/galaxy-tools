#!/usr/bin/python

# Predict metagenomic genes using metagenemark or ...

# Need to sort out from where galaxy user reads metagenemark key file! Currently in /home/gerrit/.gm_key how it does that not sure! Consult metagenmark documentation.

import sys, re, string, os, commands, glob
from optparse import OptionParser
from Bio import SeqIO

def main():
  
    usage = "usage: -s SEQ_FILE -p PREDICT_ALGORITHM -t PREDICT_TYPE -o PREDICT_SEQ_FILE -a PREDICT_SUMMARY_FILE -g PREDICT_GFF_FILE -w PROG_WORK_DIR -d"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--seq", dest="seq_file", help="Nucleotide FASTA sequences")
    parser.add_option("-p", "--predict_algorithm", dest="predict_algorithm", help="The algorithm/software to use to predict the peptide sequences (currently metagenemark only).")
    parser.add_option("-t", "--predict_type", dest="predict_type", help="The type of output. Predicted protein sequences (p) or nucleotide sequences (n).")
    parser.add_option("-o", "--predict_seq", dest="predict_seq_file", help="The predicted FASTA nucleotide or protein sequences.")
    parser.add_option("-a", "--predict_summary", dest="predict_summary_file", help="A list of assembly names and their member genes")
    parser.add_option("-g", "--predict_gff", dest="predict_gff_file", help="If specified a GFF file of the predicted genes will be generated")
    parser.add_option("-w", "--prog_work", dest="prog_work_dir", help="Program working directory, contains all processed files.")
    parser.add_option("-d", "--delete_program_work", action="store_true", dest="delete_program_work", default=False, help="Delete program working directory after program has completed.")
    (options, args) = parser.parse_args()
    
    if not options.seq_file:
        print "Please specify the nucleotide FASTA sequence file (-s SEQ_FILE)"
        return - 1
    if not options.predict_algorithm:
        print "Please specify the algorithm/software to use to predict the peptide sequences (currently metagenemark only). (-p PREDICT_ALGORITHM)"
        return - 3
    if not options.predict_type:
        print "Please specify the type of ouput.  p = protein sequences or n =  nucleotide sequences. (-t PREDICT_TYPE)"
        return - 3
    if not options.predict_seq_file:
        print "Please specify the predicted sequence output file (-o PREDICT_SEQ_FILE)"
        return - 4
    if not options.predict_summary_file:
        print "Please specify the predicted summary output file (-a PREDICT_SUMMARY_FILE)"
        return - 5
    if not options.prog_work_dir:
        print "Please specify the program working directory (-w PROG_WORK_DIR)"
        return - 6
    if (len(args) > 0):
        print "Too many input arguments"
        return - 7
    
    # Do some initialization
    print "Initialize..." 
    root_dir = os.getcwd()
    
    # Get full file paths
    seq_file = os.path.abspath(options.seq_file)
    predict_seq_file = os.path.abspath(options.predict_seq_file)
    predict_summary_file = os.path.abspath(options.predict_summary_file)
    prog_work_dir = os.path.abspath(options.prog_work_dir)
    
    if options.predict_gff_file:
            predict_gff_file = os.path.abspath(options.predict_gff_file)
    else:
        predict_gff_file = ""
    
    timestamp = commands.getoutput("date +%Y-%m-%d_%H_%M_%S_%N")
    base_dir = prog_work_dir + "/predict_peptide_seq_metagenomic_" + timestamp 

    if not os.path.isdir(prog_work_dir):
        os.system("mkdir " + prog_work_dir)
    if os.path.isdir(prog_work_dir):
        if os.path.isdir(base_dir):
            os.system("rm -rf " + base_dir)
            os.system("mkdir " + base_dir)
        else:
            os.system("mkdir " + base_dir);
    else:
       print "Program working directory does not exist."
       return - 8
            
    predict_algorithm_dir = base_dir + "/" + options.predict_algorithm
    os.mkdir(predict_algorithm_dir)
    
    predict_algorithm_options = ""
    
    # Run predict algorithm
    if (options.predict_algorithm == "metagenemark"): # phrap
        print "Run metagenemark..."
        metagenemark_model = os.environ.get("METAGENEMARK_MODEL_FILE")
        os.chdir(predict_algorithm_dir)
        metagenemark_out_file = predict_algorithm_dir + "/" + "metagenemark.txt"
        
        if not((options.predict_type == "p") or (options.predict_type == "n")):
            print "MetaGeneMark: please specify protein(p) or nucleotide(n) for the predicted output type"
            return -10
        
        if(options.predict_type == "p"):
            predict_algorithm_options = "-a"
        if(options.predict_type == "n"):
            predict_algorithm_options = "-d"
        run_metagenemark(predict_algorithm_dir, metagenemark_model, predict_algorithm_options, seq_file, metagenemark_out_file)
        parse_metagenemark(predict_algorithm_dir, metagenemark_out_file, predict_summary_file)
        os.system("cat " + predict_algorithm_dir + "/predicted_seq/*.fsa > " + predict_seq_file)
        
        if predict_gff_file: # Rerun prediction if user needs gff output
            predict_algorithm_options = "-f g"
            run_metagenemark(predict_algorithm_dir, metagenemark_model, predict_algorithm_options, seq_file, predict_gff_file)
                    
    else:
        print "Predict algorithm\/software not supported (currently metagenemark only)"
        return - 9
    
    # Prepare output
    print "Prepare output..."

    # Delete program working directory if indicated
    if(options.delete_program_work):
        print "Delete working directory"
        os.system("rm -rf " + base_dir)
        
    # Done
    print "Done."
        
    return 0

def run_metagenemark(metagenemark_dir, model, algorithm_options, seq_file, out_file):
    log_file = metagenemark_dir + "/metagenemark.log"
    os.system("gmhmmp -m " + model + " " + algorithm_options + " -o " + out_file + " " + seq_file + " >> " + log_file + " 2>& 1")
    
def parse_metagenemark(metagenemark_dir, metagenemark_file, assembly_summary_file):
    fd_in = open(metagenemark_file, "r")
    fd_as = open(assembly_summary_file, "w")
    predicted_seq_dir = metagenemark_dir + "/predicted_seq"
    if not os.path.isdir(predicted_seq_dir):
        os.system("mkdir " + predicted_seq_dir)
    start_seq = False
    seq = ""
    for line in fd_in:
        if line.startswith(">"):
            start_seq = True
        if start_seq and not line.startswith("\n"):
            seq =  seq + line
        else:
            start_seq = False
    fd_out = open("tmp.fsa", "w")
    fd_out.write(seq)
    fd_out.close()
    split_fasta_metagenemark("tmp.fsa", predicted_seq_dir, fd_as)
#    os.system("rm -rf tmp.fsa")
    fd_in.close()
    fd_as.close()
    
def split_fasta_metagenemark(seq_file, split_seq_dir, fd_as):
    
    assembly_map = {}
    fd_seq = open(seq_file, "rU")
    seqs = SeqIO.parse(fd_seq, "fasta")

    for seq in seqs:
        # Get name of original sequence from metagenemark header
        # Split on all whitespaces, returns the first string with " > " removed
        origin_seq_name = seq.description[seq.description.find(">") + 1:len(seq.description)].split()[0] 
        # Get name of gene sequence from metagenemark header
        # Split on " | ", returns the first string
        gene_seq_name = seq.name.split("|")[0]
        gene_seq_name = origin_seq_name + "_" + gene_seq_name
        gene_seq_file = split_seq_dir + "/" + gene_seq_name + ".fsa"
        fd_split_seq = open(gene_seq_file, "w")
        # remove the next 2 lines if you would like to keep the original metagenemark header
        seq.description = seq.id
        seq.id = gene_seq_name
        #########################        
        SeqIO.write([seq], fd_split_seq, "fasta")
        fd_split_seq.close()
        #os.system("sed - i \"s/> />/g\" " + new_seq_file)
        
        if(origin_seq_name not in assembly_map):
            genes = []
            genes.append(gene_seq_name)
            assembly_map[origin_seq_name] = genes
        else:
            genes = assembly_map[origin_seq_name]
            genes.append(gene_seq_name)
            assembly_map[origin_seq_name] = genes

    fd_seq.close()

    for origin_seq_name in assembly_map:
        fd_as.write(origin_seq_name)
        genes = assembly_map[origin_seq_name]
        for gene_seq_name in genes:
            fd_as.write("\t" + gene_seq_name)
        fd_as.write("\n")

if __name__ == "__main__":
#    split_qual_cap3('XHC00016.fsa.cap.contigs.qual', 'tmp')
    sys.exit(main())


