#!/usr/bin/python

# EST clustering

# Currently clustering are done using the wcd clustering algorithm. Other algorithms will later be supported.

# The wcd program does not use qual score for clustering. If seq qual scores are specified the 
# wcd clustered information will be used to create cluster qual scores. This quality scores can then be used as
# input for further processing (e.g. assembly).

import sys, re, string, os, subprocess, commands
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import UnknownSeq

def main():
    usage = "usage: %prog -s SEQ_FILE -S QUAL_FILE -c CLUSTER_SEQ_FILE -C CLUSTER_QUAL_FILE -i CLUSTER_ID -r CLUSTER_SUMMARY_FILE -p CLUSTER_ALGORITHM -w PROG_WORK_DIR -d"
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--seq", dest="seq_file", help="File with sequences to be clustered.")
    parser.add_option("-S", "--qual", dest="qual_file", help="Quality scores for the input sequences (optional).")
    parser.add_option("-c", "--cluster_seq", dest="cluster_seq_file", help="Zipped FASTA sequences file. Archive contains cluster file with corresponding clustered FASTA sequence.")
    parser.add_option("-C", "--cluster_qual", dest="cluster_qual_file", help="Zipped QUAL scores file. Archive contains cluster files with corresponding clustered QUAL scores (optional).")
    parser.add_option("-i", "--cluster_id", dest="cluster_id", help="Cluster id")
    parser.add_option("-r", "--cluster_summary", dest="cluster_summary_file", help="Cluster summary report.")
    parser.add_option("-p", "--cluster_algorithm", dest="cluster_algorithm", help="The cluster algorithm (currently wcd only).")
    parser.add_option("-w", "--prog_work", dest="prog_work_dir", help="Program working directory, contains all processed files.")
    parser.add_option("-d", "--delete_program_work", action="store_true", dest="delete_program_work", default=False, help="Delete program working directory after program has completed.")
    (options, args) = parser.parse_args()
    
    if not options.seq_file:
        print "Please specify the FASTA sequence file (-s SEQ_FILE)"
        return - 1
    if not options.cluster_seq_file:
        print "Please specify the zipped cluster FASTA sequence file (-c CLUSTER_SEQ_FILE)"
        return - 3
    if not options.cluster_id:
        print "Please specify the cluster id (-i CLUSTER_ID)"
        return - 5
    if not options.cluster_summary_file:
        print "Please specify the cluster summary report file (-r CLUSTER_SUMMARY_FILE)"
        return - 6
    if not options.cluster_algorithm:
        print "Please specify the cluster algorithm (-p CLUSTER_ALGORITHM)"
        return - 7
    if not options.prog_work_dir:
        print "Please specify the program working directory (-w PROG_WORK_DIR)"
        return - 8
    if (len(args) > 0):
        print "Too many input arguments"
        return - 9
    
    # Do some initialization
    print "Initialize..." 
    root_dir = os.getcwd()
    cluster_id = options.cluster_id
    
    # Get full file paths
    seq_file = os.path.abspath(options.seq_file)
    if options.qual_file:
        qual_file = os.path.abspath(options.qual_file)
    cluster_seq_file = os.path.abspath(options.cluster_seq_file)
    if options.cluster_qual_file:
        cluster_qual_file = os.path.abspath(options.cluster_qual_file)
    cluster_summary_file = os.path.abspath(options.cluster_summary_file)
    prog_work_dir = os.path.abspath(options.prog_work_dir)
    
        
    timestamp = commands.getoutput("date +%Y-%m-%d_%H_%M_%S_%N")
    base_dir = prog_work_dir + "/cluster_" + timestamp 

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
       return - 10
    
    # Create working directory 
    if os.path.isdir(base_dir):
        os.system("rm -rf " + base_dir)
        os.system("mkdir " + base_dir)
    else:
        os.system("mkdir " + base_dir)
        
    # Run wcd clustering
    if (options.cluster_algorithm == "wcd"):
        print "Run wcd..." 
        cluster_seq_dir = run_wcd(base_dir, seq_file, cluster_id)
    else:
        print "Cluster algorithm not supported (currently wcd only)"
        return - 10
    
    # Prepare cluster qual files
    if(options.qual_file and options.cluster_qual_file):
        print "Prepare cluster quality scores..."        
        cluster_qual_dir = prepare_cluster_qual_files(base_dir, qual_file, cluster_seq_dir)
    elif (not options.qual_file and options.cluster_qual_file):
        print "No cluster quality score will be prepared"
        print "Please specify the QUAL score file (-S QUAL_FILE)"
    elif (options.qual_file and not options.cluster_qual_file):
        print "No cluster quality score will be prepared"
        print "Please specify the zipped clustered QUAL scores file (-C CLUSTER_QUAL_FILE)"
    else:
       print "No cluster quality score will be prepared"
       print "Please specify the QUAL score file (-S QUAL_FILE)"
       print "Please specify the zipped clustered QUAL scores file(-C CLUSTER_QUAL_FILE)"
       
    # Write summary report
    print "Write summary report..."
    write_cluster_summary(cluster_seq_dir, cluster_summary_file)
    
    # Prepare output
    print "Prepare output..."
    tmp_zip = base_dir + "/tmp.zip" # Galaxy work around need to create a temporary zip archive and move to the output data set
    
    if(options.qual_file and options.cluster_qual_file):
#        os.system ("zip -j " + cluster_qual_file + " " + cluster_qual_dir + "/*")
        os.system ("zip -j " + tmp_zip + " " + cluster_qual_dir + "/*")
        os.system("mv " + tmp_zip + " " + cluster_qual_file)
        
#    os.system("zip -j " + cluster_seq_file + " " + cluster_seq_dir + "/*")
    os.system("zip -j " + tmp_zip + " " + cluster_seq_dir + "/*")
    os.system("mv " + tmp_zip + " " + cluster_seq_file)
    
    # Delete program working directory if indicated
    if(options.delete_program_work):
        print "Delete working directory"
        os.system("rm -rf " + base_dir)
        
    # Done
    print "Done."
    return 0

def run_wcd(work_dir, seq_file, cluster_id):
    cluster_seq_dir = work_dir + "/cluster_fasta"
    os.mkdir(cluster_seq_dir)
    os.chdir(work_dir)
    os.system("wcd --show_clusters -o cluster.cls " + seq_file) # get clusters
    os.chdir(cluster_seq_dir)
    os.system("wcd --init_cluster ../cluster.cls --split " + cluster_id + " " + seq_file) # prepare cluster FASTA files
    os.chdir(work_dir)
    return cluster_seq_dir
    
def prepare_cluster_qual_files(work_dir, qual_file, cluster_seq_dir):
    cluster_qual_dir = work_dir + "/cluster_qual"
    os.mkdir(cluster_qual_dir)
    # get a list of all quality scores
    fd_qual = open(qual_file, "rU");
    quals = SeqIO.to_dict(SeqIO.parse(fd_qual, "qual"));
    # get quality scores for the clusters
    for cluster_seq_file in os.listdir(cluster_seq_dir):
        if os.path.isfile(cluster_seq_dir + "/" + cluster_seq_file): # check if file, can do some more checking here e.g. is fasta file
            fd_cluster_seq = open(cluster_seq_dir + "/" + cluster_seq_file, "rU")
            cluster_seqs = SeqIO.parse(fd_cluster_seq, "fasta")
            cluster_quals = []
            for seq in cluster_seqs:
                qual = quals[seq.name]
                cluster_qual = SeqRecord(seq=UnknownSeq(len(qual.letter_annotations["phred_quality"])), id="", description=qual.description)
                cluster_qual.letter_annotations["phred_quality"] = qual.letter_annotations["phred_quality"]
                cluster_quals.append(cluster_qual)
                
            cluster_qual_file = cluster_qual_dir + "/" + cluster_seq_file.split(".")[0] + ".qual"        
            fd_cluster_qual = open(cluster_qual_file, "w")
            SeqIO.write(cluster_quals, fd_cluster_qual, "qual")
            fd_cluster_qual.close()
            os.system("sed -i \"s/> />/g\" " + cluster_qual_file)  # need to replace the space after the > in header
            fd_cluster_seq.close()
            
    fd_qual.close()
    return cluster_qual_dir

def write_cluster_summary(cluster_seq_dir, summary_file):
    fd_summary = open(summary_file, 'w')
    summary = ""
    for cluster_seq_file in sorted(os.listdir(cluster_seq_dir)):
        if os.path.isfile(cluster_seq_dir + "/" + cluster_seq_file): # check if file, can do some more checking here e.g. is fasta file
            fd_cluster_seq = open(cluster_seq_dir + "/" + cluster_seq_file, "rU")
            cluster_seqs = SeqIO.parse(fd_cluster_seq, "fasta")
            summary = summary + cluster_seq_file.split(".")[0]
            for seq in cluster_seqs:
                summary = summary + "\t" + seq.name
            summary = summary + "\n"
            fd_cluster_seq.close()
    process = subprocess.Popen("grep -v \"^[[:digit:]]*.$\" cluster.cls | wc -l", stdout=subprocess.PIPE, shell=True)
    nr_clusters = process.communicate()[0]
    process = subprocess.Popen("grep \"^[[:digit:]]*.$\" cluster.cls | wc -l", stdout=subprocess.PIPE, shell=True)
    nr_singletons = process.communicate()[0]
    header = "# nr Clusters: " + nr_clusters
    header = header + "# nr Singletons: " + nr_singletons
    header = header + "# Column 1: Cluster_id/Singleton_id" + "\n"
    header = header + "# Columns 2 to n: Member Sequences" + "\n"
    summary = header + summary
    fd_summary.write(summary)
    fd_summary.close()
    
if __name__ == "__main__":
    sys.exit(main())


