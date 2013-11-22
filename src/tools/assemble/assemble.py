#!/usr/bin/python

# Assembles clusters with phrap or cap3

import sys, re, string, os, commands, glob
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import UnknownSeq
from Bio.Sequencing import Ace

def main():
    
    PHRAP_DEFAULT_OPTIONS = "-penalty -2 -gap_init -2 -gap_ext -1 -ins_gap_ext -3 -del_gap_ext -3 -vector_bound 0 -trim_score 20 -forcelevel 0 -retain_duplicates -new_ace"
    CAP3_DEFAULT_OPTIONS = " "

    usage = "usage: -c CLUSTER_SEQ_FILE -C CLUSTER_QUAL_FILE -a CONTIG_SEQ_FILE -A CONTIG_QUAL_FILE -e CONTIG_ACE_FILE -m CONTIG_IMAGE_FILE \
    -s SINGLETON_SEQ_FILE -S SINGLETON_QUAL_FILE -p ASSEMBLY_ALGORITHM -r ASSEMBLY_SUMMARY_FILE_A -b ASSEMBLY_SUMMARY_FILE_B -w PROG_WORK_DIR -d"
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--cluster_seq", dest="cluster_seq_file", help="Zipped FASTA sequencesfile. Archive contains cluster files with corresponding clustered FASTA sequence.")
    parser.add_option("-C", "--cluster_qual", dest="cluster_qual_file", help="Zipped QUAL scores file. Archive contains cluster files with corresponding clustered QUAL scores (optional).")
    parser.add_option("-a", "--contig_seq", dest="contig_seq_file", help="Contig FASTA sequences file.")
    parser.add_option("-A", "--contig_qual", dest="contig_qual_file", help="Contig QUAL scores file") 
    parser.add_option("-e", "--contig_ace", dest="contig_ace_file", help="Zipped ACE values file. Archive contains contig files with corresponding ACE values.")
    parser.add_option("-m", "--contig_image", dest="contig_image_file", help="Zipped PNG file. Archive contains assembled images of the contigs.")
    parser.add_option("-s", "--singleton_seq", dest="singleton_seq_file", help="Singleton FASTA sequences file (optional).")
    parser.add_option("-S", "--singleton_qual", dest="singleton_qual_file", help="Singleton QUAL scores file).")
    parser.add_option("-r", "--assembly_summary_a", dest="assembly_summary_file_a", help="Assembly summary report a. All contigs and their sequences")
    parser.add_option("-b", "--assembly_summary_b", dest="assembly_summary_file_b", help="Assembly summary report b. All cluster with more then one contig and their sequences.")
    parser.add_option("-p", "--assembly_algorithm", dest="assembly_algorithm", help="The assembly algorithm (currently phrap an cap3 only).")
    parser.add_option("-o", "--assembly_algorithm_options", dest="assembly_algorithm_options", help="Command line options for the assembly algorithm. If not specified default values will be used.\
    " + "PHRAP_DEFAULT_OPTIONS=\"" + PHRAP_DEFAULT_OPTIONS + "\" CAP3_DEFAULT_OPTIONS=\"" + CAP3_DEFAULT_OPTIONS + "\"")
    parser.add_option("-w", "--prog_work", dest="prog_work_dir", help="Program working directory, contains all processed files.")
    parser.add_option("-d", "--delete_program_work", action="store_true", dest="delete_program_work", default=False, help="Delete program working directory after program has completed.")
    (options, args) = parser.parse_args()
    
    if not options.cluster_seq_file:
        print "Please specify the zipped cluster FASTA sequence file (-c CLUSTER_SEQ_FILE)"
        return - 1
#    if not options.cluster_qual_file:
#        print "Please specify the zipped cluster QUAL score file (-C CLUSTER_QUAL_FILE)"
#        return - 2
    if not options.contig_seq_file:
        print "Please specify the contig FASTA sequence file (-a CONTIG_SEQ_FILE)"
        return - 3
    if not options.contig_qual_file:
        print "Please specify the contig QUAL score file (-A CONTIG_QUAL_FILE)"
        return - 4
    if not options.contig_ace_file:
        print "Please specify the zipped contig ACE values file (-e CONTIG_ACE_FILE)"
        return - 5
    if not options.contig_image_file:
        print "Please specify the zipped PNG file (-m CONTIG_IMAGE_FILE)"
        return - 6
    if not options.singleton_seq_file:
        print "Please specify the singleton FASTA sequence file (-s SINGLETON_SEQ_FILE)"
        return - 7
#    if not options.singleton_qual_file:
#        print "Please specify the singleton contig QUAL score file (-S SINGLETON_QUAL_FILE)"
#        return - 8
    if not options.assembly_summary_file_a:
        print "Please specify the assembly summary report file a (-r ASSEMBLY_SUMMARY_FILE_A)"
        return - 9
    if not options.assembly_summary_file_b:
        print "Please specify the assembly summary report file b (-b ASSEMBLY_SUMMARY_FILE_B)"
        return - 10
    if not options.assembly_algorithm:
        print "Please specify the assembly algorithm (-p ASSEMBLY_ALGORITHM)"
        return - 11
    if not options.prog_work_dir:
        print "Please specify the program working directory (-w PROG_WORK_DIR)"
        return - 12
    if (len(args) > 0):
        print "Too many input arguments"
        return - 13
    
    # Do some initialization
    print "Initialize..." 
    root_dir = os.getcwd()
    
    # Get full file paths
    cluster_seq_file = os.path.abspath(options.cluster_seq_file)
    if options.cluster_qual_file:
        cluster_qual_file = os.path.abspath(options.cluster_qual_file)
        singleton_qual_file = os.path.abspath(options.singleton_qual_file)
    contig_seq_file = os.path.abspath(options.contig_seq_file)
    contig_qual_file = os.path.abspath(options.contig_qual_file)
    contig_ace_file = os.path.abspath(options.contig_ace_file)
    contig_image_file = os.path.abspath(options.contig_image_file)
    singleton_seq_file = os.path.abspath(options.singleton_seq_file)    
    assembly_summary_file_a = os.path.abspath(options.assembly_summary_file_a)
    assembly_summary_file_b = os.path.abspath(options.assembly_summary_file_b)
    prog_work_dir = os.path.abspath(options.prog_work_dir)
    
    timestamp = commands.getoutput("date +%Y-%m-%d_%H_%M_%S_%N")
    base_dir = prog_work_dir + "/assemble_" + timestamp 

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
       return - 14 

    algorithm_dir = base_dir + "/" + options.assembly_algorithm
    os.mkdir(algorithm_dir)
    
    # Unzip cluster fasta
    cluster_seq_dir = base_dir + "/cluster_fasta"
    os.mkdir(cluster_seq_dir)
    os.system("unzip -qq " + cluster_seq_file + " -d " + cluster_seq_dir + " > /dev/null")
    # Force all extensions to .fsa , quick and dirty
    for seq_file in os.listdir(cluster_seq_dir):
        new_seq_file = algorithm_dir + "/" + seq_file.split(".")[0] + ".fsa"
        os.rename(cluster_seq_dir + "/" + seq_file, new_seq_file)
    os.system("rm -rf " + cluster_seq_dir)
    # If cluster qual scores was specified
    if (options.cluster_qual_file):
        # Unzip cluster qual scores
        cluster_qual_dir = base_dir + "/cluster_qual"
        os.mkdir(cluster_qual_dir)
        os.system("unzip -qq " + cluster_qual_file + " -d " + cluster_qual_dir + " > /dev/null")
        # Force all extensions to .fsa.qual (format is needed for both phrap and cap3), quick and dirty
        for qual_file in os.listdir(cluster_qual_dir):
            new_qual_file = algorithm_dir + "/" + qual_file.split(".")[0] + ".fsa.qual"
            os.rename(cluster_qual_dir + "/" + qual_file, new_qual_file)
        os.system("rm -rf " + cluster_qual_dir)   
    else:
        print "Cluster quality scores file was not specified, scores will not be used in assembly"
        print "Please specify the zipped cluster QUAL score file (-C CLUSTER_QUAL_FILE)"
        print "Continue without quality scores..."
    
    # Run assembly
    if (options.assembly_algorithm == "phrap"): # phrap
        print "Run phrap..."
        if(options.assembly_algorithm_options):
            print "User PHRAP options"
            run_phrap(algorithm_dir, options.assembly_algorithm_options)
        else:
            print "Default PHRAP options"
            run_phrap(algorithm_dir, PHRAP_DEFAULT_OPTIONS)
    elif (options.assembly_algorithm == "cap3"): # cap3
        print "Run cap3..."
        if(options.assembly_algorithm_options):
            print "User CAP3 options"
            run_cap3(algorithm_dir, options.assembly_algorithm_options)
        else:
            print "Default CAP3 options"
            run_cap3(algorithm_dir, CAP3_DEFAULT_OPTIONS)
    else:
        print "Assembly algorithm not supported (currently phrap and cap3 only)"
        return - 10
    
    # Parse contigs
    print "Parse contigs..."
    contig_seq_dir = base_dir + "/contig_fasta"
    os.mkdir(contig_seq_dir)
    contig_qual_dir = base_dir + "/contig_qual"
    os.mkdir(contig_qual_dir)
    parse_contigs(contig_seq_dir, contig_qual_dir, options.assembly_algorithm, algorithm_dir)
    
    # Parse singletons
    print "Parse singletons..."
    singleton_seq_dir = base_dir + "/singleton_fasta"
    os.mkdir(singleton_seq_dir)
    singleton_qual_dir = base_dir + "/singleton_qual"
    os.mkdir(singleton_qual_dir)
    parse_singletons_fasta(singleton_seq_dir, singleton_qual_dir, options.assembly_algorithm, algorithm_dir)
    
    # Parse ace files
    print "Parse ace files..."
    contig_ace_dir = base_dir + "/contig_ace"
    os.mkdir(contig_ace_dir)
    parse_ace(contig_ace_dir, options.assembly_algorithm, algorithm_dir)
    
    parse_singletons_fasta_in_ace(contig_ace_dir, singleton_seq_dir) # singletons need to be removed from ace file, some may sneak in
    
    if options.cluster_qual_file: # if the zipped QUAL scores were provided
        tmp_dir = base_dir + "/tmp" # can now get the quality scores for the singletons
        os.mkdir(tmp_dir)
        parse_singletons_qual(singleton_seq_dir, singleton_qual_dir, algorithm_dir, tmp_dir) 
        os.system("rm -rf " + tmp_dir)
    
    # Generate contig image
    print "Generate contig images..."
    contig_image_dir = base_dir + "/contig_image"
    os.mkdir(contig_image_dir)
    generate_contig_images(contig_image_dir, contig_ace_dir)
    
    # Write summary reports
    print "Write summary reports..."
    write_contig_summary_a(contig_ace_dir, singleton_seq_dir, assembly_summary_file_a)
    write_contig_summary_b(contig_ace_dir, singleton_seq_dir, assembly_summary_file_b)
    
    # Prepare output
    print "Prepare output..."
    tmp_zip = base_dir + "/tmp.zip" # Galaxy work around need to create a temporary zip archive and move to the output data set
 
    # contigs
    os.system ("cat " + contig_seq_dir + "/*.fsa > " + contig_seq_file)
    os.system ("cat " + contig_qual_dir + "/*.qual > " + contig_qual_file)
    # singletons
    if (len(glob.glob(singleton_seq_dir + "/*.fsa")) > 0):
        os.system ("cat " + singleton_seq_dir + "/*.fsa > " + singleton_seq_file)
    else:
        os.system ("touch " + singleton_seq_file)
    if options.cluster_qual_file: # if the zipped QUAL scores were provided
        if (len(glob.glob(singleton_qual_dir + "/*.qual")) > 0):        
             os.system ("cat " + singleton_qual_dir + "/*.qual > " + singleton_qual_file)
        else:
            os.system ("touch " + singleton_qual_file)
        
#    os.system("zip -j " + contig_ace_file + " " + contig_ace_dir + "/*.ace")
    os.system("zip -qqj " + tmp_zip + " " + contig_ace_dir + "/*.ace")
    os.system("mv " + tmp_zip + " " + contig_ace_file)
#    os.system("zip -j " + contig_image_file + " " + contig_image_dir + "/*.png")
    os.system("zip -qqj " + tmp_zip + " " + contig_image_dir + "/*.png")
    os.system("mv " + tmp_zip + " " + contig_image_file)
    
    # Delete program working directory if indicated
    if(options.delete_program_work):
        print "Delete working directory"
        os.system("rm -rf " + base_dir)
        
    # Done
    print "Done."
        
    return 0

def run_phrap(phrap_dir, phrap_options):
    log_file = phrap_dir + "/phrap.log"
    os.chdir(phrap_dir)
    for seq_file in glob.glob("*.fsa"):
        os.system("phrap " + phrap_options + " " + seq_file + " >> " + log_file + " 2>&1")

def run_cap3(cap3_dir, cap3_options):
    log_file = cap3_dir + "/cap3.log"
    os.chdir(cap3_dir)
    for seq_file in glob.glob("*.fsa"):
        os.system("cap3 " + cap3_options + " " + seq_file + " >> " + log_file + " 2>&1")
        
def parse_contigs(contig_seq_dir, contig_qual_dir, algorithm, algorithm_dir):
    os.chdir(algorithm_dir)
    # cap3
    if (algorithm == "cap3"):
        for seq_file in glob.glob("*.fsa.cap.contigs"):
            cluster_name = seq_file.split(".")[0]
            if(os.path.getsize(algorithm_dir + "/" + seq_file) > 0): # if contigs exists
                os.system("sed -i \"s/.*Contig/>" + cluster_name + "_/g\" " + algorithm_dir + "/" + seq_file)
                split_fasta(algorithm_dir + "/" + seq_file, contig_seq_dir)
                os.system("sed -i \"s/.*Contig/>" + cluster_name + "_/g\" " + algorithm_dir + "/" + seq_file + ".qual")
                split_qual_cap3(algorithm_dir + "/" + seq_file + ".qual", contig_qual_dir) # needs a different qual split
    # phrap            
    elif (algorithm == "phrap"):
        for seq_file in glob.glob("*.fsa.contigs"):
            cluster_name = seq_file.split(".")[0]
            if(os.path.getsize(algorithm_dir + "/" + seq_file) > 0): # if contigs exists
                os.system("sed -i \"s/.*Contig/>" + cluster_name + "_/g\" " + algorithm_dir + "/" + seq_file)
                split_fasta(algorithm_dir + "/" + seq_file, contig_seq_dir)
                os.system("sed -i \"s/.*\.fsa\.Contig/>" + cluster_name + "_/g\" " + algorithm_dir + "/" + seq_file + ".qual")
                split_qual(algorithm_dir + "/" + seq_file + ".qual", contig_qual_dir)
    else:
        print "Unknown assembly algorithm"
 
def parse_singletons_fasta(singleton_seq_dir, singleton_qual_dir, algorithm, algorithm_dir):
    os.chdir(algorithm_dir)
    # cap3
    if (algorithm == "cap3"):
        for seq_file in glob.glob("*.fsa.cap.singlets"):
            cluster_name = seq_file.split(".")[0]
            if(os.path.getsize(algorithm_dir + "/" + seq_file) > 0): # if singletons exists
                split_fasta(algorithm_dir + "/" + seq_file, singleton_seq_dir)
    # phrap            
    elif (algorithm == "phrap"):
        for seq_file in glob.glob("*.fsa.singlets"):
            cluster_name = seq_file.split(".")[0]
            if(os.path.getsize(algorithm_dir + "/" + seq_file) > 0): # if singletons exists
                split_fasta(algorithm_dir + "/" + seq_file, singleton_seq_dir)
    else:
        print "Unknown assembly algorithm"

        
# Singletons can be reported in the ace file and not in the *.singlets 
def parse_singletons_fasta_in_ace(contig_ace_dir, singleton_seq_dir):
    # get contig info
    os.chdir(contig_ace_dir)
    for ace_file in sorted(glob.glob("*.ace")):
            ace_record = Ace.read(open(ace_file))
            contigs = ace_record.contigs
            for contig in contigs:
                if contig.nreads == 1:
                    singleton_name = contig.reads[0].rd.name
                    singleton_seq = Seq(contig.reads[0].rd.sequence)
                    singleton_record = SeqRecord(seq=singleton_seq, id="", name="", description=singleton_name)
                    singleton_file = singleton_seq_dir + "/" + singleton_name + ".fsa"
                    singleton_fd = open(singleton_file, "w")
                    SeqIO.write([singleton_record], singleton_fd, "fasta")
                    singleton_fd.close()
                    os.system("sed -i \"s/> />/g\" " + singleton_file)


def parse_singletons_qual(singleton_seq_dir, singleton_qual_dir, algorithm_dir, tmp_dir):
    os.chdir(tmp_dir)
    os.system("cat " + algorithm_dir + "/*.fsa.qual " + " > tmp.qual")
    print ("cat " + algorithm_dir + "/*.fsa.qual " + " > tmp.qual")
    split_qual_cap3("tmp.qual", tmp_dir) # easiest to split with split_qual_cap3 easiest, works for all
    os.chdir(singleton_seq_dir)
    for singleton_fsa_file in sorted(glob.glob("*.fsa")):
        singleton_qual_file = singleton_fsa_file.split(".")[0] + ".qual"
        os.system("mv " +  tmp_dir + "/" + singleton_qual_file + " " + singleton_qual_dir) # just assume it is there, otherwise something went wrong else where
    
    
def parse_ace(contig_ace_dir, algorithm, algorithm_dir):
    os.chdir(algorithm_dir)
     # cap3
    if (algorithm == "cap3"):
        for seq_file in glob.glob("*.fsa.cap.contigs"):
            cluster_name = seq_file.split(".")[0]
            if(os.path.getsize(algorithm_dir + "/" + seq_file) > 0): # if contig exists
                ace_file = seq_file.rstrip('contigs') + 'ace'
                os.system("sed -i \"s/CO Contig/CO " + cluster_name + "_/g\" " + algorithm_dir + "/" + ace_file)
                os.rename(algorithm_dir + "/" + ace_file, contig_ace_dir + "/" + cluster_name + ".ace")
    # phrap            
    elif (algorithm == "phrap"):
        for seq_file in glob.glob("*.fsa.contigs"):
            cluster_name = seq_file.split(".")[0]
            if(os.path.getsize(algorithm_dir + "/" + seq_file) > 0): # if contigs exists
                ace_file = seq_file.rstrip('contigs') + 'ace'
                os.system("sed -i \"s/CO Contig/CO " + cluster_name + "_/g\" " + algorithm_dir + "/" + ace_file)
                os.rename(algorithm_dir + "/" + ace_file, contig_ace_dir + "/" + cluster_name + ".ace")
                # HACK - sort out later
                # BioPyton does not read phrap generated ace file correctly
                # Need to replace "DS" tag string  in ACE file with a "DS " string
                os.system("sed -i \"s/DS/DS /g\" " + contig_ace_dir + "/" + cluster_name + ".ace")
                 

def generate_contig_images(contig_image_dir, contig_ace_dir):
    os.chdir(contig_ace_dir)
    for ace_file in glob.glob("*.ace"):
        contig = ace_file.split(".")[0]
        os.system("contigimage.pl " + ace_file + " " + contig_image_dir)
        
# A summary of all the contigs and their sequences        
def write_contig_summary_a(contig_ace_dir, singleton_seq_dir, summary_file):
    fd_summary = open(summary_file, 'w')
    nr_contigs = 0
    nr_singletons = 0
    summary = ""
    # get contig info
    os.chdir(contig_ace_dir)
    for ace_file in sorted(glob.glob("*.ace")):
            ace_record = Ace.read(open(ace_file))
            contigs = ace_record.contigs
            for contig in contigs:
                # do not write singletons if found in ace file
                if contig.nreads == 1:
                    continue
                summary = summary + contig.name
                for read in contig.reads:
                    summary = summary + "\t" + read.rd.name
                summary = summary + "\n"
                nr_contigs+=1
    # get singleton info
    os.chdir(singleton_seq_dir)
    for singleton_file in sorted(glob.glob("*.fsa")):
        summary = summary + singleton_file.split(".")[0]
        summary = summary + "\n"
        nr_singletons+=1
    header = "# nr Contigs: " + str(nr_contigs) + "\n"
    header = header + "# nr Singletons: " + str(nr_singletons) +  "\n"
    header = header + "# Column 1: Contig_id or Singleton_id" + "\n"
    header = header + "# Columns 2 to n: Member Sequences" + "\n"
    summary = header + summary
    fd_summary.write(summary)
    fd_summary.close()

# A summary of all cluster with more then one contig and their sequences    
def write_contig_summary_b(contig_ace_dir, singleton_seq_dir, summary_file):
    fd_summary = open(summary_file, 'w')
    nr_contigs = 0
    summary = ""
    # get contig info
    os.chdir(contig_ace_dir)
    for ace_file in sorted(glob.glob("*.ace")):
            ace_record = Ace.read(open(ace_file))
            if (ace_record.ncontigs > 1):
                contigs = ace_record.contigs
                for contig in contigs:
                    # do not write singletons if found in ace file
                    if contig.nreads == 1:
                        continue
                    summary = summary + contig.name
                    for read in contig.reads:
                        summary = summary+ "\t" + read.rd.name
                    summary = summary + "\n"
                    nr_contigs+=1
    header = "# nr Contigs: " + str(nr_contigs) + "\n"
    header = header + "# Column 1: Contig_id" + "\n"
    header = header + "# Columns 2 to n: Member Sequences" + "\n"
    summary = header + summary
    fd_summary.write(summary)
    fd_summary.close()
    
def split_fasta(seq_file, split_seq_dir):
    fd_seq = open(seq_file, "rU")
    seqs = SeqIO.parse(fd_seq, "fasta")
    for seq in seqs:
        new_seq_file = split_seq_dir + "/" + seq.name + ".fsa"
        fd_split_seq = open(new_seq_file, "w")
        SeqIO.write([seq], fd_split_seq, "fasta")
        fd_split_seq.close()
        os.system("sed -i \"s/> />/g\" " + new_seq_file)
    fd_seq.close()

def split_qual(qual_file, split_qual_dir):
    fd_qual = open(qual_file, "rU")
    quals = SeqIO.parse(fd_qual, "qual")
    for qual in quals:
        new_qual_file = split_qual_dir + "/" + qual.name + ".qual"
        fd_split_qual = open(new_qual_file, "w")     
        SeqIO.write([qual], fd_split_qual, "qual")
        fd_split_qual.close()
        os.system("sed -i \"s/> />/g\" " + new_qual_file)
    fd_qual.close()
    
# Cap3 qual scores are not in the BioPython specified range
# Need to split differently
def split_qual_cap3(qual_file, split_qual_dir):
    fd_qual = open(qual_file, "rU")
    start = True
    header = ""
    seq_name = ""
    header_pattern = re.compile('^>.*')
    for line in fd_qual:
        if re.search(header_pattern, line):
            if not start:
                new_qual_file = split_qual_dir + "/" + seq_name + ".qual"
                fd_split_qual = open(new_qual_file, "w")
                fd_split_qual.write(header)
                fd_split_qual.write(scores)
                fd_split_qual.close()  
            header = line
            seq_name = line.lstrip('>').rstrip('\n').split(' ')[0]
            scores = ""
            start = False
        else:
            scores = scores + line
    # Write last file
    new_qual_file = split_qual_dir + "/" + seq_name + ".qual"
    fd_split_qual = open(new_qual_file, "w")
    fd_split_qual.write(header)
    fd_split_qual.write(scores)
    fd_split_qual.close()  
    
    fd_qual.close()

if __name__ == "__main__":
#    split_qual_cap3('XHC00016.fsa.cap.contigs.qual', 'tmp')
    sys.exit(main())


