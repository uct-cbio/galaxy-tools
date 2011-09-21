#!/usr/bin/python

# Wrapper for the console version of InterProScan

# /usr/local/share/paeano/software/interproscan/release/iprscan/bin/iprscan -cli -i Galaxy39-peptide.fasta -o Galaxy39-peptide.fasta.xml -seqtype p -goterms -iprlookup -format xml -altjobs > log.txt 2>&1

import sys, re, string, commands, os, glob
from optparse import OptionParser
from Utils import SystemUtils

def main():
    
    usage = "usage: -i SEQ_FILE -t SEQ_TYPE -j JOB_ID_FILE -o INTERPROSCAN_XML_FILE -w PROG_WORK_DIR -d"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--seq", dest="seq_file", help="A FASTA sequence file.")
    parser.add_option("-t", "--seq_type", dest="seq_type", help="Type of sequence. n for DNA/RNA, p for protein (default). ")
    parser.add_option("-j", "--job_id", dest="job_id_file", help="The InterProScan job id of the run. Stored in a file.")
    parser.add_option("-o", "--interproscan_xml", dest="interproscan_xml_file", help="InterProScan XML results file.")
    parser.add_option("-w", "--prog_work", dest="prog_work_dir", help="Program working directory, contains all processed files.")
    parser.add_option("-d", "--delete_program_work", action="store_true", dest="delete_program_work", default=False, help="Delete program working directory after program has completed.")
    (options, args) = parser.parse_args()
    
    if not options.seq_file:
        print "Please specify the FASTA sequence file (-i SEQ_FILE)"
        return - 1
#    if not options.seq_type:
#        print "Please specify the type of sequence. n for DNA/RNA, p for protein (default) (-t SEQ_TYPE)"
#        return - 2
    if not options.job_id_file:
        print "Please specify the InterProScan job id file (-o JOB_ID_FILE)"
        return - 3
    if not options.interproscan_xml_file:
        print "Please specify the InterProScan XML results file (-o INTERPROSCAN_XML_FILE)"
        return - 3
    if not options.prog_work_dir:
        print "Please specify the program working directory (-w PROG_WORK_DIR)"
        return - 4
    if (len(args) > 0):
        print "Too many input arguments"
        return - 5
    
    # Do some initialization
    print "Initialize..." 
    root_dir = os.getcwd()
    
    # Get full file paths
    seq_file = os.path.abspath(options.seq_file)
    if options.seq_type:
        seq_type = options.seq_type
    else:
        seq_type = "p" # default
    job_id_file = os.path.abspath(options.job_id_file)    
    interproscan_xml_file = os.path.abspath(options.interproscan_xml_file)    
    prog_work_dir = os.path.abspath(options.prog_work_dir)
    
    # Create working and base directory 
    timestamp = commands.getoutput("date +%Y-%m-%d_%H_%M_%S_%N")
    base_dir = prog_work_dir + "/interproscan_" + timestamp 

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

    # Run InterProScan
    print "Run InterProScan..."
    os.chdir(base_dir)
    interproscan_exec = SystemUtils.find_program("iprscan")
    os.system(interproscan_exec + " -cli  -i " +  seq_file + " -seqtype " + seq_type + " -goterms -iprlookup -format xml -altjobs -o " + interproscan_xml_file + " 2> job.id")
    print (interproscan_exec + " -cli -i " +  seq_file + " -seqtype " + seq_type + " -goterms -iprlookup -format xml -altjobs -o " + interproscan_xml_file + " 2> job.id")
    os.system("sed \"s/SUBMITTED\s//\" job.id > " + job_id_file)
    
    
    # Delete program working directory if indicated
    if(options.delete_program_work):
        print "Delete working directory"
        os.system("rm -rf " + prog_work_dir)
        
    # Done
    print "Done."
    return 0

if __name__ == "__main__":
    sys.exit(main())


