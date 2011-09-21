#!/usr/bin/python

# Split a combined xml InterProScan file into individual files. This files can then be imported into the GUI version of Blast2GO.

import sys, re, string, commands, os, glob
from optparse import OptionParser
from Utils import SystemUtils

def main():
    
    usage = "usage: -i INTERPROSCAN_XML_COMBINED_FILE -o INTERPROSCAN_XML_SINGLEL_FILE -w PROG_WORK_DIR -d"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--interproscan_combined_xml", dest="interproscan_combined_xml", help="InterProScan XML file with combined results.")
    parser.add_option("-o", "--interproscan_single_xml", dest="interproscan_single_xml", help="Zip archive containing individual InterProScan XML results.")
    parser.add_option("-w", "--prog_work", dest="prog_work_dir", help="Program working directory, contains all processed files.")
    parser.add_option("-d", "--delete_program_work", action="store_true", dest="delete_program_work", default=False, help="Delete program working directory after program has completed.")
    (options, args) = parser.parse_args()
    
    if not options.interproscan_combined_xml:
        print "Please specify the InterProScan XML file with combined results (-i INTERPROSCAN_XML_COMBINED_FILE)"
        return - 1
    if not options.interproscan_single_xml:
        print "Please specify the zip archive containing individual InterProScan XML results (-o INTERPROSCAN_XML_SINGLE_FILE)"
        return - 2
    if not options.prog_work_dir:
        print "Please specify the program working directory (-w PROG_WORK_DIR)"
        return - 3
    if (len(args) > 0):
        print "Too many input arguments"
        return - 4
    
    # Do some initialization
    print "Initialize..." 
    root_dir = os.getcwd()
    
    # Get full file paths
    interproscan_combined_xml = os.path.abspath(options.interproscan_combined_xml)
    interproscan_single_xml = os.path.abspath(options.interproscan_single_xml)    
    prog_work_dir = os.path.abspath(options.prog_work_dir)
    
    # Create working and base directory 
    timestamp = commands.getoutput("date +%Y-%m-%d_%H_%M_%S_%N")
    base_dir = prog_work_dir + "/split_interproscan_results" + timestamp 

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

    # Run split InterProScan results
    print "Run split InterProScan results..."
    os.chdir(base_dir)
    interproscan_dir = base_dir + "/interproscan"
    split_iterproscan_blast2go_exec = SystemUtils.find_program("split_interproscan_blast2go.pl")
    os.system(split_iterproscan_blast2go_exec + " -interproscan_xml " + interproscan_combined_xml + " -blast2go_xml " +  interproscan_dir)
    
    # Prepare output
    print "Prepare output..."
    tmp_zip = base_dir + "/tmp.zip" # Galaxy work around need to create a temporary zip archive and move to the output data set
    
#    os.system("zip -j " + interproscan_single_xml + " " + interproscan_dir + "/*")
    os.system("zip -j " + tmp_zip + " " + interproscan_dir + "/*")
    os.system("mv " + tmp_zip + " " + interproscan_single_xml)
    
    # Delete program working directory if indicated
    if(options.delete_program_work):
        print "Delete working directory"
        os.system("rm -rf " + prog_work_dir)
        
    # Done
    print "Done."
    return 0

if __name__ == "__main__":
    sys.exit(main())


