#!/usr/bin/python

# Run console version of Blast2G0 (b2g4pipe)

# nohup time java -jar ~/blast2go/b2g4pipe/latest/pipeline/blast2go.jar -in blastp_nr.xml -out all_blastp_nr -prop blast2go.properties -annot -dat -v -ips interproscan/ -img -wiki ~/blast2go/b2g4pipe/latest/pipeline/wiki_arrayannot.txt > log.txt 2>&1 &

import sys, re, string, commands, os, glob
from optparse import OptionParser
from Utils import SystemUtils

def main():
    
    usage = "usage: -b BLAST_REPORT_FILE -i INTERPRO_XML_FILE -a BLAST2GO_ANNOTATION_FILE -g BLAST2GO_DAT_FILE -m BLAST2GO_IMAGE_FILE -w PROG_WORK_DIR -d"
    parser = OptionParser(usage=usage)
    parser.add_option("-b", "--blast_report", dest="blast_report_file", help="A BLAST report file containing concatenated BLAST XML results.")
    parser.add_option("-i", "--interproscan_xml", dest="interproscan_xml_file", help="An InterProScan XML results file. Optional. \
    If included results will be used in annotation.")
    parser.add_option("-g", "--blast2go_dat", dest="blast2go_dat_file", help="Blast2GO output dat file. This file can be imported into the Blast2GO GUI version.")
    parser.add_option("-a", "--blastgo_annotation", dest="blastgo_annotation_file", help="Blast2GO annotation output file.")
    parser.add_option("-m", "--blast2go_image", dest="blastgo_image_file", help="Blast2GO zipped output file containing statistic charts.")
    parser.add_option("-w", "--prog_work", dest="prog_work_dir", help="Program working directory, contains all processed files.")
    parser.add_option("-d", "--delete_program_work", action="store_true", dest="delete_program_work", default=False, help="Delete program working directory after program has completed.")
    (options, args) = parser.parse_args()
    
    if not options.blast_report_file:
        print "Please specify the file containing concatenated BLAST XML results (-b BLAST_REPORT_FILE)"
        return - 1
#    if not options.interproscan_xml_file:
#        print "Please specify the InterProScan XML results file (-i INTERPRO_XML_FILE)"
#        return - 1
    if not options.blast2go_dat_file:
        print "Please specify the Blast2GO dat output file (-g BLAST2GO_DAT_FILE)"
        return - 3
    if not options.blastgo_annotation_file:
        print "Please specify the Blast2Go annotation output file (-a BLAST2GO_ANNOTAION_FILE)"
        return - 4
    if not options.blastgo_image_file:
        print "Please specify the Blast2GO zipped output file containing statistic charts  (-m BLAST2GO_IMAGE_FILE)"
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
    blast_report_file = os.path.abspath(options.blast_report_file)
    if options.interproscan_xml_file:
        interproscan_xml_file = os.path.abspath(options.interproscan_xml_file)
    blast2go_dat_file = os.path.abspath(options.blast2go_dat_file)
    blastgo_annotation_file = os.path.abspath(options.blastgo_annotation_file)
    blastgo_image_file = os.path.abspath(options.blastgo_image_file)
    prog_work_dir = os.path.abspath(options.prog_work_dir)
    
    # Create working and base directory 
    timestamp = commands.getoutput("date +%Y-%m-%d_%H_%M_%S_%N")
    base_dir = prog_work_dir + "/b2g4pipe_" + timestamp
    b2g4pipe_log_file = base_dir + "/b2g4pipe_log.txt"

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
    blast2go_jar = os.environ.get("BLAST2GO_JAR_FILE")
    blast2go_properties = os.environ.get("BLAST2GO_PROPERTIES_FILE")

    # Run Blast2GO
    print "Run blast2go..."
    os.chdir(base_dir)
    # Check if InterProScan results are provided
    if options.interproscan_xml_file:
        interproscan_dir = base_dir + "/interproscan"
        split_iterproscan_blast2go_exec = SystemUtils.find_program("split_interproscan_blast2go.pl")
        os.system(split_iterproscan_blast2go_exec + " -interproscan_xml " + interproscan_xml_file + " -blast2go_xml " +  interproscan_dir)
        print (split_iterproscan_blast2go_exec + " -interproscan_xml " + interproscan_xml_file + " -blast2go_xml " +  interproscan_dir)
        os.system("java -jar " +  blast2go_jar + " -prop " +   blast2go_properties + " -in " + blast_report_file + " -ips " + interproscan_dir + "/ -annot -dat -v -img -out blast2go > " + b2g4pipe_log_file) # need to add a "/" after the interproscan_dir
    else:
        os.system("java -jar " +  blast2go_jar + " -prop " +   blast2go_properties + " -in " + blast_report_file + " -annot -dat -v -img -out blast2go > " + b2g4pipe_log_file)

    # Prepare output
    print "Prepare output..."
    os.system("mv " + base_dir + "/blast2go.dat" + " " + blast2go_dat_file)    
    os.system("mv " + base_dir + "/blast2go.annot" + " " + blastgo_annotation_file)    
    tmp_zip = base_dir + "/tmp.zip" # Galaxy work around need to create a temporary zip archive and move to the output data set
    os.system("zip -qqj " + tmp_zip + " " + base_dir + "/*.png")
    os.system("mv " + tmp_zip + " " + blastgo_image_file)
    
    # Delete program working directory if indicated
    if(options.delete_program_work):
        print "Delete working directory"
        os.system("rm -rf " + prog_work_dir)
        
    # Done
    print "Done."
        
    return 0

if __name__ == "__main__":
    sys.exit(main())

