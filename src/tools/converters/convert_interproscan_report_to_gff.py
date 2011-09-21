#!/usr/bin/python

# The script calls the converter.pl script included with the iprscan installation and generates a gff3 file for the specified job id

# Can later include the conversion to txt and xml as well 

import sys
import re
import string
from optparse import OptionParser
import os
import glob
import commands
from Utils import SystemUtils

def main():
    usage = "usage: %prog -j JOB_ID -o OUT_FILE"
    parser = OptionParser(usage=usage)
    parser.add_option("-j", "--job_id", dest="job_id", help="JobId of the InterProScan run.")
    parser.add_option("-o", "--out", dest="out_file", help="GFF3 output file.")
    (options, args) = parser.parse_args()
    
    
    if not options.job_id:
        print "Please specify the InterProScan JobId (-j JOB_ID)"
        return -1
    if not options.out_file:
        print "Please specify the GFF3 output file (-o OUT)"
        return -2
    if (len(args) > 0):
        print "Too many input arguments"
        return -3
    
    out_file = os.path.abspath(options.out_file)
    iprscan_output_dir = os.environ.get("IPRSCAN_OUTPUT_DIR")
    date = options.job_id.split("-")[1]
    merged_raw_file = iprscan_output_dir + "/" + date + "/" + options.job_id + "/merged.raw"
    converter_exec = SystemUtils.find_program("converter.pl")
    os.system(converter_exec + " -format gff3  -input " + merged_raw_file  + " -output " + out_file)
    return 0

if __name__ == "__main__":
    sys.exit(main())


