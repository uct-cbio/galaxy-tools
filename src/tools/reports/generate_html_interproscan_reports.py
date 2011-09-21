#!/usr/bin/python

# The script goes through 3 main steps
# 1) The iprscan tool generates html report output pages pages for each chunk (the job id is needed to do this!)
# 2) The html pages are cleaned up
# 3) A job index page is created which links to the report pages

# The job index page will contain links to the report file names. Not to the absolute path! This is specifically done so that the output page and links 
# will work within galaxy. 

# The job html page is generate with the Python html template library, HTMLTemplate-1.5.0. Download it from here: http://pypi.python.org/pypi/HTMLTemplate/

import sys
import re
import string
from optparse import OptionParser
import os
import glob
import commands
from HTMLTemplate import Template

def main():
    usage = "usage: %prog -j JOB_ID -p JOB_HTML_PAGE -o OUT_DIR -w PROG_WORK_DIR -d"
    parser = OptionParser(usage=usage)
    parser.add_option("-j", "--job_id", dest="job_id", help="JobId of the InterProScan run.")
    parser.add_option("-p", "--job_html_page", dest="job_html_page", help="InterProScan job html output page.")
    parser.add_option("-o", "--out", dest="out_dir", help="Output directory with single html IterProScan result.")
    parser.add_option("-w", "--prog_work", dest="prog_work_dir", help="Program working directory, contains all processed files.")
    parser.add_option("-d", "--delete_program_work", action="store_true", dest="delete_program_work", default=False, help="Delete program working directory after program has completed.")
    (options, args) = parser.parse_args()
    
    if not options.job_id:
        print "Please specify the InterProScan JobId (-j JOB_ID)"
        return -1
    if not options.job_html_page:
        print "Please specify InterProScan job html output page (-p JOB_HTML_PAGE)"
        return -2
    if not options.out_dir:
        print "Please specify the output directory (-o OUT_DIR)"
        return -4
    if not options.prog_work_dir:
        print "Please specify the program working directory (-w PROG_WORK_DIR)"
        return -4
    if (len(args) > 0):
        print "Too many input arguments"
        return -5
    
    iprscan_output_dir = os.environ.get("IPRSCAN_OUTPUT_DIR")
    iprscan_images_dir = os.environ.get("IPRSCAN_IMAGES_DIR")
    out_dir = os.path.abspath(options.out_dir)
    job_html_page = os.path.abspath(options.job_html_page)
    date = options.job_id.split("-")[1]
    iprscan_output_dir = iprscan_output_dir + "/" + date + "/" + options.job_id
    nr_chunks = len(glob.glob(iprscan_output_dir + "/chunk*"))
    
    prog_work_dir = os.path.abspath(options.prog_work_dir)
    timestamp = commands.getoutput("date +%Y-%m-%d_%H_%M_%S_%N")
    base_dir = prog_work_dir + "/generate_interproscan_html_report_" + timestamp
    
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
       return -5
   
    if not os.path.isdir(out_dir):
        os.system("mkdir " + out_dir)
    else:
        os.system("rm -rf " + out_dir)
        os.system("mkdir " + out_dir) 
    
    iprscan_raw_html_file = os.path.abspath(base_dir + "/iprscan_raw.html")
    
    for i in range(1,int(nr_chunks) + 1):
        # create picture output
        os.system("iprscan tool=iprscan jobid=" + options.job_id + " cnk=" + str(i) + " view=picture  >> " + iprscan_raw_html_file)
        # create table output
        os.system("iprscan tool=iprscan jobid=" + options.job_id + " cnk=" + str(i) + " view=Table  >> " + iprscan_raw_html_file)    


    fd_in = open(iprscan_raw_html_file, "r")
    
    ips_report_start = False
    ips_report = ""
    
    ips_report_start_string = "<table width=\"100%\" cellSpacing=\"0\" cellPadding=\"3\" border=\"1\" bordercolor=\"#999999\" bgcolor=\"#eeeeee\" class=\"tabletop\">"
    ips_report_end_string = "</tr>\n</table>\n</td></tr>\n</table>\n<br>"
    
    ips_report_header = "<html><head><title>InterProScan Results - interhelp@ebi.ac.uk</title></head><link rel=\"SHORTCUT ICON\" href=\"https://foobar.com/images/bookmark.gif\">\n<br>\n<span class=\"text\">"
    ips_report_footer = "</span>"
    
    ips_reports = {}
    ips_report_id = "";

    for line in fd_in:
        if (line.find(ips_report_start_string) > -1):
            if(ips_report_start): # Do not save results for first round. ips_report_start = False initially
                if (ips_report_id in ips_reports):
                    report = ips_reports[ips_report_id]
                    report = report + ips_report
                    ips_reports[ips_report_id] = report 
                else:    
                    ips_reports[ips_report_id] = ips_report
                
            ips_report_start = True
            ips_report = ""    
            
        if(ips_report_start):
            ips_report = ips_report + line
            if (line.find("SEQUENCE") > -1):
                ips_report_id = get_ips_id(line)
            
    # Get last report
    if (ips_report_id in ips_reports):
        report = ips_reports[ips_report_id]
        report = report + ips_report
        ips_reports[ips_report_id] = report 
    else:    
        ips_reports[ips_report_id] = ips_report
        

    # Add html header and footer
    for ips_report_id in ips_reports:
        report = ips_reports[ips_report_id]
        report = ips_report_header + report + ips_report_footer
        ips_reports[ips_report_id] = report
        
    # Clean html
    for ips_report_id in ips_reports:
        report = ips_reports[ips_report_id]
        report = clean_html(report, ips_report_id)
        ips_reports[ips_report_id] = report
    
    
    interproscan_reports = []
    report_ids = ips_reports.keys()
    sorted_report_ids = sorted(report_ids)
    
    # Print all reports
    for ips_report_id in sorted_report_ids:
        fd_out = open(out_dir + "/" + ips_report_id + "_interproscan.html", "w")
        fd_out.write(ips_reports[ips_report_id])
        fd_out.close()
        interproscan_report = InterProScanReport(ips_report_id, ips_report_id + "_interproscan.html")
        interproscan_reports.append(interproscan_report)

    fd_in.close()
    
# Cannot access html pages or links in galaxy output data subfolders. Need to put all file and images in root directory. Messy! find better way.
# Need to do for now!
#    os.system("mkdir " + out_dir + "/images")
#    os.system("cp " +  iprscan_images_dir + "/*.gif " + out_dir + "/images")
    os.system("cp " +  iprscan_images_dir + "/*.gif " + out_dir )
    
    # Generate InterProScan job index page
    interproscan_job = InterProScanJob(options.job_id, interproscan_reports)
    intproscan_job_html_template_str =  get_intproscan_job_html_template()
    intproscan_job_html_template = Template(render_interproscan_job_template, intproscan_job_html_template_str)
    intproscan_job_html = intproscan_job_html_template.render(interproscan_job)       
    fd_out = open(job_html_page,"w")
    fd_out.write(intproscan_job_html)
    fd_out.close()
    
    # Delete program working directory if indicated
    if(options.delete_program_work):
        print "Delete working directory"
        os.system("rm -rf " + base_dir)
    
    return 0

###########################################################################
# Methods to cleanup iprscan html reports
def get_ips_id(line):
    ips_id = ""
    p = re.compile('.*target=\"newwindow\">(.*)<\/a>')
    m = re.search(p,line)
    if(m):
        ips_id = m.group(1)
    else:
        print "Cannot get sequence id from InterProScan report"
        sys.exit()
    return ips_id

def clean_html(report, id):
    report = report.replace("InterProScan Results - interhelp@ebi.ac.uk",id + " InterProScan Results") # add title
    report = report.replace("https://foobar.com/images","./") # this replacement is specifically for galaxy. Need all files in root folder.
    report = report.replace("https://foobar.com/","./") # all foobar links
    
    report = report.replace("http://www.ebi.ac.uk/InterProScan/images","./") # this replacement is specifically for galaxy. Need all files in root folder.
    report = report.replace("http://www.ebi.ac.uk/InterProScan/","./") # all www.ebi.ac.uk/InterProScan links
    
    # Clean up table headers    
    p = re.compile('<td colspan=\"2\"><b>SEQUENCE:</b>.*</b>') # 1st table header if view=picture in iprscan
    m = re.search(p,report)
    if(m):
        line = m.group(0)
        report = report.replace(line,"<td colspan=\"2\"><b>SEQUENCE:</b>&nbsp;" + id)
    p = re.compile('<td colspan=\"4\"><b>SEQUENCE:</b>.*</b>') # 2nd table header if view=Table in iprscan 
    m = re.search(p,report)
    if(m):
        line = m.group(0)
        report = report.replace(line,"<td colspan=\"4\"><b>SEQUENCE:</b>&nbsp;" + id)
    # Disable links found in some reports
    p = re.compile('&nbsp;&nbsp;<a href=\".*</a></td>') # 
    m = re.search(p,report)
    if(m):
        line = m.group(0)
        report = report.replace(line,"")
    
    # Disable form submit ant retrieve buttons
    p = re.compile(".*<input type=\"submit.*\">")
    m = re.findall(p,report)
    if (m):
        for match in m:
            report = report.replace(match,"")
    
    return report
###########################################################################

###########################################################################
# Classes and methods to create the InterProScan job index page
class InterProScanJob:
    def __init__(self, id, reports):
        self.id = id
        self.reports = reports
        
class InterProScanReport:
    def __init__(self, id, report):
        self.id = id
        self.report = report
        
def get_intproscan_job_html_template():
    
    interproscan_job_html_template = '''<html>
    <head>
        <title node="con:ips_job_id">JOB ID</title>
    </head>
    <body>

        <h1 node="con:ips_job_header">JOB ID</h1>
           
        <table cellpadding="3" border="1" cellspacing="1" width="100%">
            <thead>
                <tr>
                    <th>InterProScan Report</th>
                </tr>
            </thead>
            <tbody>
                <tr node="rep:ips_report">
                    <td><a node="con:report" href="http://interproscan_report.html">InterProScan Report</a></td>
                </tr>
            </tbody>
        </table>
    
    </body>
</html>'''

    return interproscan_job_html_template

def render_interproscan_job_template(node, ips_job):
    node.ips_job_id.content = ips_job.id
    node.ips_job_header.content = "InterProScan reports for job id: " + ips_job.id 
    node.ips_report.repeat(render_interproscan_report, ips_job.reports)
    
def render_interproscan_report(node, ips_report):
    node.report.content = ips_report.id
    node.report.atts['href'] = ips_report.report
###########################################################################

if __name__ == "__main__":
    sys.exit(main())


