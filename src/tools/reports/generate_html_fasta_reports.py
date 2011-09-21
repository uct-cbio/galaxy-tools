#!/usr/bin/python

# Program takes as input a file with FASTA sequences. An HTML output page with the links to the FASTA sequences are generated.

# The output HTML page will contain links to the report file names. Not to the absolute path! This is specifically done so that the output page and links 
# will work within galaxy. 

# The output HTML page is generate with the Python HTML template library, HTMLTemplate-1.5.0. Download it from here: http://pypi.python.org/pypi/HTMLTemplate/

import sys, re, string, commands, os, glob
from optparse import OptionParser
from Utils import SystemUtils
import imghdr
from HTMLTemplate import Template
from Bio import SeqIO

def main():
    usage = "usage: %prog -f FASTA -p OUT_PAGE -o OUT_DIR -w PROG_WORK_DIR -d"
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--fasta", dest="fasta_file", help="File with FASTA sequences.")
    parser.add_option("-p", "--out_page", dest="out_page", help="HTML output page with links to FASTA sequences.")
    parser.add_option("-o", "--out", dest="out_dir", help="Output directory containing the FASTA sequences.")
    (options, args) = parser.parse_args()
    
    if not options.fasta_file:
        print "Please specify the file with FASTA sequences (-f FASTA)"
        return -1
    if not options.out_page:
        print "Please specify HTML output page (-p OUT_PAGE)"
        return -2
    if not options.out_dir:
        print "Please specify the output directory (-o OUT_DIR)"
        return -3
    if (len(args) > 0):
        print "Too many input arguments"
        return -4
    
    fasta_file = os.path.abspath(options.fasta_file)
    out_page = os.path.abspath(options.out_page)
    out_dir = os.path.abspath(options.out_dir)

    if not os.path.isdir(out_dir):
        os.system("mkdir " + out_dir)
    else:
        os.system("rm -rf " + out_dir)
        os.system("mkdir " + out_dir) 
        
    split_fasta(fasta_file, out_dir)
    os.chdir(out_dir)
    files = os.listdir('.')
    html_fastas = []
    for file in files:
        print file
        html_fasta = HTMLFasta(file,file)
        html_fastas.append(html_fasta)
            
    # Generate output page
    html_page_with_fastas = HTMLPageWithFastas("HTML FASTA report", html_fastas)
    html_page_with_fastas_template_str = get_html_with_fastas_template()
    html_page_with_fastas_template = Template(render_html_with_fastas_template, html_page_with_fastas_template_str)
    html_page_with_fastas_html = html_page_with_fastas_template.render(html_page_with_fastas)       
    fd_out = open(out_page,"w")
    fd_out.write(html_page_with_fastas_html)
    fd_out.close()
        
    # Done
    print "Done."
    return 0
###########################################################################
# Split FASTA file into individual sequences
def split_fasta(seq_file, split_seq_dir):
    fd_seq = open(seq_file, "rU")
    seqs = SeqIO.parse(fd_seq, "fasta")
    for seq in seqs:
#        new_seq_file = split_seq_dir + "/" + seq.name + ".fsa"
        new_seq_file = split_seq_dir + "/" + seq.name + ".txt"
        fd_split_seq = open(new_seq_file, "w")
        SeqIO.write([seq], fd_split_seq, "fasta")
        fd_split_seq.close()
        os.system("sed -i \"s/> />/g\" " + new_seq_file)
    fd_seq.close()


###########################################################################
# Classes and methods to create the HTML file with links to FASTA files
class HTMLPageWithFastas:
    def __init__(self, title, html_fastas):
        self.title = title
        self.html_fastas = html_fastas
        
class HTMLFasta:
    def __init__(self, title, path):
        self.title = title
        self.path = path
        
def get_html_with_fastas_template():
    html_with_fastas_template = '''<html>
    <head>
        <title node="con:hpwf_title">HTML page with FASTA files</title>
    </head>
    <body>

        <h1 node="con:hpwf_header">HTML page with FASTA file header</h1>
           
        <table cellpadding="3" border="1" cellspacing="1" width="100%">
            <thead>
                <tr>
                    <th>Fasta Sequences</th>
                </tr>
            </thead>
            <tbody>
                <tr node="rep:hpwf_fasta">
                    <td><a node="con:html_fasta" href="http://seq.fsa">FASTA Sequences</a></td>
                </tr>
            </tbody>
        </table>
    
    </body>
</html>'''
    return html_with_fastas_template

def render_html_with_fastas_template(node, hpwf_fastas):
    node.hpwf_title.content = hpwf_fastas.title
    node.hpwf_header.content = hpwf_fastas.title
    node.hpwf_fasta.repeat(render_html_fasta, hpwf_fastas.html_fastas)
    
def render_html_fasta(node, html_fasta):
    node.html_fasta.content = html_fasta.title
    node.html_fasta.atts['href'] = html_fasta.path
###########################################################################

if __name__ == "__main__":
    sys.exit(main())


