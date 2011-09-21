#!/usr/bin/python

# Program takes as input a zip archive of images. An HTML output page with the images are generated.

# The output HTML page will contain links to the report file names. Not to the absolute path! This is specifically done so that the output page and links 
# will work within galaxy. 

# The output HTML page is generate with the Python HTML template library, HTMLTemplate-1.5.0. Download it from here: http://pypi.python.org/pypi/HTMLTemplate/

import sys, re, string, commands, os, glob
from optparse import OptionParser
from Utils import SystemUtils
import imghdr
from HTMLTemplate import Template

def main():
    usage = "usage: %prog -z ZIP_ARCHIVE -p OUT_PAGE -o OUT_DIR -w PROG_WORK_DIR -d"
    parser = OptionParser(usage=usage)
    parser.add_option("-z", "--zip_archive", dest="zip_archive", help="Zip archive consisting of images.")
    parser.add_option("-p", "--out_page", dest="out_page", help="HTML output page with images.")
    parser.add_option("-o", "--out", dest="out_dir", help="Output directory with images.")
    parser.add_option("-w", "--prog_work", dest="prog_work_dir", help="Program working directory, contains all processed files.")
    parser.add_option("-d", "--delete_program_work", action="store_true", dest="delete_program_work", default=False, help="Delete program working directory after program has completed.")
    (options, args) = parser.parse_args()
    
    if not options.zip_archive:
        print "Please specify the zip archive consisting of images (-z ZIP_ARCHIVE)"
        return -1
    if not options.out_page:
        print "Please specify HTML output page (-p OUT_PAGE)"
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
    
    zip_archive = os.path.abspath(options.zip_archive)
    out_page = os.path.abspath(options.out_page)
    out_dir = os.path.abspath(options.out_dir)

    prog_work_dir = os.path.abspath(options.prog_work_dir)
    timestamp = commands.getoutput("date +%Y-%m-%d_%H_%M_%S_%N")
    base_dir = prog_work_dir + "/generate_html_with_images_report_" + timestamp
    
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
        
    os.chdir(base_dir)
    os.system("unzip " + zip_archive)
    
    # Simple test to check if the files are images. Add only the images to the HTML page to be rendered.
    files = os.listdir('.')
    html_images = []
    for file in files:
        if(imghdr.what(file)):
            print file
            html_image = HTMLImage(file,file)
            html_images.append(html_image)
            os.system("cp " +  file + " " + out_dir )
        else:
            print file + " is not an image file!"
            
    # Generate output page
    html_page_with_images = HTMLPageWithImages("HTML with images report",html_images)
    html_page_with_images_template_str = get_html_with_images_template()
    html_page_with_images_template = Template(render_html_with_images_template, html_page_with_images_template_str)
    html_page_with_images_html = html_page_with_images_template.render(html_page_with_images)       
    fd_out = open(out_page,"w")
    fd_out.write(html_page_with_images_html)
    fd_out.close()
        
    # Delete program working directory if indicated
    if(options.delete_program_work):
        print "Delete working directory"
        os.system("rm -rf " + base_dir)
        
    # Done
    print "Done."
    return 0

###########################################################################
# Classes and methods to create the HTML file with images index page
class HTMLPageWithImages:
    def __init__(self, title, images):
        self.title = title
        self.images = images
        
class HTMLImage:
    def __init__(self, title, path):
        self.title = title
        self.path = path
        
def get_html_with_images_template():
    html_with_images_template = '''<html>
    <head>
        <title node="con:hpwi_title">HTML page with images title</title>
    </head>
    <body>

        <h1 node="con:hpwi_header">HTML page with images header</h1>
           
        <table cellpadding="3" border="1" cellspacing="1" width="100%">
            <thead>
                <tr>
                    <th>Images</th>
                </tr>
            </thead>
            <tbody>
                <tr node="rep:hpwi_image">
                    <td><img node="con:image" src="image.png" alt="image.png" /></td>
                </tr>
            </tbody>
        </table>
    
    </body>
</html>'''
    return html_with_images_template

def render_html_with_images_template(node, hpwi_images):
    node.hpwi_title.content = hpwi_images.title
    node.hpwi_header.content = hpwi_images.title
    node.hpwi_image.repeat(render_html_image, hpwi_images.images)
    
def render_html_image(node, html_image):
    node.image.content = html_image.title
    node.image.atts['src'] = html_image.title
    node.image.atts['alt'] = html_image.path
###########################################################################

if __name__ == "__main__":
    sys.exit(main())


