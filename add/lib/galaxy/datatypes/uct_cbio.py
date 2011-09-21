"""
UCT CBIO classes
"""
import logging, os, sys, time, sets, tempfile, shutil, subprocess, data, binascii, zipfile, gzip, urllib
from cgi import escape
from galaxy import util
from galaxy.datatypes.sniff import *
from bx.intervals.io import *
from galaxy.datatypes import metadata
from galaxy.datatypes.metadata import MetadataElement
from galaxy.datatypes.binary import Binary

# If a binary extension needs support add it to unsniffable_binary_formats in binary.py 
# This is necessary to sucessfully load the binary using upload.py tool
# e.q ['fa', 'fsa', 'fasta'] for a zip file containing fasta sequences

# Generic binary class
class GenBin( Binary ):
    """Binary file"""
    file_ext = "dat"

    def __init__(self, **kwd):
        """Initialize binary datatype"""
        Binary.__init__(self, **kwd)

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            dataset.peek  = "Binary file"
            dataset.blurb = data.nice_size( dataset.get_size() )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'
    def display_peek( self, dataset ):
        try:
            return dataset.peek
        except:
            return "Binary file (%s)" % ( data.nice_size( dataset.get_size() ) )

# Generic zip class
class GenZip( Binary ):
    """Binary zip file containing any data"""
    file_ext = "zip"

    def __init__(self, **kwd):
        """Initialize zip datatype"""
        Binary.__init__(self, **kwd)

    def set_peek( self, dataset, is_multi_byte=False ):
        if not dataset.dataset.purged:
            zip_file = zipfile.ZipFile( dataset.file_name, "r" )
            num_files = len( zip_file.namelist() )
            dataset.peek  = "Archive of %s files" % ( str( num_files ) )
            dataset.blurb = data.nice_size( dataset.get_size() )
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'
    def display_peek( self, dataset ):
        try:
            return dataset.peek
        except:
            return "Binary file archive (%s)" % ( data.nice_size( dataset.get_size() ) )
    def get_mime( self ):
        """Returns the mime type of the datatype"""
        return 'application/zip'

