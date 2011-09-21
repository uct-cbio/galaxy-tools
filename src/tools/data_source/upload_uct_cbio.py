#!/usr/bin/env python
#Processes uploads from the user.

# WARNING: Changes in this tool (particularly as related to parsing) may need
# to be reflected in galaxy.web.controllers.tool_runner and galaxy.tools

import urllib, sys, os, gzip, tempfile, shutil, re, gzip, zipfile, codecs, binascii
from galaxy import eggs
# need to import model before sniff to resolve a circular import dependency
import galaxy.model
from galaxy.datatypes import sniff
from galaxy.datatypes.binary import *
from galaxy.datatypes.images import Pdf
from galaxy.datatypes.registry import Registry
from galaxy import util
from galaxy.util.json import *

try:
    import bz2
except:
    bz2 = None

assert sys.version_info[:2] >= ( 2, 4 )

def stop_err( msg, ret=1 ):
    sys.stderr.write( msg )
    sys.exit( ret )
def file_err( msg, dataset, json_file ):
    json_file.write( to_json_string( dict( type = 'dataset',
                                           ext = 'data',
                                           dataset_id = dataset.dataset_id,
                                           stderr = msg ) ) + "\n" )
    # never remove a server-side upload
    if dataset.type in ( 'server_dir', 'path_paste' ):
        return
    try:
        os.remove( dataset.path )
    except:
        pass
def safe_dict(d):
    """
    Recursively clone json structure with UTF-8 dictionary keys
    http://mellowmachines.com/blog/2009/06/exploding-dictionary-with-unicode-keys-as-python-arguments/
    """
    if isinstance(d, dict):
        return dict([(k.encode('utf-8'), safe_dict(v)) for k,v in d.iteritems()])
    elif isinstance(d, list):
        return [safe_dict(x) for x in d]
    else:
        return d
def check_html( temp_name, chunk=None ):
    if chunk is None:
        temp = open(temp_name, "U")
    else:
        temp = chunk
    regexp1 = re.compile( "<A\s+[^>]*HREF[^>]+>", re.I )
    regexp2 = re.compile( "<IFRAME[^>]*>", re.I )
    regexp3 = re.compile( "<FRAMESET[^>]*>", re.I )
    regexp4 = re.compile( "<META[^>]*>", re.I )
    regexp5 = re.compile( "<SCRIPT[^>]*>", re.I )
    lineno = 0
    for line in temp:
        lineno += 1
        matches = regexp1.search( line ) or regexp2.search( line ) or regexp3.search( line ) or regexp4.search( line ) or regexp5.search( line )
        if matches:
            if chunk is None:
                temp.close()
            return True
        if lineno > 100:
            break
    if chunk is None:
        temp.close()
    return False
def check_binary( temp_name ):
    is_binary = False
    temp = open( temp_name, "U" )
    chars_read = 0
    for chars in temp:
        for char in chars:
            chars_read += 1
            if ord( char ) > 128:
                is_binary = True
                break
            if chars_read > 100:
                break
        if chars_read > 100:
            break
    temp.close()
    return is_binary
def check_bam( temp_name ):
    return Bam().sniff( temp_name )
def check_sff( temp_name ):
    return Sff().sniff( temp_name )
def check_pdf( temp_name ):
    return Pdf().sniff( temp_name )
def check_bigwig( temp_name ):
    return BigWig().sniff( temp_name )
def check_bigbed( temp_name ):
    return BigBed().sniff( temp_name )
def check_gzip( temp_name ):
    # This method returns a tuple of booleans representing ( is_gzipped, is_valid )
    # Make sure we have a gzipped file
    try:
        temp = open( temp_name, "U" )
        magic_check = temp.read( 2 )
        temp.close()
        if magic_check != util.gzip_magic:
            return ( False, False )
    except:
        return ( False, False )
    # We support some binary data types, so check if the compressed binary file is valid
    # If the file is Bam, it should already have been detected as such, so we'll just check
    # for sff format.
    try:
        header = gzip.open( temp_name ).read(4)
        if binascii.b2a_hex( header ) == binascii.hexlify( '.sff' ):
            return ( True, True )
    except:
        return( False, False )
    CHUNK_SIZE = 2**15 # 32Kb
    gzipped_file = gzip.GzipFile( temp_name, mode='rb' )
    chunk = gzipped_file.read( CHUNK_SIZE )
    gzipped_file.close()
    # See if we have a compressed HTML file
    if check_html( temp_name, chunk=chunk ):
        return ( True, False )
    return ( True, True )
def check_bz2( temp_name ):
    try:
        temp = open( temp_name, "U" )
        magic_check = temp.read( 3 )
        temp.close()
        if magic_check != util.bz2_magic:
            return ( False, False )
    except:
        return( False, False )
    CHUNK_SIZE = 2**15 # reKb
    bzipped_file = bz2.BZ2File( temp_name, mode='rb' )
    chunk = bzipped_file.read( CHUNK_SIZE )
    bzipped_file.close()
    # See if we have a compressed HTML file
    if check_html( temp_name, chunk=chunk ):
        return ( True, False )
    return ( True, True )
def check_zip( temp_name ):
    if zipfile.is_zipfile( temp_name ):
        return True
    return False
def parse_outputs( args ):
    rval = {}
    for arg in args:
        id, files_path, path = arg.split( ':', 2 )
        rval[int( id )] = ( path, files_path )
    return rval
def add_file( dataset, registry, json_file, output_path ):
    data_type = None
    line_count = None
    converted_path = None
    stdout = None
    link_data_only = dataset.get( 'link_data_only', 'copy_files' )

    try:
        ext = dataset.file_type
    except AttributeError:
        file_err( 'Unable to process uploaded file, missing file_type parameter.', dataset, json_file )
        return

    if dataset.type == 'url':
        try:
            temp_name, dataset.is_multi_byte = sniff.stream_to_file( urllib.urlopen( dataset.path ), prefix='url_paste' )
        except Exception, e:
            file_err( 'Unable to fetch %s\n%s' % ( dataset.path, str( e ) ), dataset, json_file )
            return
        dataset.path = temp_name
    # See if we have an empty file
    if not os.path.exists( dataset.path ):
        file_err( 'Uploaded temporary file (%s) does not exist.' % dataset.path, dataset, json_file )
        return
    if not os.path.getsize( dataset.path ) > 0:
        file_err( 'The uploaded file is empty', dataset, json_file )
        return
    if not dataset.type == 'url':
        # Already set is_multi_byte above if type == 'url'
        try:
            dataset.is_multi_byte = util.is_multi_byte( codecs.open( dataset.path, 'r', 'utf-8' ).read( 100 ) )
        except UnicodeDecodeError, e:
            dataset.is_multi_byte = False
    if not data_type:
        # See if we have a zip archive
        is_zipped = check_zip( dataset.path )
        if is_zipped:
            data_type = 'zip'
    if not data_type:
        if check_binary( dataset.path ):
            # We have a binary dataset, but it is not Bam, Sff or Pdf
            data_type = 'binary'
            ext = "dat"
            
    datatype = registry.get_datatype_by_extension( ext )
    shutil.move( dataset.path, output_path )

    # Write the job info
    stdout = stdout or 'uploaded %s file' % data_type
    info = dict( type = 'dataset',
                 dataset_id = dataset.dataset_id,
                 ext = ext,
                 stdout = stdout,
                 name = dataset.name,
                 line_count = line_count )
    json_file.write( to_json_string( info ) + "\n" )
    if datatype.dataset_content_needs_grooming( output_path ):
        # Groom the dataset content if necessary
        datatype.groom_dataset_content( output_path )

def add_composite_file( dataset, registry, json_file, output_path, files_path ):
        if dataset.composite_files:
            os.mkdir( files_path )
            for name, value in dataset.composite_files.iteritems():
                value = util.bunch.Bunch( **value )
                if dataset.composite_file_paths[ value.name ] is None and not value.optional:
                    file_err( 'A required composite data file was not provided (%s)' % name, dataset, json_file )
                    break
                elif dataset.composite_file_paths[value.name] is not None:
                    dp = dataset.composite_file_paths[value.name][ 'path' ]
                    isurl = dp.find('://') <> -1 # todo fixme
                    if isurl:
                       try:
                           temp_name, dataset.is_multi_byte = sniff.stream_to_file( urllib.urlopen( dp ), prefix='url_paste' )
                       except Exception, e:
                           file_err( 'Unable to fetch %s\n%s' % ( dp, str( e ) ), dataset, json_file )
                           return
                       dataset.path = temp_name
                       dp = temp_name
                    if not value.is_binary:
                        if dataset.composite_file_paths[ value.name ].get( 'space_to_tab', value.space_to_tab ):
                            sniff.convert_newlines_sep2tabs( dp )
                        else:
                            sniff.convert_newlines( dp )
                    shutil.move( dp, os.path.join( files_path, name ) )
        # Move the dataset to its "real" path
        shutil.move( dataset.primary_file, output_path )
        # Write the job info
        info = dict( type = 'dataset',
                     dataset_id = dataset.dataset_id,
                     stdout = 'uploaded %s file' % dataset.file_type )
        json_file.write( to_json_string( info ) + "\n" )

def __main__():

    if len( sys.argv ) < 4:
        print >>sys.stderr, 'usage: upload.py <root> <datatypes_conf> <json paramfile> <output spec> ...'
        sys.exit( 1 )

    output_paths = parse_outputs( sys.argv[4:] )
    json_file = open( 'galaxy.json', 'w' )

    registry = Registry( sys.argv[1], sys.argv[2] )

    for line in open( sys.argv[3], 'r' ):
        dataset = from_json_string( line )
        dataset = util.bunch.Bunch( **safe_dict( dataset ) )
        try:
            output_path = output_paths[int( dataset.dataset_id )][0]
        except:
            print >>sys.stderr, 'Output path for dataset %s not found on command line' % dataset.dataset_id
            sys.exit( 1 )
        if dataset.type == 'composite':
            files_path = output_paths[int( dataset.dataset_id )][1]
            add_composite_file( dataset, registry, json_file, output_path, files_path )
        else:
            add_file( dataset, registry, json_file, output_path )
    # clean up paramfile
    try:
        os.remove( sys.argv[1] )
    except:
        pass

if __name__ == '__main__':
    __main__()