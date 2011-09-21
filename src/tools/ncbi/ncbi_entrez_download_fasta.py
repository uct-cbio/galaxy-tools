#!/usr/bin/python

# Retrieve FASTA sequence from NCBI based on an Entrez query.

# to start a new download
# entrez_download_fasta.py -m test@test.com -d nucleotide -e "Cypripedioideae[Orgn] AND matK[Gene]" -f out.fasta -i download_info.txt

# to continue a previous download
# entrez_download_fasta.py -m test@test.com -d nucleotide -c -w "NCID_1_38065753_130.14.18.52_9001_1300878409_78339627" -q 1 -s 11 -n 38 -f out_continue.fasta -i info_continue.txt

import sys
import os
import string
import re
from optparse import OptionParser
from Bio import Entrez

def main():

    usage = "usage: %prog -m EMAIL -d DATABASE -e ENTREZ_QUERY -c -w WEB_ENV -q QUERY_KEY -s START_POS -n NR_SEQ_QUERY -f FASTA_FILE -i INFO_FILE"
    parser = OptionParser(usage=usage)
    parser.add_option("-m", "--email", dest="email", help="Email address. Need to provide this when doing an Entrez search and fetch")
    parser.add_option("-d", "--database", dest="database", help="Database e.g. nucleotide")
    parser.add_option("-e", "--entrez_query", dest="entrez_query", help="Entrez query e.g. \"Bacteria\"[Organism] OR \"Archaea\"[Organism] OR prokaryotes[All Fields] not \"Escherichia coli\"[Organism]")
    parser.add_option("-c", "--continue_download", action="store_true", dest="continue_download", default=False, help="If flag is specified program will continue a previous download. User need to provide the WEB_ENV, QUERY_KEY and SEQ_START")
    parser.add_option("-w", "--web_env", dest="web_env", help="Please provide the previous web_env.")
    parser.add_option("-q", "--query_key", dest="query_key", help="Please provide the previous query_key.")
    parser.add_option("-s", "--start_pos", dest="start_pos", help="Please provide position of the sequence to start downloading from. E.g. the position where the previous download failed.")
    parser.add_option("-n", "--nr_seq_query", dest="nr_seq_query", help="Please provide the number of sequences found in the original query.")
    parser.add_option("-f", "--fasta_file", dest="fasta_file", help="FASTA output file")
    parser.add_option("-i", "--info_file", dest="info_file", help="Information related to the download. Contains the web_env and query_key for the search.")
    
    (options, args) = parser.parse_args()
    
    if not options.continue_download:
        print "Start a new download..."
        if not options.email:
            print "Please specify an email address. Need to provide this when doing an Entrez search or fetch (-m EMAIL)"
            return - 1
        if not options.database:
            print "Please specify the database to fetch info from (-d DATABASE)"
            return - 2
        if not options.entrez_query:
            print "Please specify an entrez query (-e ENTREZ_QUERY)"
            return - 3
        if not options.fasta_file:
            print "Please specify the FASTA output file (-f FASTA_FILE)"
            return - 4
        if not options.info_file:
            print "Please specify the download info file (-i INFO_FILE)"
            return - 5
         # Need to to some checking on the on the length of the arguments provided. Currently not working because we have 2 options
         # (1) Start new download and (2) Continue previous download. Need to handle these separately. Sort it out later. Not to crucial.
    else:
        print "Continue a previous download..."
        if not options.email:
            print "Please specify an email address. Need to provide this when doing an Entrez search or fetch (-m EMAIL)"
            return - 6
        if not options.database:
            print "Please specify the database to fetch info from (-d DATABASE)"
            return - 7
        if not options.web_env:
            print "Please specify the previous web_env (-w WEB_ENV)"
            return - 8
        if not options.query_key:
            print "Please specify the previous query_key (-q QUERY_KEY)"
            return - 9
        if not options.start_pos:
            print "Please specify the position of the sequence to start downloading from (-s START_POS)"
            return - 10
        if not options.nr_seq_query:
            print "Please specify the number of sequences in original query (-n NR_SEQ_QUERY)"
            return - 11
        if not options.fasta_file:
            print "Please specify the FASTA output file (-f FASTA_FILE)"
            return - 12
        if not options.info_file:
            print "Please specify the download info file (-i INFO_FILE)"
            return - 13

    if (len(args) > 0):
        print "Too many arguments"
        return - 14

    BATCH_SIZE = 100
    
    # Input strings generated by browser/galaxy needs to be replaced
    mapped_chars = { '>' :'__gt__',
                 '<' :'__lt__',
                 '\'' :'__sq__',
                 '"' :'__dq__',
                 '[' :'__ob__',
                 ']' :'__cb__',
         '{' :'__oc__',
                 '}' :'__cc__',
                 '@' :'__at__',
    }

    # Start a new download
    if(not options.continue_download):
        email = options.email
        database = options.database
        entrez_query = options.entrez_query
        fasta_file = options.fasta_file
        info_file = options.info_file
        
        for key, value in mapped_chars.items():
            database = database.replace(value, key)
            entrez_query = entrez_query.replace(value, key)
            email = email.replace(value, key)
        
        Entrez.email = email
        
        # Open info_file for writing
        info_file_fd = open(info_file, "w")
        info_file_fd.write('Email address: %s\n' % Entrez.email)
        info_file_fd.write('Database: %s\n' % database)
        info_file_fd.write('Entrez query: %s\n' % entrez_query)
        try:
            handle = Entrez.esearch(db=database,term=entrez_query, usehistory='y')
            results = Entrez.read(handle)
            handle.close()
        except Exception, e:
            info_file_fd.write( "Error raised when trying do an Entrez.esearch: %s\n" % str(e))
            info_file_fd.close()
            sys.exit( "Error raised! Exiting now!")
            
        gi_list = results["IdList"]
        nr_seq_query = int(results["Count"])
    #    assert count == len(gi_list) # do not do this test gi_list it is always 20
        # Get web_env and query_key from search results
        web_env = results["WebEnv"]
        query_key = results["QueryKey"]
        
        # Write query specific info to info_file 
        info_file_fd.write('Number of sequences in query: %d\n' % nr_seq_query)
        info_file_fd.write('Number of sequences to be dowloaded: %d\n' % nr_seq_query)
        info_file_fd.write('Download sequences from position: %d\n' % 0)
        info_file_fd.write('web_env: %s\n' % web_env)
        info_file_fd.write('query_key: %s\n'% query_key)
        info_file_fd.write('Downloading sequences in batches of %d\n' % BATCH_SIZE)
        info_file_fd.close()
        
        # Now retrieve the FASTA sequences in batches of 5
        fasta_file_fd = open(fasta_file, "w")
        for start in range(0,nr_seq_query, BATCH_SIZE):
            end = min(nr_seq_query, start + BATCH_SIZE)
            print "Dowloading sequence %i to %i" % (start+1, end)
            try:
                fetch_handle = Entrez.efetch(db=database, rettype="fasta", retstart=start, retmax=BATCH_SIZE, webenv=web_env, query_key=query_key)
                data = fetch_handle.read()
            except Exception, e:
                info_file_fd = open(info_file, "a")
                info_file_fd.write( "Error raised when trying do an Entrez.efind: %s\n" % str(e))
                info_file_fd.close()
                sys.exit( "Error raised! Exiting now!")
            fetch_handle.close()
            fasta_file_fd.write(data)
            info_file_fd = open(info_file, "a")
            info_file_fd.write('Downloaded sequence %i to %i\n' % (start+1, end))
            info_file_fd.close()
            
        fasta_file_fd.close()
        
    else: # Continue a previous download
        email = options.email
        database = options.database
        web_env = options.web_env
        query_key = options.query_key
        start_pos = int(options.start_pos); # should check if start position is a integer
        nr_seq_query = int(options.nr_seq_query); # should check if nr_seq_query is a integer
        fasta_file = options.fasta_file
        info_file = options.info_file
        
        for key, value in mapped_chars.items():
            database = database.replace(value, key)
            email = email.replace(value, key)
            web_env = web_env.replace(value, key)
            query_key = query_key.replace(value, key)
        
        Entrez.email = email
        
        # Open info_file for writing
        info_file_fd = open(info_file, "w")
        info_file_fd.write('Email address: %s\n' % Entrez.email)
        info_file_fd.write('Database: %s\n' % database)
         # Write query specific info to info_file
        info_file_fd.write('Number of sequences in original query: %d\n' % nr_seq_query)
        info_file_fd.write('Number of sequences to be dowloaded: %d\n' % (int(nr_seq_query) - int(start_pos) + 1))
        info_file_fd.write('Download sequences from position: %d\n' % start_pos)
        info_file_fd.write('web_env: %s\n' % web_env)
        info_file_fd.write('query_key: %s\n' % query_key)
        info_file_fd.write('Downloading sequences in batches of %d\n' %  BATCH_SIZE)
        info_file_fd.close()
        info_file_fd = open(info_file, "a")

        # Now retrieve the FASTA sequences in batches of BATCH_SIZE
        fasta_file_fd = open(fasta_file, "w")
        for start in range(start_pos - 1,nr_seq_query, BATCH_SIZE):
            end = min(nr_seq_query, start + BATCH_SIZE)
            print "Dowloading sequence %i to %i" % (start+1, end)
            try:
                fetch_handle = Entrez.efetch(db=database, rettype="fasta", retstart=start, retmax=BATCH_SIZE, webenv=web_env, query_key=query_key)
                data = fetch_handle.read()
            except Exception, e:
                info_file_fd = open(info_file, "a")
                info_file_fd.write( "Error raised when trying do an Entrez.efind: %s\n" % str(e))
                info_file_fd.write( "Retrying...")
                info_file_fd.close()
#                sys.exit( "Error raised! Exiting now!") # do not exit anymore
            fetch_handle.close()
            fasta_file_fd.write(data)
            info_file_fd = open(info_file, "a")
            info_file_fd.write('Downloaded sequence %i to %i\n' % (start+1, end))
            info_file_fd.close()
        
        fasta_file_fd.close()
    
if __name__ == "__main__":
    sys.exit(main())

  
