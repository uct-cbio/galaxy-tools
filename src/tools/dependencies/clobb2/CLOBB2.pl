#!/usr/bin/perl -w

##############################################################################
##############################################################################
##############################################################################
# CLOBB2.pl - the 'new improved' clustering algorithm
# Usage CLOBB2.pl <cluster_name>

# Requisities - a directory named 'sequences' which
# contains fasta files of individual sequences. Other
# directories and files are created as needed by the
# program. Note cluster name is typically a 3 letter
# identifier (e.g. DRC - for danio rerio clusters).

# You may also have to change the location of the
# ncbi blast executables in the lines :

# system ("megablast -d $DATABASE -i $DATABASE.$n -e 1 -W 16 -v 1000 -b 1000 -F F -D 2 -o $rootdir/$CLUSTERNAME$n.out >> $log.$s");
# system ("blastall -p blastn -d $DATABASE -i $DATABASE.$n -G 2 -e 1 -v 100 -F F -o $rootdir/$CLUSTERNAME$n.out >> $log.$s");

# CLOBB produces the following output :
# <cluster_name>EST - a fasta database of the sequences
# with cluster assignments - <cluster_name>EST00001 -> <cluster_name>EST99999
# A directory 'OUT' which contains a logfile detailing the treatment of each sequence.
# A file called 'merge' which contains a list of all the clusters
# which have been merged into another cluster.
# A file called 'superclusters' which lists related clusters
#
#
##############################################################################
# Last Updated: December 17, 2003. by Andrew Yam
# Improvements :
# Use of a single BLAST (no db reformatting required)
# Adoption of MegaBLAST
# (although blastall
# Incorporation of strict :-)
# Benchmarking has also been included
# CLOBB now behaves pseudo-linearly w.r.t time taken to run
# and is upto 10 times faster than the previous version
#
# some modifications by Ralf Schmid (R.Schmid@ed.ac.uk)
# pattern match to catch AB123456 format, trimming down of logfile (diskspace),
# rearrange blast and blast parsing (diskspace), removed Time:HiRes stuff (not in
# Bio-Linux), added some error handling
#
##############################################################################
##############################################################################
##############################################################################
# Contributors :
# John Parkinson (jparkin@sickids.ca)
# Mark Blaxter (mark.blaxter@ed.ac.uk)
# Andrew Yam (andrew_yam@hotmail.com)
# David Guiliano (d.guiliano@ucl.ac.uk)

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License Version 2
# as published by the Free Software Foundation.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

use File::stat;
use Bio::SeqIO;
use strict;

# check for argument
if ( !@ARGV ) {
    print "\nERROR: No argument provided\n\n";
    print "USAGE: CLOBB2.pl [ABC]\n";
    print "The argument is required as cluster identifier\n\n";
    exit;
}

# Variables initialization
our ( $CLUSTERNAME, $DATABASE, $DATABASE_LOCATION, $INDEX_NUM );
our ( $query_accession, $sl_flag,   $seq_count );
our ( $rootdir,         $seq_files, $seq_files_done );
our (%clusID_table);
our ($clustername);
my ( @s_length, @s_start, @s_end, @q_start, @q_end,    @chimname, @sbjfile, @clustermatch );
my ( $line,     $in,      $seq,   $file,    $filesize, $dbsize,   $q_files, $NoHit );
my ( $i,        $n,       $k,     $kk,      $chimera,  $local,    $tmp,     $count, $string );
my ( $accession, $query_length, $sbjct_length, $homology, $gaps, $total );
my ( $evalue, $num_match, $identity, $query_strand, $sbjct_strand, $query_start, $query_end, $sbjct_start, $sbjct_end );
my ( $match, $mismatch, $percent );
my ( $overlap_length, $overlap_error );
my ( $slarger, $ssmaller, $qlarger, $qsmaller );
my ( $start_time, $end_time, $elapsed );

$CLUSTERNAME       = $ARGV[0];
$DATABASE          = $CLUSTERNAME . "EST";
$DATABASE_LOCATION = "./$DATABASE";
$|                 = 1;                      # make STDOUT print buffer flush immediately

$rootdir        = "./OUT";
$seq_files      = "./sequences";             # Directory to find the sequence files
$seq_files_done = "./sequences_done";        # Directory to place used sequence files

# Set up the directories and relevant files if not already present
$INDEX_NUM = 0;
if ( !stat("sequences") ) {
    print "DIRECTORY \"sequences\" NOT FOUND!\n";
    print "Please make sure you are in the right directory.\n";
    exit;
}
if ( !stat("OUT") ) { system("mkdir OUT"); }
else {
    if ( -e "OUT/*out" ) { system("rm OUT/*.out"); }
}
if ( !stat("$seq_files_done") ) { system("mkdir $seq_files_done"); }
if ( !stat("superclusters") )   { system("touch superclusters"); }
if ( !stat("merge") )           { system("touch merge"); }
if ( stat("database") )         { system("rm database"); }

#Deal with common file ending .fsa
my @renamed = glob("$seq_files/*.fsa");
foreach my $rename (@renamed) {
    $rename =~ s/\.fsa//g;
    system("mv $rename.fsa $rename");
}

# Find the old INDNUMBER from the database if database exists
# Store the cluster ID in a hash table corresponding to the sequence.
if ( !stat($DATABASE) ) {
    system("touch $DATABASE");
} else {
    $in = Bio::SeqIO->new( -file => "$DATABASE", '-format' => 'Fasta' );
    while ( $line = $in->next_seq() ) {
        $seq = $line->display_id();
        $seq =~ s/^.+\|([A-Z]+\d+).*$/$1/;
        $clusID_table{$seq} = "";
        $clustername = $line->desc();
        $clustername =~ /$CLUSTERNAME(\d+)/;
        $clusID_table{$seq} = $CLUSTERNAME . $1;

        if ( $1 > $INDEX_NUM ) {
            $INDEX_NUM = $1;
            $INDEX_NUM =~ s/^0+//;
        }
    }
    system("cp $DATABASE ./database");
}

# Generate DB containing all sequences to be blasted.
# Query files are splitted by 1000 sequences per file.
opendir( SEQDIR, "$seq_files" );
$seq_count = 0;
$q_files   = 1;
system("touch $DATABASE.$q_files");
my $my_count = 0;

while ( defined( $file = readdir(SEQDIR) ) ) {
    $my_count = $my_count + 1;
    if ( $file !~ /^\./ ) {

        $filesize = stat("$seq_files/$file")->size;

        if ( $filesize > 100 ) {    # Minimum size of file for clustering

            # Split up the database if the size of the DB exceeds 1000 sequences.
            $seq_count++;

            if ( $seq_count % 1000 == 0 ) {
                $q_files++;
                system("touch $DATABASE.$q_files");
            }
            system "cat $seq_files/$file >> $DATABASE_LOCATION.$q_files";
            system "echo '\n' >> $DATABASE_LOCATION.$q_files";
            $clusID_table{$file} = "";
        } else {
            system("rm $seq_files/$file");
        }
    }
}
close(SEQDIR);

# Concatenating BLAST DB
for ( $n = 1 ; $n <= $q_files ; $n++ ) {
    system "cat $DATABASE_LOCATION.$n >> $DATABASE_LOCATION";
}

system "cp $DATABASE_LOCATION tmp_database";
system "formatdb -i $DATABASE_LOCATION -p F";

# Open a log file for the session.  Use the date to give each logfile a distinctive name.
my ( $min, $hr, $mday, $yr, $mon, $log, $s );
$hr   = "";
$min  = "";
$mday = "";
$mon  = "";
$yr   = "";
$s    = 1;

( $min, $hr, $mday, $yr ) = (localtime)[ 1, 2, 3, 5 ];
$mon = ( "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" )[ (localtime)[4] ];
$min = "0" . $min if $min < 10;
$hr  = "0" . $hr  if $hr < 10;
$log = "$rootdir/logfile_$hr:${min}_$mday$mon$yr";
open( LOG, ">$log.$s" );
$| = 1;

my $seqCount = 1;
print LOG "Blasting $CLUSTERNAME\n";

# BLAST the query against each of the splitted DB.
# changed order to BLAST -> parse -> rm blasts -> nextBLAST... (again disk space issue)
for ( $n = 1 ; $n <= $q_files ; $n++ ) {

#    system(
#"megablast -d $DATABASE -i $DATABASE.$n -e 1 -W 16 -v 1000 -b 1000 -F F -D 2 -o $rootdir/$CLUSTERNAME$n.out >> $log.$s"
#    );

    # lower case regions will be filterred in blast (e.g can be used when low complexity regions are flagged/masked)
    system(
"megablast -d $DATABASE -i $DATABASE.$n -e 1 -W 16 -v 1000 -b 1000 -F F -U T -D 2 -o $rootdir/$CLUSTERNAME$n.out >> $log.$s"
    );

    &blast_parsing($n);
    system "rm $rootdir/$CLUSTERNAME$n.out";
}

print LOG "\nOK\n";
system("rm $DATABASE.*");

opendir( SEQDIR, "$seq_files" );
while ( defined( $file = readdir(SEQDIR) ) ) {
    if ( $file !~ /^\./ ) {
        $clustername = $CLUSTERNAME . &index( ++$INDEX_NUM );

        #strip file endings, otherwise things go berserk
        $file =~ s/\.fsa//g;
        $file =~ s/\.seq//g;
        &newfile( $file, $clustername );
        $clusID_table{$file} = $clustername;
        print LOG "$file assigned $clustername\n";
        system "cat $seq_files_done/$file >> ./database";
        system "echo '\n' >> ./database";
        $seq_count--;
        printf( "\r%5d sequences remaining.", $seq_count );
    }
}
closedir(SEQDIR);

system("mv ./database $DATABASE_LOCATION");
print LOG "the current $CLUSTERNAME is $INDEX_NUM\n";
close(LOG);

##############################################################################
################################  subroutines ################################
##############################################################################

#############################################################################################
# Parsing the BLAST output file.
#############################################################################################
sub blast_parsing {

    my $i = $_[0];

    #for ($i = 1; $i <= $q_files; $i++) {

    open( BLAST, "<$rootdir/$CLUSTERNAME$i.out" );

    $filesize = stat("$log.$s")->size;
    if ( $filesize > 100000000 ) {
        close(LOG);
        $| = 1;
        system("gzip $log.$s");
        $s++;
        open( LOG, ">$log.$s" );
    }

    $query_accession = "";
    $accession       = "";
    while ( my $line = <BLAST> ) {

        # If no hit found.....
        $NoHit = 0;
        if ( $line !~ /^\s*$/ ) {
            if ( $line =~ /No hits found/ ) {
                $NoHit = 1;
                last;
            }

            # For a new Result....,
            if ( $line =~ /^Query= / ) {
                if ( $line =~ /^Query=.*\|(\w+\d+)\s+.*/ || $line =~ /^Query= ([\w\.-]+)\s+.*/ ) {

                    # If all information is obtained, start cluster.
                    if ( $query_accession ne "" ) {

                        &DoCluster(
                                    $k,               $kk,            $chimera,
                                    $query_accession, \@clustermatch, \@chimname,
                                    \@sbjfile,        \@q_start,      \@q_end
                        );

                        # Assign the read to a cluster old or new.
                        &newfile( $query_accession, $clustername );
                        $clusID_table{$query_accession} = $clustername;
                        print LOG "$query_accession assigned $clustername\n";
                        system "cat $seq_files_done/$query_accession >> ./database";
                        system "echo '\n' >> ./database";

                        $seq_count--;
                        printf( "\r%5d sequences remaining.", $seq_count );
                    }

                    # Reset variables at each new result.
                    $query_accession = "";
                    $query_length    = 0;
                    $accession       = "";
                    $sbjct_length    = 0;
                    $evalue          = "";
                    $num_match       = 0;
                    $identity        = 0;
                    $gaps            = 0;
                    $homology        = "";
                    $query_strand    = "";
                    $sbjct_strand    = "";
                    $query_start     = 0;
                    $query_end       = 0;
                    $sbjct_start     = 0;
                    $sbjct_end       = 0;
                    $total           = 0;

                    $k               = 0;
                    $kk              = 0;
                    $s_length[0]     = 0;
                    $s_start[0]      = 0;
                    $s_end[0]        = 0;
                    $q_start[0]      = 0;
                    $q_end[0]        = 0;
                    $query_accession = $1;

                    print LOG "Query: $query_accession\n";
                }

                else { print "Couldn't extract Accession ID from $line"; exit; }    # try to catch queryname problems
            } elsif ( $line =~ /^\s+\((\d+)\sletters\)/ ) {
                $query_length = $1;                                                 #total query length;
            }

            # For each hit for the query...
            elsif ( $line =~ /^>/ ) {
                if ( $line =~ /^>gi.*\|(\w+\d+).*/ || $line =~ /^>([\w\.-]+).*/ ) {
                    $accession   = $1;
                    $clustername = $clusID_table{$accession};

                    # Reset variable at each new hit.
                    $sbjct_length = 0;
                    $evalue       = "";
                    $num_match    = 0;
                    $identity     = 0;
                    $gaps         = 0;
                    $query_strand = "";
                    $sbjct_strand = "";
                    $homology     = "";
                    $total        = 0;
                    $query_start  = 0;
                    $query_end    = 0;
                    $sbjct_start  = 0;
                    $sbjct_end    = 0;
                } else {
                    print "Couldn't extract Accession ID from $line";
                    exit;
                }    # try to catch queryname problems
            } elsif ( $line =~ /^\s+Length\s=\s(\d+)/ ) {
                $sbjct_length = $1;    #total subject length;
            }

            # If the subject sequence has previously assigned a clusterID, then proceed....
            if ( ( $accession ne "" ) && ( $query_accession ne "" ) ) {
                if ( ( $clusID_table{$accession} ne "" ) && ( $query_accession ne $accession ) ) {

                    # For each HSP for the hit...
                    if ( $line =~ /^\sScore.*Expect\s=\s(.*)/ ) {
                        $chimera      = 0;
                        $local        = 0;
                        $total        = 0;
                        $query_start  = 0;
                        $query_end    = 0;
                        $sbjct_start  = 0;
                        $sbjct_end    = 0;
                        $num_match    = 0;
                        $total        = 0;
                        $identity     = 0;
                        $gaps         = 0;
                        $query_strand = "";
                        $sbjct_strand = "";
                        $homology     = "";
                        $evalue       = $1;
                    }

                    elsif ( $line =~ /^\sIdentities\s=\s(\d+)\/(\d+)\s\((\d+)\%\).*$/ ) {
                        $num_match = $1;
                        $total     = $2;
                        $identity  = $3;
                        if ( $line =~ /Gaps\s=\s(\d+)\// ) {
                            $gaps = $1;
                        }
                    }

                    elsif ( $line =~ /^\sStrand\s=\s(\w+)\s\/\s(\w+)/ ) {
                        $query_strand = $1;
                        $sbjct_strand = $2;
                    }

                    elsif ( $line =~ /^Query:\s(\d+).*\s(\d+)$/ ) {

                        #          if ($query_strand eq "Plus") {
                        if ( $query_start == 0 ) {
                            $query_start = $1;
                        }
                        $query_end = $2;

                        #          }
                        #          else {
                        #            if ( $query_end == 0) {
                        #              $query_end = $1;
                        #            }
                        #            $query_start = $2;
                        #          }
                    }

                    elsif ( $line =~ /\|{2,}/ ) {
                        chomp($line);
                        $line =~ s/\s{11}//;
                        $homology .= $line;
                    }

                    elsif ( $line =~ /^Sbjct:\s(\d+).*\s(\d+)$/ ) {

                        #          if ($sbjct_strand eq "Plus") {
                        if ( $sbjct_start == 0 ) {
                            $sbjct_start = $1;
                        }
                        $sbjct_end = $2;

                        #          }
                        #          else {
                        #            if ( $sbjct_end == 0) {
                        #              $sbjct_end = $1;
                        #            }
                        #            $sbjct_start = $2;
                        #          }
                    }

                    if    ( $sbjct_end >= $sbjct_start ) { $slarger = $sbjct_end;   $ssmaller = $sbjct_start; }
                    elsif ( $sbjct_start > $sbjct_end )  { $slarger = $sbjct_start; $ssmaller = $sbjct_end; }
                    if    ( $query_end >= $query_start ) { $qlarger = $query_end;   $qsmaller = $query_start; }
                    elsif ( $sbjct_start > $sbjct_end )  { $qlarger = $query_start; $qsmaller = $query_end; }

                    # The condition that satisfies the end of each HSP.
                    if ( ( $total * 2 ) == ( ( $slarger - $ssmaller + 1 ) + ( $qlarger - $qsmaller + 1 ) + $gaps ) ) {
                        if ( $query_start < 10 && $sbjct_start < 10 ) {
                            $homology =~ s/\s{10}//;
                        }

#print "QA:$query_accession QL:$query_length SA:$accession SL:$sbjct_length E:$evalue NM:$num_match T:$total I:$identity ";
#print "QS:$query_strand SS:$sbjct_strand QST:$query_start QE:$query_end SST:$sbjct_start SE:$sbjct_end\n";
#print "H:$homology\n";

                        # First check the sequence that was read in for low 'Identity', high E value hits
                        if ( $num_match > 30 ) {

                          # The E value score is high enough to warrant looking for a region of high % identity.
                          # If a stretch of 30 bases with > 95% identity is found then regard this HSP as a 'pukka' hit.
                            if ( ( $evalue =~ /e-(\d+)/ || $evalue == 0.0 ) && $identity <= 94 ) {
                                $count  = 0;
                                $string = "";

                                while ( $homology =~ /(.)/g ) {
                                    $string .= $1;
                                    $count++;
                                    if ( $count > 30 ) {
                                        $string =~ s/.//;
                                        $count--;
                                    }

                                    if ( $count == 30 ) {
                                        $match    = 0;
                                        $mismatch = 0;

                                        while ( $string =~ /(.)/g ) {
                                            if    ( $1 eq "|" ) { $match++; }
                                            elsif ( $1 eq " " ) { $mismatch++; }
                                        }
                                        $percent = $match / ( $match + $mismatch );
                                        if ( $percent > 0.94 ) {
                                            $local = 1;
                                            last;
                                        }
                                    }
                                }
                            }

                            if ( $identity > 94 || $local == 1 ) {

                              # We have a match of sorts now check to see if seqs overlap without matching properly.
                              # NB If accepted one HSP, must ignore next HSPs which may be from the same subject.
                              # This is acheived using the chimera flag.
                              # An overlap of 10% of the query/subject length is regarded as OK (since the start/ends of
                              # reads tends to be low quality.
                              # If the overlap is larger - then number of N's in the overlap are counted.
                              # See count_ns sub routine for criteria for regarding overlaps as pukka.
                                if ( $chimera == 0 ) {
                                    $sl_flag = 0;

                                    if ( $sbjct_start < $sbjct_end ) {    # plus strand
                                        $overlap_length = $query_length - $query_end;
                                        $overlap_error  = $query_length / 10;

                                        if ( ( $sbjct_length - $sbjct_end ) < $overlap_length ) {
                                            $overlap_length = $sbjct_length - $sbjct_end;
                                            $overlap_error  = $sbjct_length / 10;
                                        }

                                        if ( $overlap_length > $overlap_error ) {
                                            $sl_flag = 1;
                                            &count_ns( $overlap_length, $sbjct_end, $query_end, 1, 1, $accession );
                                        }

                                        $overlap_length = $query_start;
                                        $overlap_error  = $query_length / 10;

                                        if ( $sbjct_start < $overlap_length ) {
                                            $overlap_length = $sbjct_start;
                                            $overlap_error  = $sbjct_length / 10;
                                        }

                                        if ( $overlap_length > $overlap_error ) {
                                            $sl_flag = 1;
                                            &count_ns( $overlap_length, $sbjct_start, $query_start, -1, -1,
                                                       $accession );
                                        }
                                    }

                                    if ( $sbjct_start > $sbjct_end ) {    # minus strand
                                        $overlap_length = $query_length - $query_end;
                                        $overlap_error  = $query_length / 10;

                                        if ( $sbjct_end < $overlap_length ) {
                                            $overlap_length = $sbjct_end;
                                            $overlap_error  = $sbjct_length / 10;
                                        }

                                        if ( $overlap_length > $overlap_error ) {
                                            $sl_flag = 1;
                                            &count_ns( $overlap_length, $sbjct_end, $query_end, -1, 1, $accession );
                                        }

                                        $overlap_length = $query_start;
                                        $overlap_error  = $query_length / 10;

                                        if ( ( $sbjct_length - $sbjct_start ) < $overlap_length ) {
                                            $overlap_length = $sbjct_length - $sbjct_start;
                                            $overlap_error  = $sbjct_length / 10;
                                        }

                                        if ( $overlap_length > $overlap_error ) {
                                            $sl_flag = 1;
                                            &count_ns( $overlap_length, $sbjct_start, $query_start, 1, -1, $accession );
                                        }
                                    }

                                # Its the middle which overlaps and it should be placed in a new cluster.
                                # Can also tag it to show it has similarity to reject other matches to the same cluster.
                                    if ( $sl_flag == 1 ) {
                                        $chimera = 1;
                                        print( LOG "$query_accession is similar to $clustername ($accession)\n" );
                                        $chimname[$kk] = $clustername;
                                        $kk++;
                                    }

                                    # The overlap looks genuine and should be assigned to the same cluster.
                                    if ( $sl_flag == 0 ) {
                                        print( LOG "$query_accession appears part of $clustername ($accession)\n" );
                                        $chimera          = 2;
                                        $s_length[$k]     = $sbjct_length;
                                        $s_start[$k]      = $sbjct_start;
                                        $s_end[$k]        = $sbjct_end;
                                        $q_start[$k]      = $query_start;
                                        $q_end[$k]        = $query_end;
                                        $clustermatch[$k] = $clustername;
                                        $sbjfile[$k]      = $accession;
                                        $k++;
                                    }
                                }    # End of loop dealing with a pukka sequence
                            }
                        }
                    }
                }
            }    # End of section checking high E value HSPs for %identity
        }
    }

    if ( $NoHit == 0 ) {
        &DoCluster( $k, $kk, $chimera, $query_accession, \@clustermatch, \@chimname, \@sbjfile, \@q_start, \@q_end );

        # Assign the read to a cluster old or new.
        &newfile( $query_accession, $clustername );
        $clusID_table{$query_accession} = $clustername;
        print LOG "$query_accession assigned $clustername";
        system "cat $seq_files_done/$query_accession >> ./database";
        system "echo '\n' >> ./database";

        $seq_count--;
        printf( "\r%5d sequences remaining.", $seq_count );

    }
}

########################################################################################################################
########################################################################################################################
sub DoCluster {

    my ( $k, $kk, $chimera, $query_accession, $clustermatch, $chimname, $sbjfile, $q_start, $q_end ) = @_;
    my ( $merge_flag, $pukka_name_flag, $cluster_num, $low_cluster, $chim_flag, $chim_cluster );
    my ( $n, $m, $super_line, $new_line, $line_no );

    $merge_flag = 1;
    if ( $k == 0 ) {    # never got any matches
        $clustername = $CLUSTERNAME . &index( ++$INDEX_NUM );
        $clusID_table{$query_accession} = $clustername;
    }

    elsif ( $k == 1 ) {    # got a single match
        $clustername = @$clustermatch[0];
        if ( $kk > 0 ) {
            $pukka_name_flag = 1;
            for ( $n = 0 ; $n < $kk ; $n++ ) {
                print LOG "CM: @$clustermatch[0]; @$chimname[$n]\n";
                if ( @$clustermatch[0] eq @$chimname[$n] ) {
                    $pukka_name_flag = 0;
                    last;
                }
            }

            if ( $pukka_name_flag == 0 ) {
                $clustername = $CLUSTERNAME . &index( ++$INDEX_NUM );
            }
        }
    }

    elsif ( $k > 1 ) {    # got more than one match
                          # Find lowest cluster number.
        $cluster_num = 10000000;
        for ( $n = 0 ; $n < $k ; $n++ ) {
            @$clustermatch[$n] =~ /$CLUSTERNAME(\d+)/;
            if ( $cluster_num > $1 ) {
                $low_cluster = $n;
                $cluster_num = $1;
            }
        }

        # Need to find out if the two subject hits overlap by > 30,
        # if yes then reject the read -> problem read.
        # Could actually treat this as a new cluster - or stick in cluster with highest identity/E value ?
        # Otherwise merge the clusters
        for ( $n = 0 ; $n < $k ; $n++ ) {
            for ( $m = $n + 1 ; $m < $k ; $m++ ) {
                if ( @$clustermatch[$n] ne @$clustermatch[$m] ) {
                    if ( @$q_start[$n] == @$q_start[$m] || @$q_end[$n] == @$q_end[$m] ) { $merge_flag = 0; }
                    if ( @$q_start[$n] > @$q_start[$m] && @$q_start[$n] + 30 < @$q_end[$m] ) { $merge_flag = 0; }
                    if ( @$q_start[$m] > @$q_start[$n] && @$q_start[$m] + 30 < @$q_end[$n] ) { $merge_flag = 0; }
                }
            }

            # The next section marks a cluster as a problem since the read matched part of a cluster well
            # but another part of the same cluster badly, possibly an alternative splice ?
            # Also these clusters should not therefore be merged
            if ( $kk > 0 ) {
                for ( $m = 0 ; $m < $kk ; $m++ ) {
                    if ( @$clustermatch[$n] eq @$chimname[$m] ) {
                        $chimera    = 3;
                        $merge_flag = 0;
                    }
                }
            }
        }

        # This read hits > 1 clusters which are obviously different clusters,
        # so we'll give this read the same cluster number as the highest hit,
        # provided that there weren't any clashes.
        if ( $merge_flag == 0 ) {
            $clustername = $CLUSTERNAME . &index( ++$INDEX_NUM );

            for ( $n = 0 ; $n < $k ; $n++ ) {
                $pukka_name_flag = 1;

                for ( $m = 0 ; $m < $kk ; $m++ ) {
                    if ( @$clustermatch[$n] eq @$chimname[$m] ) {
                        $pukka_name_flag = 0;
                        last;
                    }
                }

                if ( $pukka_name_flag == 1 ) {
                    $clustername = @$clustermatch[$n];
                    last;
                } else {
                    @$clustermatch[$k] = $clustername;
                    $k++;
                }
            }

            if ( $pukka_name_flag == 0 ) {
                @$clustermatch[$k] = $clustername;
                $k++;
            }

            #&problem($query_accession);
            print LOG "***** $query_accession hits more than 1 cluster *****\n";
        }

        # merge the clusters !
        elsif ( $merge_flag == 1 ) {
            for ( $n = 0 ; $n < $k ; $n++ ) {
                if ( @$clustermatch[$n] ne @$clustermatch[$low_cluster] ) {
                    system("sed 's/@$clustermatch[$n]/@$clustermatch[$low_cluster]/g' ./database > temp_cluster");
                    system("mv temp_cluster ./database");
                    print LOG "Merging @$clustermatch[$n] into @$clustermatch[$low_cluster]\n";
                    open( MERGE, ">>merge" );
                    print MERGE "@$clustermatch[$n] => @$clustermatch[$low_cluster]\n";
                    close(MERGE);
                    $clusID_table{ @$sbjfile[$n] } = @$clustermatch[$low_cluster];
                }
            }
            $clustername = @$clustermatch[$low_cluster];
        }

        # The same cluster was hit several times, but was a 'bad' hit in at least one case
        if ( $chimera == 3 ) {
            open( CHIM_FILE, "superclusters" );
            $chim_flag = 0;

            while ( $super_line = <CHIM_FILE> ) {
                if ( $super_line =~ />@$clustermatch[$low_cluster]/ ) { $chim_flag = 1; last; }
                $super_line =~ />(.+)/;
                $chim_cluster = $1;
                if ( $super_line =~ /@$clustermatch[$low_cluster]/ ) { $chim_flag = 2; last; }
            }

            if ( $chim_flag == 2 ) {
                close(CHIM_FILE);
                open( CHIM_FILE, "superclusters" );
                $chim_flag = 0;
                while ( $super_line = <CHIM_FILE> ) {
                    if ( $super_line =~ />$chim_cluster/ ) { $chim_flag = 1; last; }
                }
            }

            # This cluster already recorded as part of a supercluster
            if ( $chim_flag == 1 ) {
                $line_no    = $.;
                $super_line = <CHIM_FILE>;
                chomp($super_line);
                $new_line = $super_line;

                for ( $n = 0 ; $n < $k ; $n++ ) {
                    if ( $super_line !~ /@$clustermatch[$n]/ ) {
                        $new_line .= @$clustermatch[$n];
                        $new_line .= ' ';
                        $super_line = $new_line;
                    }
                }
                close(CHIM_FILE);

                open( CHIM_FILE, "superclusters" );
                if ( !stat("tmpchim") ) { system("touch tmpchim"); }
                open( TMP_CHIM, ">tmpchim" );
                select(TMP_CHIM);
                while (<CHIM_FILE>) {
                    if   ( $. == $line_no + 1 ) { print TMP_CHIM "$new_line\n"; }
                    else                        { print TMP_CHIM $_; }
                }
                close(TMP_CHIM);
                close(CHIM_FILE);
                select(STDOUT);
                system("mv -f tmpchim superclusters");
            }

            if ( $chim_flag == 0 ) {
                $super_line = '';
                $new_line   = '';

                for ( $n = 0 ; $n < $k ; $n++ ) {
                    if ( $super_line !~ /@$clustermatch[$n]/ ) {
                        $new_line .= @$clustermatch[$n];
                        $new_line .= ' ';
                    }
                    $super_line = $new_line;
                }
                close(CHIM_FILE);

                system("echo '>@$clustermatch[$low_cluster]\n$super_line' >> superclusters");
                print LOG "***** @$clustermatch[$low_cluster] appears to be a supercluster *****\n";
            }
        }
    }
}

#################################################################################################
#################################################################################################
# Put zeros in front of the index to make it four digits long.
# Can add an extra line to enable clusters > 100000
sub index {
    my ($i) = @_;
    my ($index);
    if    ( $i < 10 )     { $index = "0000$i"; }
    elsif ( $i < 100 )    { $index = "000$i"; }
    elsif ( $i < 1000 )   { $index = "00$i"; }
    elsif ( $i < 10000 )  { $index = "0$i"; }
    elsif ( $i < 100000 ) { $index = "$i"; }
    else                  { die "index out of range\n"; }
}

# takes a file and copies it and its "out" file to the error directory
#sub problem {
#  my ($badfile) = @_;
#my @badarray = split('/', $badfile);
#my $badroot = pop( @badarray );
#my $badout = "$seq_files_done/$badroot.out";
#system "cp $badout $rootdir/$CLUSTERNAME.error";
#print LOG "Copying problem file $query_accession to error directory\n";
#}

# This is redundant but could be added to seperate different
# types of problem reads
#sub problemchim {
#    my ($badfile) = @_;
#    my @badarray = split('/', $badfile);
#    my $badroot = pop( @badarray );
#    my $badout = "$rootdir/$badroot.out";
#    system "mv $badfile $badout $rootdir/$CLUSTERNAME.chim";
#    #print LOG "Moving problem file $query_accession to Chim directory\n";
#}

###############################################################################
###############################################################################

sub newfile {
    my ( $file, $newstring ) = @_;
    my $original = "";
    system("mv $seq_files/$file $seq_files_done/$file.old");
    open( OLDFILE, "$seq_files_done/$file.old" );
    open( FILE,    ">$seq_files_done/$file" );
    while (<OLDFILE>) {
        print FILE unless /^>/;
        if (/^>/) {
            chop( $original = $_ );
            print FILE "$original $newstring\n";
        }
    }

    close(OLDFILE);
    close(FILE);
}

##################################################################################
##################################################################################

sub count_ns {
    my $n_max       = shift;
    my $sfile_start = shift;
    my $qfile_start = shift;
    my $s_add       = shift;
    my $q_add       = shift;
    my $accession   = shift;

    my ( $n_seq, $n_num, $line, $n );
    my (@seq_array);

    #print"$n_max $sfile_start $qfile_start $s_add $q_add\n";
    if ( -f "$seq_files/$accession" ) {
        open( SUBJECT_FILE, "<$seq_files/$accession" );
    } else {
        open( SUBJECT_FILE, "<$seq_files_done/$accession" );
    }

    $n_seq = '';
    $n_num = 0;
    while ( $line = <SUBJECT_FILE> ) {
        chop $line;
        if ( substr( $line, 0, 1 ) ne ">" ) {
            $n_seq .= $line;
        }
    }
    close(SUBJECT_FILE);

    @seq_array = split( '', $n_seq );
    for ( $n = 0 ; $n < $n_max ; $n++ ) {

        if ( $seq_array[ $sfile_start + ( $n * $s_add ) ] eq "N" ) {
            $n_num++;
        }
    }

    if ( $n_num / $n_max > 0.05 ) {
        $sl_flag = 0;

        #print LOG "$accession Overlap N's high ($n_num N's)\n";
    } else {
        if ( -f "$seq_files/$query_accession" ) {
            open( QUERY_FILE, "<$seq_files/$query_accession" );
        } else {
            open( QUERY_FILE, "<$seq_files_done/$query_accession" );
        }
        $n_seq = '';
        $n_num = 0;
        while ( $line = <QUERY_FILE> ) {
            chop $line;
            if ( substr( $line, 0, 1 ) ne ">" ) {
                $n_seq .= $line;
            }
        }
        @seq_array = split( '', $n_seq );
        for ( $n = 0 ; $n < $n_max ; $n++ ) {

            if ( $seq_array[ $qfile_start + ( $n * $q_add ) ] eq "N" ) {
                $n_num++;
            }
        }
        if ( $n_num / $n_max > 0.05 ) {
            $sl_flag = 0;

            #print LOG "$query_accession Overlap N's high ($n_num N's)\n";
        }
        close(QUERY_FILE);
    }
}

#############################################################################
#############################################################################

=head1 NAME

CLOBB.pl - A program for clustering sequences on the basis of BLAST similarity

=head1 SYNOPSIS

CLOBB.pl <cluster_id>

=head1 DESCRIPTION

This program takes a set of DNA sequences and clusters them into groups which
putatively derive from the same gene. In order to operate, the user must have
NCBI's blastall executable in their path. The output is a blastable fasta file named
<cluster_id>EST which contains a list of the sequences and a cluster identifier 
<cluster_id>EST00001 -> <cluster_id>EST99999. A number of other files and
directories are also created.

=head1 HOW IT WORKS

The program works by performing a BLAST for each sequence against the growing
cluster database. The BLAST output is then examined for High Scoring Pairs (HSPs) from the BLAST output
which demonstrate near identical regions of sequence similarity  
(>=95% sequence identity over >30 bases - these figures are modifiable in the   
source code to increase/decrease stringency) between the traget sequence and the
growing cluster database. These matches are classified as type I matches.

The next stage of the process goes through the list of type I matches to 
identify any problems associated with the match.  This is achieved by parsing 
the beginning and end positions of the query and subject sequences from the 
blast output.  If these positions overlap beyond the HSP by more than 30 bases  
(i.e. the HSP does not extend through the full overlap of the sequences), a further check is performed to ensure that this is not due to the presence of poor quality sequence (determined by the number
of bases assigned `N' in the overlap regions). Type I matches which do not have high quality overlaps of more than 30 bases beyond the HSP are designated as type II matches.  Other type I matches which possess high quality overlaps of greater than 30 bases which are not part of a HSP are designated as type III matches.

The next stage of cluster assignment then involves checking through the lists of type II and type III matches to ensure that no conflicts arise.  Given a cluster in which some members are type II matches, if there are other members of the same cluster which have been designated as type III matches, then this indicates that the query sequence matches some but not all members of a cluster and is therefore assigned a new cluster number.  The inclusion of this feature in the algorithm prevents the rapid expansion of chimeric clusters and can result in a splitting of related sequences into many different (related) clusters.  However, when such events occur, the program catalogues the clusters involved, identifying them as `similar to' the type III match for subsequent post-process analysis (typically performed by manual curation).

Another complication occurs when two or more type II matches arise from different clusters.  Firstly the blast output is reanalysed to determine whether the HSPs of the matches occur in overlapping regions.  If they do not, the query effectively links the clusters and they are merged into the cluster with the lowest index - a separate note is recorded to indicate that such a merge operation has occurred.  If they do overlap, this may indicate that they are either alternatively spliced variants of one gene or closely related members of a gene family, and the query sequence is assigned the cluster number of the type II match with which it had the highest blast score, providing that said cluster did not contain a type III match, and an annotation added to indicate that these clusters
may be members of a `Supercluster'.  Once the query sequence has been assigned to a cluster, it is added to the growing cluster database which is then reformatted to allow the next search.



=head2 Created files/Directories

=over 4

=item * supercluster

This file contains a list of clusters, members of which were found to be similar
to each other, however, other members of the cluster where found to be bad matches

=item * merge

This file contains a list of clusters which have been merged to form a single
cluster. This occurs when a sequence matches well with two clusters which do not overlap.

=item * master

A file containing a comma delimited list of the sequences.

=item * sequences_done

A directory listing all sequences which have been processed by CLOBB. Two
sequences are created for each original, the original sequence is renamed as
<sequence>.old, the new sequence file has the assigned cluster id annotated to
its header.

=item * OUT

A directory containing a logfile detailing the CLOBB process for each sequence.

=back

=head1 AUTHOR

John Parkinson (john.parkinson@ed.ac.uk)

=head1 REFERENCE

Parkinson J, Guiliano DB, Blaxter M.
Making sense of EST sequences by CLOBBing them.
BMC Bioinformatics. 2002 Oct 25;3(1):31.

=head1 COPYRIGHT

This program is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

=head1 SEE ALSO

perl(1).

=cut

