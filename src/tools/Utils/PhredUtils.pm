#! /usr/bin/perl -w
package Utils::PhredUtils;
use strict;
use warnings;

########### Split file with a list of phred info into individual clusters ###########
############ e.g. split_phred( "all_phred.txt", "phred") ############################
sub split_phred() {
    my $file = shift(@_);
    my $dir  = shift(@_);
    my $flag = 1;
    open( PHREDFILE, "$file" ) || die "Can't open $file\n";
    while ( my $line = <PHREDFILE> ) {
        if ( $line =~ /BEGIN_SEQUENCE\s(.+)/ ) {
            open( OUTFILE, ">$dir/$1.phd.1" );
            $flag = 0;
        }
        if ( $line =~ /END_SEQUENCE/ ) {
            print OUTFILE $line;
            close OUTFILE;
            $flag = 1;
        }
        if ( $flag == 0 ) {
            print OUTFILE $line;
        }
    }
    close PHREDFILE;
}
############################ End split_phred() ######################################

1;