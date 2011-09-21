#! /usr/bin/perl -w
package Utils::QualScoreUtils;
use strict;
use warnings;

########### Split file with a list of quality scores into individual files ##########
################### e.g. split_qual( "all.qual", "qual" ) ###########################
sub split_qual() {
    my $file   = shift(@_);
    my $dir    = shift(@_);
    my $first  = 1;
    my $header = "";
    my $data   = "";
    open( QUALFILE, "$file" ) || die "Can't open $file\n";
    while ( my $line = <QUALFILE> ) {

        if ( $line =~ /^>(\S+)/ ) {
            if ( $first == 1 ) {

                # Fist sequence score
                $header = $1;
                open( OUTFILE, ">$dir/$header.qual" );
                $first = 0;
            } else {

                # Write scores of previous sequence to file
                print OUTFILE ">$header\n";
                my @qual_scores = split( ' ', $data );
                my $n = 0;
                for ( my $i = 0 ; $i < @qual_scores ; $i++ ) {
                    print OUTFILE "$qual_scores[$i] ";
                    $n++;
                    if ( $n % 30 == 0 ) { print OUTFILE "\n"; }
                }
                if ( $n % 30 != 0 ) {
                    print OUTFILE "\n";
                }
                close OUTFILE;

                # Start reading new sequence data
                $header = $1;
                $data   = "";
                open( OUTFILE, ">$dir/$header.qual" );
            }
        }
        if ( $line !~ /^>/ ) {
            chomp($line);
            $data = $data . " " . $line;
        }
    }

    # Write last record
    print OUTFILE ">$header\n";
    my @qual_scores = split( ' ', $data );
    my $n = 0;
    for ( my $i = 0 ; $i < @qual_scores ; $i++ ) {
        print OUTFILE "$qual_scores[$i] ";
        $n++;
        if ( $n % 30 == 0 ) { print OUTFILE "\n"; }
    }
    if ( $n % 30 != 0 ) {
        print OUTFILE "\n";
    }
    close OUTFILE;
    close QUALFILE;
}
############################ End split_qual() #######################################

1;