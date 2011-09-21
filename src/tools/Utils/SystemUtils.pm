#! /usr/bin/perl -w
package Utils::SystemUtils;
use strict;
use warnings;

##################### Find a specific program ######################################
sub find_program() {
    my @PATH     = split( ":", "$ENV{'PATH'}" );
    my $prog     = $_[0];
    my $pathflag = 0;
    my $path;
    my $finalpath;
    foreach $path (@PATH) {
        if ( -f "$path/$prog" ) {
            $pathflag  = 1;
            $finalpath = $path;
            last;
        }
    }
    if ( $pathflag == 0 ) {
        print "\nCan't find the $prog utility on your system\n";
        print "Please ensure it is installed and in your path\n\n";
        exit();
    } else {
        return "$finalpath/$prog";
    }
}
######################## End find_program() #########################################

1;