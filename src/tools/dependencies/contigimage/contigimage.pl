#!/usr/bin/perl

# @(#) $RCSfile: contigimage,v $ $Revision: 1.4 $ $Date: 2002/10/29 19:10:58 $ $Name:  $
# 
# Center for Computational Genomics and Bioinformatics,
# Academic Health Center, University of Minnesota
# Copyright (c) 2002. The Regents of the University of Minnesota
# 
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, see
# 
#     http://www.gnu.org/copyleft/gpl.html
# 
# or write to the Free Software Foundation, Inc., 59 Temple Place, 
# Suite 330, Boston, MA  02111-1307  USA
# 

use GD;

local ($new_dir, $af, $coid);
local ($scaleh, $consensush, $csh, $colorhash);
local ($width, $height, $theight, $minbase, $maxbase);
local (@COARRAY, @BQARRAY, @RDID, @RDARRAY, %SENSE, %OFFSET);

$af = $ARGV[0];

# Put the created images into a directory, this could be a cmd line option
$new_dir = $ARGV[1];
if (! -d $new_dir) {
    unlink $new_dir;			# in case it exists, but is a file
    if (!mkdir $new_dir) {		# setting mode doesn't work here
	die "unable to make dir $new_dir $!\n";
    }
    chmod 0775, $new_dir;
}

&do_ace($af, $new_dir);

#exit;

# Global variables used 
#
#		COARRAY BQARRAY are same length, almost same index
#		There is no BQ entry if the CO entry is "*", meaning gap
#     @COARRAY			contains consensus sequence 1 char per index
#     @BQARRAY			contains quality array 1 number per index.
#
#		RDID RDARRAY same length, indexed by read number
#     @RDID			names of the reads
#     @RDARRAY			array of array references, each element of
#				RDARRAY points to an array which contains 
#				the read sequence 1 char per index
#
#		SENSE and OFFSET hash tables indexed by rdid (read name)
#     %SENSE			sense of each read
#     %OFFSET			offset of read wrt consensus, may be negative

#
# do_ace gathers all the contig info into the arrays and hashes
# described above.
#
sub do_ace {
    my $af = shift;
    my $new_dir = shift;
    open(ACE, "<$af") or die "couldn't open ace file: $!\n";

    # need to be set every time through
    $minbase =  100000;
    $maxbase = -100000;

    $/ = "";				# break on blank lines, not newlines
    while (<ACE>) {
	chomp;
	if (/^CO\s+(\w+)\s+(\d+)\s+(\d+)/) {
	    &ace_array;			# process previous contig
	    my $coline = $_;
	    $coid = $1;			# contig ID
	    my $colen = $2;		# length of consensus seq
	    my $coreadn = $3;		# number of reads in contig
	    $coline =~ s/^CO.*\n?//;	# remove the first line (^CO)
	    $coline =~ s/\s+//go;		# remove whitespace
	    @COARRAY = split (//, $coline);
	}
	if (/^BQ/) {			# base quality info
	    my $bqline = $_;		# save line for editing
	    $bqline =~ s/^BQ.*\n?//;	# remove the first line (^BQ)
	    $bqline =~ s/^\s+//;	# remove leading whitespace
	    @BQARRAY = split (/\s+/, $bqline);
	} 
	if (/^RD\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
	    my $rdid = $1;
	    my $rdlen = $2;
	    push (@RDID, $rdid);

	    my $end = $OFFSET{$rdid} + $rdlen;
	    if ($end > $maxbase) {$maxbase = $end;}

	    my $rdline = $_;
	    $rdline =~ s/^RD.*\n?//;		# remove the first line (^RD)
	    $rdline =~ s/\s+//go;		# remove whitespace
	    my @rdarray = split (//, $rdline);
	    push(@RDARRAY, [@rdarray]);		# array of arrays
	}
	if (/^AF/) {
	    my @afarray = split "\n";		# may be many lines of AF
	    my $afele;
	    foreach $afele (@afarray) {
		my ($afmarker, $afid, $afsense, $afoffset);
		($afmarker, $afid, $afsense, $afoffset) = split ' ', $afele;
		last if ($afmarker ne "AF");
		# save sense and offset according to afid (read name)
		$SENSE{$afid} = $afsense;
		$OFFSET{$afid} =  $afoffset;
		if($afoffset < $minbase) { $minbase = $afoffset; }
	    }
	}
	# all other lines are ignored
    } 

    close(ACE);
    &ace_array;					# process last contig
}

sub ace_array {
	if (not defined($coid)) {return;}		# first time, nothing to do
	if (0) {					# could have debug flag
	    print "#### BEGIN ace_array ####\n";
	    print "CO info: $coid $colen $coreadn \n";
	    print "COARRAY size = $#COARRAY\n";		#, COARRAY = @COARRAY\n\n";
	    print "BQARRAY size = $#BQARRAY\n";		#, BQARRAY = @BQARRAY\n\n";
	    print "RDID size = $#RDID, RDID = @RDID\n\n";
	    print "RDARRAY size = $#RDARRAY\n";		#, RDARRAY = @RDARRAY\n\n";
	    foreach $rdarray (@RDARRAY) {
		print "rdarray = @$rdarray\n\n";
	    }
	    for ($i = 0; $i <= $#RDARRAY; $i++) {
		$rar  = $RDARRAY[$i];
		@rdarray = @$rar;
		# for now generate gel on left
		for ($j = 0; $j <= $#rdarray; $j++) {
		    $l = $rdarray[$j];
		    print "$l ";
		}
		print "\n\n";
	    }
	
	    @sensearray = sort keys %SENSE;
	    foreach $sense (@sensearray) {
		print "sense{$sense} = $SENSE{$sense}\n";
	    }
	    @offsetarray = sort keys %OFFSET;
	    foreach $offset (@offsetarray) {
		print "offset{$offset} = $OFFSET{$offset}\n";
	    }
	    print "minbase, maxbase = $minbase, $maxbase\n";
	    print "#### END ace_array ####\n";
	}
	
	&gd_contig ($new_dir);			# create the contig image
	# reinitialize all arrays and hashes to null
	@COARRAY = ();
	@BQARRAY = ();
	@RDARRAY = ();
	@RDID = ();
	%SENSE = ();
	%OFFSET = ();
	
	# reset for the next contig
	$minbase =  100000;
	$maxbase = -100000;
}

#
# gd_contig creates a new contig image 
#
sub gd_contig {
    $new_dir = shift;

    # These are needed now to figure the height of the image
    $scaleh = 30;					# scale height
    $consensush = 30;					# consensus height
    $csh = $scaleh + $consensush;			# bottom of consensus

    $width = $maxbase - $minbase;		# minbase can be negative
    $height = $csh + 10 * ($#RDID + 1) + 4;	# scale + consensus + each read
    $theight = @RDID + 1;			# height of thumbnail image
    my $im = new GD::Image($width + 100, $height);
    my $tim = new GD::Image($width/4, $theight);

    # allocate some colors
    # [A = green, C = cyan, G = orange, T=red]
    # [N = gray, X = darkgray, * = darkgray]

    my $ip;
    foreach $ip ($im, $tim) {
	# first color allocated will be background
	$black = $ip->colorAllocate(0, 0, 0);
	$white = $ip->colorAllocate(255, 255, 255);
	$blue = $ip->colorAllocate(0, 0, 255);		# not used

	$Acolor = $green = $ip->colorAllocate(0, 255, 0);
	$Ccolor = $cyan = $ip->colorAllocate(0, 255, 255);
	$Gcolor = $orange = $ip->colorAllocate(255, 165, 0);
	$Tcolor = $red = $ip->colorAllocate(255, 0, 0);
	$Ncolor = $gray = $ip->colorAllocate(165, 165, 165);
	$Ucolor = $Xcolor = $dgray = $ip->colorAllocate(127, 127, 127);
    }

    # initialize letter to color lookup hash table
    $colorhash{"A"} = $Acolor;
    $colorhash{"a"} = $Acolor;
    $colorhash{"C"} = $Ccolor;
    $colorhash{"c"} = $Ccolor;
    $colorhash{"G"} = $Gcolor;
    $colorhash{"g"} = $Gcolor;
    $colorhash{"T"} = $Tcolor;
    $colorhash{"t"} = $Tcolor;
    $colorhash{"N"} = $Ncolor;
    $colorhash{"n"} = $Ncolor;
    $colorhash{"X"} = $Xcolor;
    $colorhash{"x"} = $Xcolor;
    $colorhash{"*"} = $Ucolor;

    gd_consensus($im);				# draw scale and consensus
    gd_reads($im, $tim);			# draw all reads

    # Convert the image to PNG and print it to file
    # perhaps add a verbose option, lots of output here
    #print "Creating png file $new_dir/$coid.png\n";
    open(PNG, ">$new_dir/$coid.png") || die "open $new_dir/$coid.png $!\n";
    # make sure we are writing to a binary stream
    binmode PNG;
    print PNG $im->png;
    close PNG;


	# The only modification made to the original contigimage.pl 
	# was to remove the printing of the thumbnail pictures
	# in the next few lines
	

    # only write thumbnail if the contig is deep enough
    #if ($theight > 4) {
        #print "Creating png file $new_dir/T$coid.png\n";
	#open(PNG, ">$new_dir/T$coid.png") || 
	#    die "open $new_dir/T$coid.png $!\n";
	# make sure we are writing to a binary stream
	#binmode PNG;
	#print PNG $tim->png;
	#close PNG;
    #}


}

#
# creates consensus and scale at top of image
#
sub gd_consensus {
    $im = shift;

    # Add scale to top of image and vertical bars every 100 bp
    for ($i = 0; $i < $maxbase; $i+=100) {
	# allow 100 pixels for the name field
	$x = 100 - $minbase + $i - 2;			# why 2?
	#print "i = $i, x = $x\n";
	my $fonthere = GD::Font->Small;
	$im->stringUp($fonthere, $x + 3, $scaleh - 5, $i, $white);
	$im->line($x + 2, $scaleh - 9, $x + 2, $height, $gray);
    }
    $im->line(100 - $minbase, $scaleh - 2, 100 + $maxbase, $scaleh - 2, $white);

    # Add consensus name to upper left corner
    my $fonthere = GD::Font->MediumBold;
    $im->string($fonthere, 6, $scaleh + 4, $coid, $white);
    $im->string($fonthere, 6, $scaleh + 14, "CONSENSUS", $white);

    # Add consensus image underneath scale
    $q = 0;						# BQ array index
    for ($i = 0; $i <= $#COARRAY; $i++) {
	$l = $COARRAY[$i];				# get letter
	$colorindex = $colorhash{$l};			# get color
	# allow 100 pixels for the name field
	$x = 100 + $i - $minbase;
	$y = $csh - 5;
	if ($l eq "*") {				# gap has no BQ
	    $bq = 0;
	} else {
	    $bq = $BQARRAY[$q];				# get the bq value
	}
	# calculate the height of the bar based on the quality in BQARRAY
	# minimum height is 5, maximum 30
	if ($bq < 5) {$y = $csh - 5;}
	if ($bq > 30) {$y = $csh - $consensush;}
	if (($bq >= 5) && ($bq <= 30) ) {
	    $y = $csh - $bq;
	}
	#print "i = $i, q = $q, CO = $COARRAY[$i], BQ = $bq, y = $y\n";
	$im->line($x, $y, $x, $csh - 2, $colorindex);	# paint stripe
	$q++ unless ($l eq "*");
    }
}

#
# Go through RDARRAY arrays to generate the read images
#
sub gd_reads {
    my $im = shift;
    my $tim = shift;
    my ($i, $j);

    # Add reads to image
    for ($i = 0; $i <= $#RDARRAY; $i++) {
	my ($rar, $rdid, $y, $fg, $offset, @rdarray);
	# put the name on the left, allowing 100 pixels
	$rar  = $RDARRAY[$i];		# this needs two statements
	@rdarray = @$rar;
	$rdid = $RDID[$i];			# name (ID) of THIS read
	$y = $csh +2 + $i * 10;
	# use white fg for normal sequences, cyan fg for reverse complemented
	$fg = $white;
	if ($SENSE{$rdid} eq "C") {$fg = $cyan;}
	my $fonthere = GD::Font->Small;
	$im->string($fonthere, 1, $y, $rdid, $fg);

	#
	# draw read image itself
	# lower case indicates lower quality, use half height
	# each time thru paints one vertical stripe, one pixel wide
	# color is base, height is quality
	#
	$offset = 100 + $OFFSET{$rdid};		# offset for THIS read
	$toffset = $OFFSET{$rdid};		# thumb offset for THIS read
	for ($j = 0; $j <= $#rdarray; $j++) {
	    my ($l, $y, $y1, $y2, $x, $colorindex);
	    $l = $rdarray[$j];				# get letter
	    $y = $csh + 5;				# lower case half high
	    if ($l =~ /[ACGT]/) {$y = $csh + 1;}	# full height
	    $y1 = $y + $i * 10;
	    $y2 = $csh + 9 + $i * 10;
	    $x = $j - $minbase + $offset - 1;
	    $colorindex = $colorhash{$l};		# get color
	    $im->line($x, $y1, $x, $y2, $colorindex);	# paint stripe
	    $tim->setPixel(($x -100)/4, $i, $colorindex);	# thumb pixel
	}
    }
}
