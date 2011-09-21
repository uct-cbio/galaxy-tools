#! /usr/bin/perl -w
package Utils::FastaUtils;
use strict;
use warnings;

########### Split file with a list of fasta sequences into individual files #########
################### e.g. split_fasta( "all.n", "fasta" ) ############################
sub split_fasta() {
	my $file   = shift(@_);
	my $dir    = shift(@_);
	my $flag   = 0;
	my $in_seq = "";
	my @heading;
	my $header;
	my $filename;
	open( FSAFILE, "$file" ) || die "Can't open $file\n";

	while ( my $line = <FSAFILE> ) {
		if ( $flag == 1 && $line !~ /^>/ ) {
			chop $line;
			$in_seq .= $line;
		}
		if ( $flag == 1 && $line =~ /^>/ ) {    # print to file
			$flag = 0;
			open( OUTFILE, ">$dir/$filename.fsa" )
			  || die "Can't open $dir/$filename.fsa\n";
			print OUTFILE "$header\n";
			my $n = 0;
			while ( $in_seq =~ /(.)/g ) {
				print OUTFILE "$1";
				$n++;
				if ( ( $n % 60 ) == 0 ) { print OUTFILE "\n"; }
			}
			if ( ( $n % 60 ) != 0 ) { print OUTFILE "\n"; }
			$in_seq = "";
			close OUTFILE;
		}
		if ( $flag == 0 && substr( $line, 0, 1 ) eq ">" ) {  # an unusual regex.
			chop $line;
			$flag    = 1;
			$header  = $line;
			@heading = split( ' ', $line );
			$heading[0] =~ s/>//;
			$filename = $heading[0];
			next;
		}
	}

	# Print last record to file
	open( OUTFILE, ">>$dir/$filename.fsa" )
	  || die "Can't open $filename.fsa\n";
	print OUTFILE "$header\n";
	my $n = 0;
	while ( $in_seq =~ /(.)/g ) {
		print OUTFILE "$1";
		$n++;
		if ( ( $n % 60 ) == 0 ) { print OUTFILE "\n"; }
	}
	if ( ( $n % 60 ) != 0 ) { print OUTFILE "\n"; }
	close OUTFILE;
}
############################ End split_fasta() ######################################

############ Slit a fasta file. split file with several sequences ###################
################### into files containing only x sequences ##########################
# new file = x_y.fsa, where x = file/split number and y number of sequences in file #
######### for last file number of sequences can be less then specified ##############
################### e.g. split_n_fasta( "all.n", "fasta" ) ##########################
sub split_n_fasta() {
	my $file   = shift(@_);
	my $dir    = shift(@_);
	my $x      = shift(@_);
	my $flag   = 0;
	my $in_seq = "";
	my @heading;
	my $header;
	my $filename;
	my $count   = 0;
	my $file_nr = 0;
	my $content = "";
	open( FSAFILE, "$file" ) || die "Can't open $file\n";

	while ( my $line = <FSAFILE> ) {
		if ( $flag == 1 && $line !~ /^>/ ) {
			chop $line;
			$in_seq .= $line;
		}
		if ( $flag == 1 && $line =~ /^>/ ) {    # print to file
			$flag = 0;

			$content = $content . "$header\n";
			my $n = 0;
			while ( $in_seq =~ /(.)/g ) {
				$content = $content . "$1";
				$n++;
				if ( ( $n % 60 ) == 0 ) { $content = $content . "\n"; }
			}
			if ( ( $n % 60 ) != 0 ) { $content = $content . "\n"; }
			$in_seq = "";

			$count++;

			if ( $count eq $x ) {
				$file_nr++;
				# file number + number of sequences in file
				open( OUTFILE, ">$dir/$file_nr\_$count\.fsa" )
				  || die "Can't open $dir/$file_nr\_$count\.fsa";
				print OUTFILE "$content";
				close OUTFILE;
				$count   = 0;
				$content = "";
			}
		}
		if ( $flag == 0 && substr( $line, 0, 1 ) eq ">" ) {  # an unusual regex.
			chop $line;
			$flag    = 1;
			$header  = $line;
			@heading = split( ' ', $line );
			$heading[0] =~ s/>//;
			$filename = $heading[0];
			next;
		}
	}

	# Print last record to file
	$count++;
	if ( $count >= 1) {
		
		$content = $content . "$header\n";
		my $n = 0;
		while ( $in_seq =~ /(.)/g ) {
			$content = $content . "$1";
			$n++;
			if ( ( $n % 60 ) == 0 ) { $content = $content . "\n"; }
		}
		if ( ( $n % 60 ) != 0 ) { $content = $content . "\n"; }

		$file_nr++;
		# file number + number of sequences in file
		open( OUTFILE, ">$dir/$file_nr\_$count\.fsa" )
		  || die "Can't open $dir/$file_nr\_$count\.fsa\n";
		print OUTFILE "$content";
		close OUTFILE;
		$count = 0;

	}

}
############################ End split_fasta() ######################################

1;
