#! /usr/bin/perl -w
# Some parts are similar to the trimming part in trace2dbest.pl but most are rewritten and new
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd qw(realpath);
use Utils::FileUtils;
use Utils::SystemUtils;
use Utils::PhredUtils;
use Utils::FastaUtils;
use Utils::QualScoreUtils;

=head1 NAME

trim.pl - Trim vector, ecoli, poly A trail sequences from a fasta file with sequences sequences. Quality scores are also trimmed to correspond with the trimming in the input sequences. 

=head1 SYNOPSIS

-project <txt> - project directory / project name [input]

-seq <txt> - fasta sequence file [input]

-qual <txt> - quality scores file [input]

-phred <txt> - phred file [input]

-vector <txt> - vector file [input]

-ecoli <txt> - ecoli file [input]

-5prime_adapter <txt> - 5 prime adapter sequence [input]

-3prime_adapter <txt> - 3 prime adapter sequence [input]

-trimseq <txt> - trimmed fasta sequence file [output]

-trimqual <txt> - trimmed quality scores file [output]

-help - Get more detailed help

=head1 RUNNING

./trim.pl -project /home/user/project -seq /home/user/all.fsa -qual /home/user/all.qual -phred /home/user/all_phred.txt -vector /home/user/vector.seq -ecoli /home/user/ecoli.seq -5prime_adapter "AATTCGGCACGAGG" -3prime_adapter "AAAAAAAAAAAACTCGAG" -trimseq /home/user/all_trimmed.fsa -trimqual /home/user/all_trimmed.qual

=head1 OUTPUT

Directory structure:
    project_dir
        trim
            fasta
            qual
            phred
            fasta_trim
            qual_trim
            process
            subfiles
            fasta_fail
            qual_fail

=head1 TODO
    
    - Sequences are trimmed 1 by 1.  Can probably be faster if a concatenated file is trimmed.

=cut

{
########################## main() ##################################################

	my $project_dir;
	my $base_dir;
	my %dirs;
	my $log_file;
	my (
		$fasta_sequence_file,        $quality_score_file,
		$phred_file,                 $vector_file,
		$ecoli_file,                 $prime5_adapter, $prime3_adapter, $trimmed_fasta_sequence_file,
		$trimmed_quality_score_file, $help
	);

	GetOptions(
		"project=s"  => \$project_dir,
		"seq=s"      => \$fasta_sequence_file,
		"qual=s"     => \$quality_score_file,
		"phred=s"    => \$phred_file,
		"vector=s"   => \$vector_file,
		"ecoli=s"    => \$ecoli_file,
		"5prime_adapter=s"    => \$prime5_adapter,
		"3prime_adapter=s"    => \$prime3_adapter,
		"trimseq=s"  => \$trimmed_fasta_sequence_file,
		"trimqual=s" => \$trimmed_quality_score_file,
		"help"       => \$help,
	);

	if (
		($help)
		|| !(
			   ($project_dir)
			&& ($fasta_sequence_file)
			&& ($quality_score_file)
			&& ($phred_file)
			&& ($vector_file)
			&& ($ecoli_file)
			&& ($prime5_adapter)
			&& ($prime3_adapter)
			&& ($trimmed_fasta_sequence_file)
			&& ($trimmed_quality_score_file)
		)
	  )
	{
		print
" Trim file with fasta sequences (vector, ecoli and poly (A,T) trail). Corresponding quality files are also trimmed.  Returns trimmed sequences (fasta file) with corresponding trimmed quality scores (qual file).\n";
		print
"Usage :\n\t trim.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print "\t -seq <txt> - fasta sequence file [input]\n";
		print "\t -qual <txt> - quality scores file [input]\n";
		print "\t -phred <txt> - phred file [input]\n";
		print "\t -vector <txt> - vector file [input]\n";
		print "\t -ecoli <txt> - ecoli file [input]\n";
		print "\t -5prime_adapter <txt> - 5 prime adapter sequence [input]\n";
		print "\t -3prime_adapter <txt> - 3 prime adapter sequence[input]\n";
		print "\t -trimseq <txt> - trimmed fasta sequence file [output]\n";
		print "\t -trimqual <txt> - trimmed quality scores file [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$project_dir                 = realpath($project_dir);
	$fasta_sequence_file         = realpath($fasta_sequence_file);
	$quality_score_file          = realpath($quality_score_file);
	$phred_file                  = realpath($phred_file);
	$vector_file                 = realpath($vector_file);
	$ecoli_file                  = realpath($ecoli_file);
	$trimmed_fasta_sequence_file = realpath($trimmed_fasta_sequence_file);
	$trimmed_quality_score_file  = realpath($trimmed_quality_score_file);

	# my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	# my $timestamp = `date +%Y-%m-%d_%H_%M_%S`; # cross_match has problems if %N added. Too long directory path?
	# chomp($timestamp);
	# $base_dir = "$project_dir/trim_$timestamp";
	
	# For now give random number. (cross-match issues). Sort this out later!! 
	my $timestamp = int(rand(10000));
	$base_dir = "$project_dir/trim_$timestamp";

	%dirs     = (
		"fasta"      => "$base_dir/fasta",
		"qual"       => "$base_dir/qual",
		"phd"        => "$base_dir/phd",
		"fasta_trim" => "$base_dir/fasta_trim",
		"qual_trim"  => "$base_dir/qual_trim",
		"fasta_fail" => "$base_dir/fasta_fail",
		"qual_fail"  => "$base_dir/qual_fail",
		"out"        => "$base_dir/out"
	);

# Adapter sequences for for Xerophyta humilis project. Was originally hardcoded. Should not be in vector sequence file!
#	$prime5_adapter   = "AATTCGGCACGAGG";
#	$prime3_adapter   = "AAAAAAAAAAAACTCGAG";
	my $high_qual_cutoff = 150;

	# Create trim directories
	&Utils::FileUtils::create_project_directories( $project_dir, $base_dir,
		%dirs );

	# Create log file
	$log_file = $base_dir . "/$timestamp.log";
	open( my $log_fd, ">$log_file" );

	# Copy files to trim directory
	system "cp $fasta_sequence_file $base_dir/sequences.fsa";
	$fasta_sequence_file = "$base_dir/sequences.fsa";

	system "cp $quality_score_file $base_dir/sequences.qual";
	$quality_score_file = "$base_dir/sequences.qual";

	system "cp $phred_file $base_dir/sequences.phred";
	$phred_file = "$base_dir/sequences.phred";

	system "cp $vector_file $base_dir/vector.fsa";
	$vector_file = "$base_dir/vector.fsa";

	system "cp $ecoli_file $base_dir/ecoli.fsa";
	$ecoli_file = "$base_dir/ecoli.fsa";

	# Split into individual files
	&Utils::FastaUtils::split_fasta( $fasta_sequence_file, $dirs{"fasta"} );
	&Utils::QualScoreUtils::split_qual( $quality_score_file, $dirs{"qual"} );
	&Utils::PhredUtils::split_phred( $phred_file, $dirs{"phd"} );

	# Process and trim
	&process(
		$base_dir, $dirs{"fasta"},      $dirs{"qual"},      $dirs{"phd"},
		$dirs{"fasta_trim"}, $dirs{"qual_trim"}, $dirs{"fasta_fail"},
		$dirs{"qual_fail"}, $prime5_adapter,
		$prime3_adapter, $vector_file, $ecoli_file
	);

	# Create trimmed sequence and quality scores output
	system("cat " . $dirs{"fasta_trim"} ."/*" . " > " . $dirs{"out"} . "/trimmed_sequences.fsa");
	system("cp " . $dirs{"out"} . "/trimmed_sequences.fsa" . " " . $trimmed_fasta_sequence_file);

	system("cat " . $dirs{"qual_trim"} ."/*" . " > " . $dirs{"out"} . "/trimmed_sequences.qual");
	system("cp " . $dirs{"out"} . "/trimmed_sequences.qual" . " " . $trimmed_quality_score_file);
}
#############################  End main() ##########################################

####################################################################################
sub process() {
	my (
		$base_dir, $fasta_dir,      $qual_dir,       $phred_dir,     $fasta_trim_dir,
		$qual_trim_dir,  $fasta_fail_dir, $qual_fail_dir, $prime5_adapter,
		$prime3_adapter, $vector_file,    $ecoli_file
	) = @_;
	my ( $vminmatch, $vminscore, $eminmatch, $eminscore, $high_qual_cutoff,
		$poly_num );
	$vminmatch        = 10;
	$vminscore        = 20;
	$eminmatch        = 20;
	$eminscore        = 30;
	$high_qual_cutoff = 150;
	$poly_num         = 12;
	
	chdir "$base_dir";
	
	if ( -e "cmlogfile" ) {
		print "Removing old cross_match logfile\n";
		unlink "cmlogfile";
	}

	opendir( FASTA_DIR, $fasta_dir );
	while ( my $fasta_file = readdir FASTA_DIR ) {
		if ( $fasta_file =~ /(seq$|fasta$|fsa$)/ ) {

# TODO to give some option of naming scheme. For now only NERC e.g. Xh_RD_01A01. Do some checks on naming and be sure quality score files are available.
			my $qual_score_file = $fasta_file;
			$qual_score_file =~ s/(seq$|fasta$|fsa$)/qual/;
			my $phred_score_file = $fasta_file;
			$phred_score_file =~ s/(seq$|fasta$|fsa$)/phd.1/;
			open( FASTA, "$fasta_dir/$fasta_file" );
			open( QUAL,  "$qual_dir/$qual_score_file" );
			my @fasta_lines      = <FASTA>;
			my @qual_score_lines = <QUAL>;
			my $sequence         = "";
			my $quality_scores   = "";
			my $header           = ">" . $fasta_file . "\n";
			$header =~ s/(\.seq|\.fasta|\.fsa)//;

			for ( my $i = 1 ; $i < @fasta_lines ; $i++ ) {
				$sequence = $sequence . $fasta_lines[$i];
				chomp $sequence;
			}
			for ( my $i = 1 ; $i < @qual_score_lines ; $i++ ) {
				$quality_scores = $quality_scores . $qual_score_lines[$i];
				chomp $quality_scores;
			}
			# do adaptor trimming, otherwise vector trimming may take off a little bit of adapter meaning that this section wouldn't recognise it
			&flag_adapter( \$sequence, $prime5_adapter, $prime3_adapter )
			  ;    # 5 prime end
			&flag_adapter( \$sequence, $prime3_adapter, $prime5_adapter )
			  ;    # repeat for 3 prime end;
			  
			&flag_vector( \$sequence, "$fasta_dir/$fasta_file", $vector_file,
				$vminmatch, $vminscore );
			system
			  "cp $fasta_dir/$fasta_file.screen $fasta_dir/$fasta_file.ecoli"
			  ;    #ready for e-coli trimming
			  
			&flag_vector( \$sequence, "$fasta_dir/$fasta_file.ecoli",
				$ecoli_file, $eminmatch, $eminscore );

# Flag (X) everything from beginning to first occurence X, the same for the reverse
			if ( $sequence =~ /^([^X]{1,150})(X.*)$/ ) {
				my $flag_x_begin = $1;
				$sequence = $2;
				$flag_x_begin =~ s/./X/g;
				$sequence = $flag_x_begin . $sequence;
			}
			my $reverse_sequence = reverse $sequence;
			if ( $reverse_sequence =~ /^([^X]{1,150})(X.*)$/
			  )    # TODO 150 was set in trace2dbest, do it here as well
			{
				my $flag_x_begin = $1;
				$reverse_sequence = $2;
				$flag_x_begin =~ s/./X/g;
				$reverse_sequence = $flag_x_begin . $reverse_sequence;
			}
			$sequence = reverse $reverse_sequence;
			&flag_polyA( \$sequence, $high_qual_cutoff, $poly_num );
			&flag_polyT( \$sequence, $high_qual_cutoff, $poly_num );
			&flag_low_quality( \$sequence, "$phred_dir/$phred_score_file" )
			  ;    # sorts out unwanted regions not flag by others
			$sequence =~ s/H/X/ig;
			&flag_N( \$sequence )
			  ; # flag X's identified in the middle can be ecoli vector or something else

			# Trim sequence and quality scores by removing X flagged nucleotides
			# Write trimmed sequence and quality scores
			open( FASTA_TRIM, ">$fasta_trim_dir/$fasta_file" );
			open( QUAL_TRIM,  ">$qual_trim_dir/$qual_score_file" );
			my @seq  = split( '',  $sequence );
			my @qual = split( ' ', $quality_scores );
			my $seq_length  = @seq;
			my $qual_length = @qual;
			print FASTA_TRIM $header;
			print QUAL_TRIM $header;
			my $n = 0;

			for ( my $i = 0 ; $i < @seq ; $i++ ) {
				if ( $seq[$i] =~ /([^X])/ ) {
					print FASTA_TRIM "$seq[$i]";
					print QUAL_TRIM "$qual[$i] ";
					$n++;
					if ( ( $n % 60 ) == 0 ) {
						print FASTA_TRIM "\n";
					}
					if ( ( $n % 30 ) == 0 ) {
						print QUAL_TRIM "\n";
					}
				}
			}
			if ( ( $n % 30 ) != 0 ) {
				print QUAL_TRIM "\n";
			}
			if ( ( $n % 60 ) != 0 ) {
				print FASTA_TRIM "\n";
			}
			close FASTA_TRIM;
			close QUAL_TRIM;
			if ( $n < $high_qual_cutoff ) {
				system("mv $fasta_trim_dir/$fasta_file $fasta_fail_dir/");
				system("mv $qual_trim_dir/$qual_score_file $qual_fail_dir/");
			}
			close FASTA;
			close QUAL;
		}
	}
	close(FASTA);
}
############################ End process ######### ##################################

################ Flag(X) regions containing adapter sequences #######################
sub flag_adapter() {
	my ( $sequence, $prime5_adapter, $prime3_adapter ) = @_;

# 5 prime
#    if ( $$sequence =~ /^(.{0,200}$prime5_adapter)(.*)/i ) TODO give some threshold, for now strip everything
	if ( $$sequence =~
		/^(.*$prime5_adapter)(.*)/i )    # TODO for now strip everything
	{
		my $trim_adapter = $1;
		$$sequence = $2;
		$trim_adapter =~ s/./X/g;
		$$sequence = $trim_adapter . $$sequence;
	}

# 3 prime
#    if ( $$sequence =~ /^(.*)($prime3_adapter.{0,200})(.*)/i ) # TODO give some threshold, for now strip everything
	if ( $$sequence =~ /^(.*)($prime3_adapter.*)/i ) {
		$$sequence = $1;
		my $trim_adapter = $2;    # TODO still trim more than 200 nucleotides
		$trim_adapter =~ s/./X/g;
		$$sequence = $$sequence . $trim_adapter;
	}
}
############################ End flag_adapter() #####################################

################ Flag(X) regions containing any vector sequence #####################
sub flag_vector() {
	my ( $sequence, $fasta_file, $vector_file, $vminmatch, $vminscore ) = @_;

	# Find cross_match executable
	my $cross_match_exec = &Utils::SystemUtils::find_program("cross_match");

	# Run cross_match
#	system(
#"$cross_match_exec  $fasta_file $vector_file -minmatch $vminmatch -minscore $vminscore -screen 2>>cmlogfile 1> /dev/null"
#	); 

	system(
"$cross_match_exec  $fasta_file $vector_file -minmatch $vminmatch -minscore $vminscore -screen 1> /dev/null 2>&1"
	); # Galaxy has trouble when another process writes to a file. Rather dump everything to /dev/null

	# read screen file into sequence
	open( FASTA, "$fasta_file.screen" );
	my @flagged_fasta_lines = <FASTA>;
	my $flagged_sequence    = "";
	for ( my $i = 1 ; $i < @flagged_fasta_lines ; $i++ ) {
		$flagged_sequence = $flagged_sequence . $flagged_fasta_lines[$i];
		chomp $flagged_sequence;
	}
	my @flagged_seq = split( '', $flagged_sequence );
	my @seq         = split( '', $$sequence );
	my $fl          = @flagged_seq;
	my $sl          = @seq;
	for ( my $i = 0 ; $i < @flagged_seq ; $i++ ) {
		if ( $flagged_seq[$i] =~ /[X]/ ) {
			if ( $seq[$i] =~ /[^X]/ ) {
				$seq[$i] = $flagged_seq[$i];
			}
		}
	}
	$$sequence = "@seq";
	$$sequence =~ s/\s//g;
}
############################ End flag_vector() ######################################

################ Flag(X) regions containing poly A tail #############################
sub flag_polyA() {
	my ( $sequence, $high_qual_cutoff, $poly_num ) = @_;
	my $polya_cut = $high_qual_cutoff - $poly_num;
	if ( $$sequence =~ /(.{$polya_cut,}?)(a{$poly_num}.*)/i ) {
		my $pa_trim = $2;
		$$sequence = $1;
		$pa_trim =~ s/./X/g;
		$$sequence = $$sequence . $pa_trim;
	}
}
############################ End flag_polyA() #######################################

################ Flag(X) regions containing poly T tail #############################
sub flag_polyT() {
	my ( $sequence, $high_qual_cutoff, $poly_num ) = @_;
	my $polya_cut = $high_qual_cutoff - $poly_num;

	# a reverse poly A algorithm
	my $reverse_sequence = reverse $$sequence;
	if ( $reverse_sequence =~ /(.{$polya_cut,}?)(t{$poly_num}.*)/i ) {
		my $pt_trim = $2;
		$reverse_sequence = $1;
		$pt_trim          = reverse $pt_trim;
		$pt_trim =~ s/./X/g;
		$$sequence = reverse $reverse_sequence;
		$$sequence = $pt_trim . $$sequence;
	}
}
############################ End flag_polyT() #######################################

####################### Flag(X) low quality regions #################################
sub flag_low_quality() {
	my ( $sequence, $phred_score_file ) = @_;
	my $p_hq_start = 0;
	my $p_hq_end   = 0;
	open( PHRED, "$phred_score_file" );
	while ( my $line = <PHRED> ) {
		if ( $line =~ /TRIM:\s(-?\d+)\s(-?\d+)/ ) {
			$p_hq_start = $1;
			$p_hq_end   = $2;
			last;
		}
	}

	# now flag low quality values
	my $lqb_length = $p_hq_start - 1;
	my $hq_length  = $p_hq_end - $p_hq_start + 1;
	if ( $$sequence =~ /(^.{$lqb_length})(.{$hq_length})(.+)/ ) {
		my $lqb_seq = $1;
		my $q_seq   = $2;
		my $lqe_seq = $3;
		$lqb_seq =~ s/./X/g;
		$lqe_seq =~ s/./X/g;
		$$sequence = $lqb_seq . $q_seq . $lqe_seq;
	}
	close PHRED;
}
############################ End flag_low_quality() #################################

############################ Flag(N) unsure regions #################################
sub flag_N() {
	my ($sequence) = @_;
	my $trim_begin = "";
	my $trim_end   = "";
	my $trim_seq   = "";
	if ( $$sequence =~ /(^X*)([^X].*)/i ) {
		$trim_begin = $1;
		$trim_seq   = $2;
	}
	my $reverse_trim_seq = reverse $trim_seq;
	if ( $reverse_trim_seq =~ /(^X*)([^X].*)/i ) {
		$trim_end         = $1;
		$reverse_trim_seq = $2;
		$trim_end         = reverse $trim_end;
		$trim_seq         = reverse $reverse_trim_seq;
	}
	$trim_seq =~ s/X/N/ig;
	$$sequence = $trim_begin . $trim_seq . $trim_end;
}
################################## End flag_N() #####################################

