#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Utils::FastaUtils;
use Utils::SystemUtils;
use Utils::FileUtils;
use Utils::LogUtils;
use Cwd;
use Cwd qw(realpath);
use threads;
use Thread::Queue;

=head1 NAME

blast.pl - Blast nucleotide/protein sequence against available databases). Types of blasts: blastx, blastp and blastn.

=head1 SYNOPSIS

-project <txt> - project directory / project name [input]

-seq <txt> - sequence file (fasta format). Contains a list of nucleotide sequences [input]

-blast_type <txt> type of blast to perform (currently supports blastx, blastp and blastn) [input]

-blast_db <txt> database to blast against [input]

-top_nr_hits <txt> - only select the top x hits [input]

-e_value <txt> - e cut off value [input]

-word_size <txt> - word size, default if -1 (blastn 11, megablast 28, all others 3) [Integer]
    
-score_matrix <txt> - score_matrix: BLOSUM62, BLOSUM45, BLOSUM80, PAM30, and PAM70 [input]

-additional_blast_options <txt> - a string with additional blast options e.g. "-I T" GI numbers are included in hit definition list this option needs to be on if annotation is done with Blast2Go [input]

-blast_format <txt> - blast report output format (currently supports the blast default pairwise output (bls), asn (asn), xml (xml) or tabular with comment lines (tab)) [output]

-blast_report <txt> - blast reports  (combined in 1 file) [ouput]

-failed_blast <txt> - list of sequences that failed the blasting [ouput]

-nr_workers <txt> - number of worker threads to carry out the blasting [input]

-help - Get more detailed help

=head1 RUNNING

./blast.pl -project /home/user/project -seq /home/user/contigs.fasta -blast_type blastx -blast_db nr -top_nr_hits 10 -e_value 0.00001 -word_size 3 -score_matrix -blast_format xml -blast_report /home/user/reports.txt -failed_blast /home/user/failed.txt -nr_workers 5
 
=head1 OUTPUT

Directory structure:
    project_dir
        blast
            seq
            reports
            out 
                            
=head1 TODO
 

1) Check if e.g blastx is blasted against a protein database. Need to validate input first. In Galaxy validation can probably be done before tool is called.

2) Add more selection options e.g 
    1) low pass filter
    2) use composition-based
    3) number of processors to use (maybe not an option in Galaxy)

4) Also alow for blast default options:
    1) number of database sequence to show alignments = 250
    2) number of database sequences to show one-line descriptions for = 500

=cut

########################## main() ##################################################
{
	my $project_dir;
	my $base_dir;
	my %dirs;
	my $log_file;
	my (
		$fasta_sequences, $blast_type,  $blast_db,   $blast_report,
		$top_nr_hits,     $e_value,     $word_size,  $score_matrix,
		$add_blast_options, $blast_format,    $failed_hits, $nr_workers, $help
	);
	my $blast_db_location;    # get from environmental variable
	my $blast_exec;
   $add_blast_options = "";

	GetOptions(
		"project=s"      => \$project_dir,
		"seq=s"          => \$fasta_sequences,
		"blast_type=s"   => \$blast_type,
		"blast_db=s"     => \$blast_db,
		"top_nr_hits=s"  => \$top_nr_hits,
		"e_value=s"      => \$e_value,
		"word_size=s"    => \$word_size,
		"score_matrix=s" => \$score_matrix,
		"add_blast_options=s" => \$add_blast_options,
		"nr_workers=s"   => \$nr_workers,
		"blast_format=s" => \$blast_format,
		"blast_report=s" => \$blast_report,
		"failed_blast=s" => \$failed_hits,
		"help"           => \$help,
	);

	if (
		($help)
		|| !(
			   ($project_dir)
			&& ($fasta_sequences)
			&& ($blast_type)
			&& ($blast_db)
			&& ($top_nr_hits)
			&& ($e_value)
			&& ($word_size)
			&& ($score_matrix)
			&& ($blast_format)
			&& ($nr_workers)
			&& ($blast_report)
			&& ($failed_hits)
		)
	  )
	{
		print
"Blast nucleotide/protein sequence against available databases. Types of blasts: blastx, blastp and blastn. Returns blast report.\n";
		print
"Usage :\n\t do_blast.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print
"\t -seq <txt> - sequence file (fasta format). Contains a list of nucleotide sequences [input]\n";
		print
"\t -blast_type <txt> - type of blast to perform (currently supports blastx, blast and blastn) [input]\n";
		print "\t -blast_db <txt> - database to blast against [input]\n";
		print
		  "\t -top_nr_hits <txt> - only select top number of hits [input]\n";
		print "\t -e_value <txt> - e cut off value [input]\n";
		print
"\t -word_size, default if -1 (blastn 11, megablast 28, all others 3)\n";
		print
"\t -score_matrix <txt> - score_matrix <txt> - score_matrix: BLOSUM62, BLOSUM45, BLOSUM80, PAM30, and PAM70\n";
		print
"\t -additional_blast_options <txt> - a string with additional blast options e.g. \"-I T\" GI numbers are included in hit definition list this option needs to be on if annotation is done with Blast2Go [input]
\n";
		print
"\t -blast_format <txt> - blast report output format (currently supports the blast default pairwise output (bls), asn (asn), xml (xml) or tabular with comments (tab)) [output]\n";
		print
"\t -nr_workers <txt> - number of worker threads to carry out the blasting [input]\n";
		print
"\t -blast_report <txt> - blast reports  (combined in 1 file) [ouput]\n";
		print
"\t -failed_blast <txt> - list of sequences that failed the blasting [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$project_dir     = realpath($project_dir);
	$fasta_sequences = realpath($fasta_sequences);
	$blast_report    = realpath($blast_report);
	$failed_hits     = realpath($failed_hits);

	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);

	$base_dir = "$project_dir/blast_$timestamp";
	%dirs     = (
		"seq"     => "$base_dir/seq",
		"reports" => "$base_dir/reports",
		"out"     => "$base_dir/out"
	);

	$blast_db_location = $ENV{"BLASTDB"};

	my $file_queue = Thread::Queue->new();

	# Create blast directories
	&Utils::FileUtils::create_project_directories( $project_dir, $base_dir,
		%dirs );

	# Create log file
	$log_file = $base_dir . "/$timestamp.log";
	open( my $log_fd, ">$log_file" );

	# Write parameter information to logfile
	print "Base directory: " . $base_dir . "\n";
	&Utils::LogUtils::log_info( $log_fd, "Base directory: " . $base_dir );
	print "Fasta sequence: " . $fasta_sequences . "\n";
	&Utils::LogUtils::log_info( $log_fd,
		"Fasta sequence: " . $fasta_sequences );
	print "Blast type: " . $blast_type . "\n";
	&Utils::LogUtils::log_info( $log_fd, "Blast type: " . $blast_type );
	print "Blast database: " . $blast_db . "\n";
	&Utils::LogUtils::log_info( $log_fd, "Blast database: " . $blast_db );
	print "Top number of hits: " . $top_nr_hits . "\n";
	&Utils::LogUtils::log_info( $log_fd,
		"Top number of hits: " . $top_nr_hits );
	print "E-cutoff value: " . $e_value . "\n";
	&Utils::LogUtils::log_info( $log_fd, "E-cutoff value: " . $e_value );
	print "Word size: " . $word_size . "\n";
	&Utils::LogUtils::log_info( $log_fd, "Word size: " . $word_size );
	print "Score matrix: " . $score_matrix . "\n";
	&Utils::LogUtils::log_info( $log_fd, "Score matrix: " . $score_matrix );
	print "Blast output format: " . $blast_format . "\n";
	&Utils::LogUtils::log_info( $log_fd,
		"Blast output format: " . $blast_format );
	print "Number of working threads: " . $nr_workers . "\n";
	&Utils::LogUtils::log_info( $log_fd,
		"Number of working threads: " . $nr_workers );
	print "Blast report: " . $blast_report . "\n";
	&Utils::LogUtils::log_info( $log_fd, "Blast report: " . $blast_report );
	print "Failed hits: " . $failed_hits . "\n";
	&Utils::LogUtils::log_info( $log_fd, "Failed hits: " . $failed_hits );

	# Check if correct blast format is specified
	unless ( $blast_format =~ /^bls$|^asn$|^xml$|^tab$/ ) {
		print
"Output format not supported, currently blast default output (bls), asn (asn), xml (xml) or tabular with comments (tab) are supported.\n";
		&Utils::LogUtils::log_error( $log_fd,
"Output format not supported, currently blast default output (bls), asn (asn), xml (xml) or tabular with comments (tab) are supported."
		);
		exit();
	}

	&Utils::FastaUtils::split_fasta( $fasta_sequences, $dirs{"seq"} );

	$blast_exec = &Utils::SystemUtils::find_program("blastall");

	# Check if blast database directory exists
	if ( !-d "$blast_db_location" ) {
		print("Can't find your BLAST database directory: $blast_db_location\n");
		&Utils::LogUtils::log_error( $log_fd,
			"Can't find your BLAST database directory $blast_db_location" );
		close $log_fd;
		exit();
	}

	# Get the databases to blast against
	$blast_db =~ s/\s+/ /g;
	$blast_db =~ s/,/ /g
	  ; # This replacement is only added for galaxy, galaxy seperates a multiselect with a ","
	my @blast_db = split( ' ', $blast_db );
	$blast_db = "\"";
	for my $db (@blast_db) {
		$blast_db = $blast_db . $blast_db_location . '/' . $db;
		$blast_db = $blast_db . ' ';
	}
	chop $blast_db;    # remove trailing space
	$blast_db = $blast_db . "\"";

	# Check for protein databases
	if ( ( $blast_type eq "blastp" ) || ( $blast_type eq "blastx" ) ) {

		my @prot_dbs = glob("$blast_db_location/*.pin");
		if ( @prot_dbs < 1 ) {
			print "No protein databases can be found\n";
			&Utils::LogUtils::log_error( $log_fd,
				"No protein databases can be found" );
			close $log_fd;
			exit();
		}
		else {
			print "Found " . @prot_dbs . " protein databases\n";
			&Utils::LogUtils::log_info( $log_fd,
				"Found " . @prot_dbs . " protein databases" );
		}

	}

	# Check for nucleotide databases
	if ( ( $blast_type eq "blastn" ) ) {
		my @nuc_dbs = glob("$blast_db_location/*.nin");
		if ( @nuc_dbs < 1 ) {
			print "No nucleotide databases can be found\n";
			&Utils::LogUtils::log_error( $log_fd,
				"No nucleotide databases can be found" );
			close $log_fd;
			exit();
		}
		else {
			print "Found " . @nuc_dbs . " nucleotide databases\n";
			&Utils::LogUtils::log_info( $log_fd,
				"Found " . @nuc_dbs . " nucleotide databases" );
		}
	}

	chdir( $dirs{"seq"} );
	my @files = glob("*.fsa");

	for my $file (@files) {
#		print "Add $file to blast queue...\n";
		&Utils::LogUtils::log_info( $log_fd, "Add $file to blast queue..." );
		$file_queue->enqueue($file);
	}

	$file_queue->enqueue(undef);

	my $blast_success_file = $dirs{"out"} . "/success.$blast_format";
	my $blast_failed_file  = $dirs{"out"} . "/failed.txt";
	system "touch $blast_failed_file";

	# Start blasting
	for my $i ( 1 .. $nr_workers ) {
		if ( $word_size == -1 ) {
			$word_size = 0;
		}
		threads->create(
			\&blast,            $blast_exec,   $blast_type,
			$blast_db_location, $blast_db,     $dirs{"seq"},
			$dirs{"reports"},   $top_nr_hits,  $e_value,
			$word_size,         $score_matrix, $blast_format, $add_blast_options,
			$blast_failed_file, $file_queue
		);
	}

	for my $worker ( threads->list() ) {
		$worker->join();
	}

	# Create blast reports
	system(
		"cat " . $dirs{"reports"} . "/*.$blast_format > $blast_success_file" );
	system("cp $blast_success_file $blast_report");
	system("cp $blast_failed_file $failed_hits");
	close $log_fd;
}

#############################  End main() ###########################################

##############################  blast() #############################################
sub blast() {
	my (
		$blast_exec,   $blast_type,   $blast_db_location,
		$blast_db,     $seq_dir,      $report_dir,
		$top_nr_hits,  $e_value,      $word_size,
		$score_matrix, $blast_format, $add_blast_options, $blast_failed_file,
		$file_queue
	) = @_;

	my $format = "";

	if ( $blast_format eq "bls" ) {
		$format = "-m 0";

	}
	elsif ( $blast_format eq "asn" ) {
		$format = "-m 10";

	}
	elsif ( $blast_format eq "xml" ) {
		$format = "-m 7";

	}
	elsif ( $blast_format eq "tab" ) {
		$format = "-m 9";

	}
	else {
		print
"Output format not supported, currently blast default output (pairwise), asn (asn), xml (xml) or tabular with comments (tab) are supported\n";
		exit();
	}

	my $count = 0;
	while ( my $file = $file_queue->dequeue() ) {
#		print "Worker " . threads->tid() . ": blasting $file ...\n";
		my $flag   = 0;
		my $report = $file;
		$report =~ s/\.seq|\.fsa|\.fasta/\.$blast_format/;
		while ( $flag < 1 ) {

			system( "$blast_exec -p $blast_type -d $blast_db -i " 
				  . "\"" . $seq_dir
				  . "/$file\" $format -o "
				  . "\"" . $report_dir
				  . "/$report\" -v $top_nr_hits -e $e_value -b $top_nr_hits -W $word_size -M $score_matrix $add_blast_options 2> /dev/null"
			);

			if ( -s $report_dir . "/$report" > 500 ) {
#				print "$blast_type report generated for $file\n";
				last;
			}
			$flag++;
		}
		if ( $flag == 1 ) {
#			print "$blast_type cannot be done on $file\n";
			system
"echo \"$blast_type cannot be done on $file\" > $blast_failed_file";
		}
		system "mv \"$file\" \"$file.worker." . threads->tid() . "\"";
		$count++;
	}

	$file_queue->enqueue(undef);
#	print "Worker " . threads->tid() . " blasted $count files\n";
}
############################## end blast() ##########################################

