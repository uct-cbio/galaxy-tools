#!/usr/bin/perl

=head1 NAME

url_blast.pl - Blast nucleotide/protein sequence against available databases at NCBI using the URL blast structure. Types of blasts: blastx, blastp, blastn, tblastn and tblastx.

=head1 SYNOPSIS

-seq <txt> - sequence file (fasta format). Contains a list of nucleotide sequences [input]

-blast_type <txt> type of blast to perform (currently supports blastx, blastp, blastn, tblastn and tblastx) [input]

-blast_db <txt> database to blast against (e.g. nt, nr) [input]

-blast_format <txt> - blast output format (bls(pairwise), xml, asn and html) [input]

-e_value <txt> - e cut off value [input]

-word_size <txt> - default if -1 (blastn 11, megablast 28, all others 3) (LOOKS LIKE -1 DOES NOT WORK. SPECIFY MANUALY!) [input]

-score_matrix <txt> - score_matrix <txt> - score_matrix: BLOSUM62, BLOSUM45, BLOSUM80, PAM30, and PAM70 [input]. If program gives \"Search expire\" "
		  . "error the default values for the scoring matrix probably needs to be set. Add this to the -other_advanced option. "
		  . "Default values can be found here http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node91.html#sub:Protien-Scoring-Matrices. [input]
		  
-hitlist_size <txt> - Sets the upper limit for number of alignments and descriptions available in the formatting and reformatting step [input]

-lc_filter <txt> - filter query sequence (DUST with blastn, SEG with others) (T/F) [input]

-cbased_stats <txt> - composition-based score adjustments for blastp or tblastnt (T = composition based, F = no composition based)[input]

-entrez <txt> - entrez query (e.g. \"Zea mays[Organism] OR Oryza sativa[Organism] OR Arabidopsis thaliana[Organism]\")[input]

-other_advanced <txt> - add other advanced options to url query e.g. \"-G 9 -E 1\"[input]

-number_seq_query <txt> - number of sequences per blast query [input]

-poll_results <txt> - poll for results every x seconds [input]

-blast_report <txt> - output blast report. A concatenated file of all single single reports. [output]

-successful_blasts <txt> - sequences that have been successfully blasted [output]

-failed_blasts <txt> - sequences that failed the blast. Rerun url_blast with these sequences as input. It is possible to retrieve more results. [output]

-prog_work <txt> - program work directory. Can be removed by user after results have been retrieved with the -d flag. [output]

-d - remove program work directory with this flag [input]

-help - Get more detailed help

=head1 RUNNING

~/projects/uct_cbio_galaxy_tools/src/tools/blast/url_blast.pl -seq small_set_2.fsa -blast_type blastp -blast_db nr -blast_format xml -e_value 20000 -word_size 2 -score_matrix "PAM30" -other_advanced "-G 9 -E 1" -hitlist_size 10 -lc_filter F -cbased_stats T  -entrez "Mycobacterium tuberculosis H37Rv[Organism]" -number_seq_query 10  -poll_results 10  -blast_report small_set_blast_2.xml -successful_blasts successful_blast.fsa -failed_blasts failed_blast.fsa -prog_work out 

=head1 OUTPUT

Directory structure:
    project_dir
        url_blast
            seq
            reports
            out 
                            
=head1 TODO

=cut

use URI::Escape;
use LWP::UserAgent;
use HTTP::Request::Common qw(POST);
$ua = LWP::UserAgent->new;
use Cwd qw(realpath);
use Getopt::Long;
use Utils::FastaUtils;
use Utils::FileUtils;

{

	my (
		$blast_type, $blast_db,  $blast_format,
		$e_value,    $word_size, $score_matrix,
		$hitlist_size,
		$lc_filter, $cbased_stats,
		$entrez_query, $other_advanced, $fsa_file, $number_seq_query,
		$poll_results, $blast_report, $successful_blasts_fsa,
		$failed_blasts_fsa,
		$prog_work_dir, $delete_program_work

	);
	my ( $format_type, $successful_blasts, $failed_blasts, $blast_error );
	my %dirs;
	my $base_dir;

	$entrez_query   = "";
	$other_advanced = "";

	$hitlist_size = 100;    #default

	GetOptions(
		"seq=s"               => \$fsa_file,
		"blast_type=s"        => \$blast_type,
		"blast_db=s"          => \$blast_db,
		"blast_format=s"      => \$blast_format,
		"e_value=s"           => \$e_value,
		"word_size=s"         => \$word_size,
		"score_matrix=s"      => \$score_matrix,
		"hitlist_size=s"      => \$hitlist_size,
		"lc_filter=s"         => \$lc_filter,
		"cbased_stats=s"      => \$cbased_stats,
		"entrez=s"            => \$entrez_query,
		"other_advanced=s"    => \$other_advanced,
		"number_seq_query=i"  => \$number_seq_query,
		"poll_results=i"      => \$poll_results,
		"blast_report=s"      => \$blast_report,
		"successful_blasts=s" => \$successful_blasts_fsa,
		"failed_blasts=s"     => \$failed_blasts_fsa,
		"prog_work=s"         => \$prog_work_dir,
		"d"                   => \$delete_program_work,
		"help"                => \$help,
	);

	if (
		($help)
		|| !(
			   ($fsa_file)
			&& ($blast_type)
			&& ($blast_db)
			&& ($blast_format)
			&& ($e_value)
			&& ($word_size)
			&& ($score_matrix)
			&& ($lc_filter)
			&& ($cbased_stats)
			&& ($number_seq_query)
			&& ($poll_results)
			&& ($hitlist_size)
			&& ($blast_report)
			&& ($successful_blasts_fsa)
			&& ($failed_blasts_fsa)
			&& ($prog_work_dir)
		)
	  )
	{
		print
"Blast nucleotide/protein sequence against available databases at NCBI using the URL blast. Types of blasts: blastx, blastp, blastn, tblastn and tblastx.\n";
		print "Usage :\n\t url_blast.pl <list of arguments> \n\n";
		print "\t -seq <txt> - fasta sequence file\n";
		print
"\t -blast_type <txt> - blast type (e.g. blastn, blastx, blastp)[input]\n";
		print "\t -blast_db <txt> - blast database (e.g. nt, nr)[input]\n";
		print
"\t -blast_format <txt> - blast output format (bls(pairwise), xml, asn and html) [input]\n";
		print "\t -e_value <txt> - -e_value <txt> - e cut off value[input]\n";
		print
"\t -word_size <txt> - default if -1 (blastn 11, megablast 28, all others 3) (LOOKS LIKE -1 DOES NOT WORK SPECIFY MANUALY) [input]\n";
		print
"\t -score_matrix <txt> - score_matrix <txt> - score_matrix: BLOSUM62, BLOSUM45, BLOSUM80, PAM30, and PAM70 [input]. If program gives \"Search expire\" "
		  . "error the default values for the scoring matrix probably needs to be set. Add this to the -other_advanced option. "
		  . "Default values can be found here http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node91.html#sub:Protien-Scoring-Matrices.\n";
		print
"\t -hitlist_size <txt> - Sets the upper limit for number of alignments and descriptions available in the formatting and reformatting step [input]\n";
		print
"\t -lc_filter <txt> - filter query sequence (DUST with blastn, SEG with others) (T/F) [input]\n";
		print
"\t -cbased_stats <txt> - composition-based score adjustments for blastp or tblastnt (T = composition based, F = no composition based)[input]\n";
		print
"\t -entrez <txt> - entrez query (e.g. \"Zea mays[Organism] OR Oryza sativa[Organism] OR Arabidopsis thaliana[Organism]\")[input]\n";
		print
"\t -other_advanced <txt> - add other advanced options to url query e.g. \"-G 9 -E 1\"[input]\n";
		print
		  "\t -number_seq_query <txt> - number of sequences per blast query\n";
		print "\t -poll_results <txt> - poll for results every x seconds\n";
		print
"\t -blast_report <txt> - output blast report. A concatenated file of all single single reports.\n";
		print
"\t -successful_blasts <txt> - sequences that have been successfully blasted\n";
		print
"\t -failed_blasts <txt> - sequences that failed the blast. Rerun url_blast with these sequences as input. It is possible to retrieve more results.\n";
		print
"\t -prog_work <txt> - program work direcotry. Can be removed by user after results have been retrieved with the -d flag.\n";
		print "\t -d - remove program work directory\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$prog_work_dir         = realpath($prog_work_dir);
	$fsa_file              = realpath($fsa_file);
	$blast_report          = realpath($blast_report);
	$failed_blasts_fsa     = realpath("$failed_blasts_fsa");
	$successful_blasts_fsa = realpath("$successful_blasts_fsa");

	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);

	$base_dir = "$prog_work_dir/url_blast_$timestamp";
	%dirs     = (
		"seq"     => "$base_dir/seq",
		"reports" => "$base_dir/reports",
		"out"     => "$base_dir/out"
	);

	# Create blast directories
	&Utils::FileUtils::create_project_directories( $prog_work_dir, $base_dir,
		%dirs );

	$successful_blasts =
	  $dirs{"out"} . "/successful_blasts.txt";    # a list of successful blasts
	$failed_blasts =
	  $dirs{"out"} . "/failed_blasts.txt";        # a list of failed blasts
	  
	if( $blast_type eq "blastn"){
		$score_matrix = "";
	}

	if ( $blast_format eq "bls" ) {
		$format_type = "Text";

	}
	elsif ( $blast_format eq "asn" ) {
		$format_type = "ASN.1";

	}
	elsif ( $blast_format eq "xml" ) {
		$format_type = "XML";

	}
	elsif ( $blast_format eq "html" ) {
		$format_type = "HTML";

	}
	else {
		print
"Output format not supported, currently blast pairwise output (bls), xml, asn and html are supported\n";
		exit();
	}

	if ( $lc_filter eq "T" ) {
		$lc_filter = "L";
	}
	elsif ( $lc_filter eq "F" ) {
		$lc_filter = ""; # if filtering is not specified nothing needs to be set
	}
	else {
		print "Wrong lc_filter option, specify T/F\n";
		exit();
	}

	if ( $cbased_stats eq "T" ) {
		$cbased_stats = "yes";
	}
	elsif ( $cbased_stats eq "F" ) {
		$cbased_stats = "no";
	}
	else {
		print "Wrong $cbased_stats option, specify T/F\n";
		exit();

	}

	if ( $word_size == -1 ) {
		$word_size = 0;
	}

	&Utils::FastaUtils::split_n_fasta( $fsa_file, $dirs{"seq"},
		$number_seq_query );

	chdir( $dirs{"seq"} );
	system "touch $successful_blasts";
	system "touch $failed_blasts";
	system "touch $successful_blasts_fsa";
	system "touch $failed_blasts_fsa";

	@dir = glob("*.fsa");

	# Need to replace characters generated by Galaxy
	my %mapped_chars = (
		'>'  => '__gt__',
		'<'  => '__lt__',
		'\'' => '__sq__',
		'"'  => '__dq__',
		'['  => '__ob__',
		']'  => '__cb__',
		'{'  => '__oc__',
		'}'  => '__cc__'
	);

#	print "before replacement\n$entrez_query\n$other_advanced\n";

	foreach my $key ( keys %mapped_chars ) {
		$entrez_query =~ s/$mapped_chars{$key}/$key/g;
		$other_advanced =~ s/$mapped_chars{$key}/$key/g;
	}
	
#	print "after replacement\n$entrez_query\n$other_advanced\n";

	$entrez_query   = uri_escape($entrez_query);
	$other_advanced = uri_escape($other_advanced);

	if (1) {
		foreach my $file (@dir) {

#			print "Query sequence = $file \n";
			$blast_error = 0;

			# read and encode the queries
			open( QUERY, $file );
			$encoded_query = "";
			while (<QUERY>) {
				$encoded_query = $encoded_query . uri_escape($_);
			}
			close QUERY;

			# build the request
			$args =
			    "CMD=Put&PROGRAM=$blast_type&DATABASE=$blast_db&QUERY="
			  . $encoded_query
			  . "&ENTREZ_QUERY="
			  . $entrez_query
			  . "&EXPECT="
			  . $e_value
			  . "&WORD_SIZE="
			  . $word_size
			  . "&MATRIX_NAME="
			  . $score_matrix
			  . "&HITLIST_SIZE="
			  . $hitlist_size
			  . "&FILTER="
			  . $lc_filter
			  . "&COMPOSITION_BASED_STATISTICS="
			  . $cbased_stats
			  . "&OTHER_ADVANCED="
			  . $other_advanced;

#			print $args . "\n";

			$req = new HTTP::Request POST =>
			  'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi';
			$req->content_type('application/x-www-form-urlencoded');
			$req->content($args);

			# get the response
			$response = $ua->request($req);

			# parse out the request id
			$response->content =~ /^    RID = (.*$)/m;
			$rid = $1;
#			print "RID: " . $rid . "\n";

			# parse out the estimated time to completion
			$response->content =~ /^    RTOE = (.*$)/m;
			$rtoe = $1;

			# wait for search to complete
			sleep $rtoe;

			# poll for results
			while (1) {
				sleep $poll_results;

				$req = new HTTP::Request GET =>
"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";
				$response = $ua->request($req);

				if ( $response->content =~ /\s+Status=WAITING/m ) {
#					print "$file: Searching...\n";
					next;
				}

				if ( $response->content =~ /\s+Status=FAILED/m ) {
#					print
#"$file: Search $rid failed; please report to blast-help\@ncbi.nlm.nih.gov.\n";
					$blast_error = 1;
					last;

				}

				if ( $response->content =~ /\s+Status=UNKNOWN/m ) {
#					print "$file: Search $rid expired.\n";
					$blast_error = 1;
					last;
				}

				if ( $response->content =~ /\s+Status=READY/m ) {
					if ( $response->content =~ /\s+ThereAreHits=yes/m ) {
#						print "$file: Search complete, retrieving results...\n";
						last;
					}
					else {
#						print "$file: No hits found.\n";
						$blast_error = 1;
						last;
					}
				}

#				print "$file, Something unexpected happened.\n";
				$blast_error = 1;
				last;
			}    # end poll loop

			# retrieve and display results
			unless ($blast_error) {
				$req = new HTTP::Request GET =>
"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=$format_type&RID=$rid";
				$response = $ua->request($req);
				my $report_file = $file;
				$report_file =~ s/\.fsa/\.$blast_format/;
				open( $report_fd, ">" . $dirs{"reports"} . "/" . $report_file );
				print $report_fd $response->content;
				close $report_fd;

				open( my $successful_blasts_fd, ">>$successful_blasts" );
				print $successful_blasts_fd $file . "\n";
				close $successful_blasts_fd;
				system("cat $file >> $successful_blasts_fsa");
			}
			else {
				open( my $failed_blasts_fd, ">>$failed_blasts" );
				print $failed_blasts_fd $file . "\n";
				close $failed_blasts_fd;
				system("cat $file >> $failed_blasts_fsa");
			}

		}
	}


	my @files = glob($dirs{"reports"} . "/*" . "\.$blast_format");
	my $nr_files = scalar(@files);
	
	if($nr_files > 0){
		system("cat " . $dirs{"reports"} . "/*" . "\.$blast_format > $blast_report" );
	}else{
		system("touch  $blast_report");
	}
	
	
	if ($delete_program_work) {
		system "rm -rf 	$base_dir";
	}

}
