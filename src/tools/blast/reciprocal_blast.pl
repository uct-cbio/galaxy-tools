#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(realpath);
use Utils::SystemUtils;
use Bio::SearchIO;

# do a reciprocal blast. A sequences are blasted against database B. B sequences is blasted against database A.
{
	my ( $seq_type, $blast_format, $sequence_file_x, $sequence_file_y,
		$prog_work_dir, $delete_program_work, $help );
	my ( $blast_report_x, $blast_report_y, $hits_report_x, $hits_report_y );
	my ( %hits_result_x, %hits_result_y );
	my $blast_type;
	my $rbh_file;
	my @dir;
	my @files;
	
	$delete_program_work = '';

	GetOptions(
		"seq_type=s"   => \$seq_type,
		"xf=s"           => \$sequence_file_x,
		"yf=s"           => \$sequence_file_y,
		"w=s"            => \$prog_work_dir,
		"rbh=s"            => \$rbh_file,
		"d" => \$delete_program_work,
		"help=s"         => \$help,
	);

	if (
		($help)
		|| !(
			   ($seq_type)
			&& ($sequence_file_x)
			&& ($sequence_file_y)
			&& ($prog_work_dir)
			&& ($rbh_file)
		)
	  )
	{
		print "Usage :\n\t reciprocal_blast.pl <list of arguments> \n\n";
		print "\t -seq_type <txt> - sequence type (n) nucleotide or (p) protein [input]\n";
		print "\t -xf <txt> - fasta file x with nucleutide/protein sequences depending on blast type\n";
		print "\t -yf <txt> - fasta file y with nucleutide/protein sequences depending on blast type\n";
		print "\t -rbh <txt> - reciprocal best hits file\n";
		print "\t -w <txt> - program work directory\n";
		print "\t -d <txt> - add flag to delete program work directory when done\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$sequence_file_x = realpath($sequence_file_x);
	$sequence_file_y = realpath($sequence_file_y);
	$rbh_file = realpath($rbh_file);
	
	$blast_format = "xml";

	-e $sequence_file_x or die "$sequence_file_x does not exits \n";
	-e $sequence_file_y
	  or die "$sequence_file_y does not exits \n";

	$prog_work_dir = realpath($prog_work_dir);

	if ( -e "$prog_work_dir" ) {
		system("rm -rf $prog_work_dir");
		system("mkdir $prog_work_dir");
	}
	else {
		system("mkdir $prog_work_dir");
	}
	# ---------------
	# Copy sequence file to program work directory.
	# Keep blast database files seperate from Galaxy data directory
	# Not time efficient, can consider to do create database in Galaxy data directory later 
	# ---------------
	system("cp $sequence_file_x $prog_work_dir/seqx.fsa");
	system("cp $sequence_file_y $prog_work_dir/seqy.fsa");
	$sequence_file_x = $prog_work_dir. "/seqx.fsa";
	$sequence_file_y = $prog_work_dir. "/seqy.fsa";
	# ---------------

	my $report_x = $sequence_file_x;
	$report_x =~ s/.*\///;
	my $report_y = $sequence_file_y;
	$report_y =~ s/.*\///;

	$blast_report_x = $prog_work_dir . "/" . $report_x . ".bls";
	$blast_report_y = $prog_work_dir . "/" . $report_y . ".bls";
	$hits_report_x  = $prog_work_dir . "/" . $report_x . "_best_hits.csv";
	$hits_report_y  = $prog_work_dir . "/" . $report_y . "_best_hits.csv";

	print "Do reciprocal blast ...\n";

	my $formatdb_type;
	if ($seq_type eq "n" ) {
		$formatdb_type = "F";
		$blast_type = "blastn";
	}
	elsif ($seq_type eq "p" ) {
		$formatdb_type = "T";
		$blast_type = "blastp";
	}
	else {
		print "Specify (n) nucleotide or (p) protein\n";
		exit();
	}

	# first create blast databases
	my $formatdb_exec = &Utils::SystemUtils::find_program("formatdb");

	system "$formatdb_exec -i $sequence_file_x -p $formatdb_type 2> /dev/null";
	system "$formatdb_exec -i $sequence_file_y -p $formatdb_type 2> /dev/null";

	my $blastall_exec = &Utils::SystemUtils::find_program("blastall");

	print $blastall_exec . "\n";

	if ( $blast_format eq "bls" ) {
		system
"$blastall_exec -p $blast_type -i $sequence_file_x  -d $sequence_file_y -m 0 -e 0.01 -o $blast_report_x 2> /dev/null"
		  ;    # blast with low pass filtering
		system
"$blastall_exec -p $blast_type -i $sequence_file_y  -d $sequence_file_x -m 0 -e 0.01 -o $blast_report_y 2> /dev/null"
		  ;    # blast with low pass filtering
		       # get hit results from blast report
		&blast_parsing( $blast_report_x, 'blast', $hits_report_x,
			\%hits_result_x );
		&blast_parsing( $blast_report_y, 'blast', $hits_report_y,
			\%hits_result_y );

	}

	if ( $blast_format eq "xml" ) {
		system
"$blastall_exec -p $blast_type -i $sequence_file_x  -d $sequence_file_y -m 7 -e 0.01 -o $blast_report_x 2> /dev/null"
		  ;    # blast with low pass filtering
		system
"$blastall_exec -p $blast_type -i $sequence_file_y  -d $sequence_file_x -m 7 -e 0.01 -o $blast_report_y 2> /dev/null"
		  ;    # blast with low pass filtering
		       # get hit results from blast report
		&blast_parsing( $blast_report_x, 'blastxml', $hits_report_x,
			\%hits_result_x );
		&blast_parsing( $blast_report_y, 'blastxml', $hits_report_y,
			\%hits_result_y );
	}

	# create reciprocal best hit file
	&create_rbh_file( \%hits_result_x, \%hits_result_y, $rbh_file );
	
	if ($delete_program_work){
		system "rm -rf 	$prog_work_dir";
	}

	print "Done.\n";
}
#############################################################################################
# Parsing the BLAST output file.
#############################################################################################
#############################################################################################
sub blast_parsing {

	my ( $file_in, $format, $hits_report_file, $hits_result ) = @_;

	my $blast_io = new Bio::SearchIO(
		-format => $format,
		-file   => $file_in
	);

	open( BLAST, "<$file_in" );
	open( HITS,  ">$hits_report_file" );
	while ( my $result = $blast_io->next_result ) {
		while ( my $hit = $result->next_hit ) {
			my $query        = $result->query_name;
			my $best_hit     = $hit->name();
			my $significance = $hit->significance();
			print HITS "$query\t$best_hit\t$significance\n";
			$$hits_result{$query} = $best_hit;
			last;    # only read first hit
		}
	}
	close(BLAST);
	close(HITS);
}

sub create_rbh_file {
	my ( $result_x, $result_y, $file ) = @_;
	my ( $query, $subject );

	open( RBH, ">$file" );

	while ( ( $query, $subject ) = each %$result_x ) {
		if ( exists $$result_y{$subject} ) {
			if ( $$result_y{$subject} eq $query ) {
				print RBH "$query\t$subject\n";
			}
		}
	}
	close(RBH);
}
