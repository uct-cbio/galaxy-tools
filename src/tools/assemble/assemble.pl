#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd qw(realpath);
use Utils::SystemUtils;
use Utils::PhredUtils;
use Utils::FastaUtils;
use Utils::QualScoreUtils;
use Utils::FileUtils;

=head1 NAME

assemble.pl - Assemble contigs using phrap and get ready for blasting.

=head1 SYNOPSIS

-project <txt> - project directory / project name [input]

-clone_seq <txt> -  clone sequences (fasta format) [input]

-clone_qual <txt> - clone quality scores (qual format) [input]

-cluster_seq <txt> - cluster sequences (fasta format) [input]

-singleton_seq <txt> - singleton sequences (fasta format) [input]

-clone_id <txt> - cluster id [input]

-cluster_id <txt> - cluster id [input]

-contig_seq <txt> - assembled/contig sequences (fasta format) [output]

-contig_qual <txt> - assembled/contig quality scores (qual format) [output]

help\t - Get more detailed help

=head1 RUNNING

./assemble.pl -project /home/user/project -clone_seq /home/user/clones.fsa -clone_qual  /home/user/clones.qual -cluster_seq /home/user/clusters.fsa -singleton_seq /home/user/singletons.fsa -clone_id Xh -cluster_id XHC -contig_seq /home/user/contig.fsa -contig_qual /home/user/contig.qual

=head1 OUTPUT

Directory structure:
    project_dir
        assemble
            clusters (contains clustered equences and singleton.fasta file)
            seq_qual (copy original sequence and quality files here)
            phrap
            contigs
            log
            out
 
=head1 TODO

1) Check header format.

2) Some variable values (DEFAULT_QUALSCORE and DEFAULT_FORCELEVEL) are still hardcoded.

3) Do logging.

=cut

my $DEFAULT_QUALSCORE =
  15;    # quality score used in case no phred quality files are present
my $DEFAULT_FORCELEVEL =
  10;    # force level used in phrap to reduce number of multiple contigs

########################## main() ##################################################
{

	my $project_dir;
	my $base_dir;
	my %dirs;
	my $log_file;

	my (
		$cluster_fasta_seq,  $singleton_fasta_seq, $clone_fasta_seq,
		$clone_qual_scores,  $cluster_id,          $contig_fasta_seq,
		$contig_qual_scores, $clone_id,            $help
	);
	my @assemble_dirs = (
		"assemble", "clusters", "fsa_qual", "phrap",
		"contigs",  "log",      "output"
	);
	my ( $phrap_exec, @pseudo_contig, $seq, $line );

	GetOptions(
		"project=s"       => \$project_dir,
		"clone_seq=s"     => \$clone_fasta_seq,
		"clone_qual=s"    => \$clone_qual_scores,
		"cluster_seq=s"   => \$cluster_fasta_seq,
		"singleton_seq=s" => \$singleton_fasta_seq,
		"clone_id=s"      => \$clone_id,
		"cluster_id=s"    => \$cluster_id,
		"contig_seq=s"    => \$contig_fasta_seq,
		"contig_qual=s"   => \$contig_qual_scores,
		"help"            => \$help,
	);

	if (
		($help)
		|| !(
			   ($project_dir)
			&& ($cluster_fasta_seq)
			&& ($singleton_fasta_seq)
			&& ($clone_fasta_seq)
			&& ($clone_qual_scores)
			&& ($contig_fasta_seq)
			&& ($contig_qual_scores)
			&& ($cluster_id)
			&& ($clone_id)
		)
	  )
	{
		print "Assemble contigs using phrap and get ready for blasting. .\n";
		print
"Usage :\n\t assemble.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print "\t -clone_seq <txt> -  clone sequences (fasta format) [input]\n";
		print
		  "\t -clone_qual <txt> - clone quality scores (qual format) [input]\n";
		print
		  "\t -cluster_seq <txt> - cluster sequences (fasta format) [input]\n";
		print
"\t -singleton_seq <txt> - singleton sequences (fasta format) [input]\n";
		print "\t -clone_id <txt> - clone id [input]\n";
		print "\t -cluster_d <txt> - cluster id [input]\n";
		print
"\t -contig_seq <txt> - assembled/contig sequences (fasta format) [output]\n";
		print
"\t -contig_qual <txt> - assembled/contig quality scores (qual format) [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$project_dir         = realpath($project_dir);
	$cluster_fasta_seq   = realpath($cluster_fasta_seq);
	$singleton_fasta_seq = realpath($singleton_fasta_seq);
	$clone_fasta_seq     = realpath($clone_fasta_seq);
	$clone_qual_scores   = realpath($clone_qual_scores);
	$contig_fasta_seq    = realpath($contig_fasta_seq);
	$contig_qual_scores  = realpath($contig_qual_scores);

	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);

	$base_dir = "$project_dir/assemble_$timestamp";
	%dirs     = (
		"seq_qual" => "$base_dir/seq_qual",
		"clusters" => "$base_dir/clusters",
		"phrap"    => "$base_dir/phrap",
		"contigs"  => "$base_dir/contigs",
		"out"      => "$base_dir/out"
	);

	# Create directories
	&Utils::FileUtils::create_project_directories( $project_dir, $base_dir,
		%dirs );

	# Split clone sequences and quality scores into individual files
	&Utils::FastaUtils::split_fasta( $clone_fasta_seq, $dirs{"seq_qual"} );
	&Utils::QualScoreUtils::split_qual( $clone_qual_scores, $dirs{"seq_qual"} );

	# Split cluster sequences
	&split_clusters( $cluster_fasta_seq, $dirs{"clusters"}, $cluster_id );

	chdir($base_dir);

	# Find phrap executable
	$phrap_exec = &Utils::SystemUtils::find_program("phrap");

	# Assembling (most of the code is from Partigene_predb.pl,
	# need to check licensing)
	# Get number of sequences in each cluster from clusters directory
	my %cluster_content;
	my $no_seqs;
	opendir( DIR, $dirs{"clusters"} )
	  || die " can't open " . $dirs{"clusters"} . "	: $!";
	my $file;
	while ( defined( $file = readdir(DIR) ) ) {
		if ( $file =~ /^$cluster_id/ ) {
			$no_seqs = 0;
			open( INPUT, "< " . $dirs{"clusters"} . "/$file" )
			  || die " can't open " . $dirs{"clusters"} . "/$file" . ": $! ";
			while (<INPUT>) {
				if (/^>/) {
					$no_seqs++;
					$cluster_content{$file} = $no_seqs;
				}
			}
			close INPUT;
		}
	}
	closedir(DIR);

	# Compare number of sequences in each Clobb Cluster (ie new)
	# with number of sequences each phrap (ie old),
	# keep the ones with different numbers for update
	my $key;
	my %cluster_upgrade = %cluster_content;
	foreach $key ( keys(%cluster_upgrade) ) {
		if ( -e $dirs{"seq_qual"} . "/$key" ) {
			my $tmp_cmd = "grep -c \"^>\" " . $dirs{"seq_qual"} . "/$key";
			$no_seqs = `$tmp_cmd`;
			chomp($no_seqs);
			if ( $cluster_upgrade{$key} == $no_seqs ) {
				delete $cluster_upgrade{$key};
			}
		}
	}

	# Create @process_list for clusters to be processed
	my @process_list;
	foreach $key ( keys(%cluster_upgrade) ) {
		push( @process_list, $key );
	}
	my $size = @process_list;
	if ( $size == 0 ) {
		print " Could not find any clusters with new sequences \n ";
		exit();
	}
	else {
		print " \n $size Clusters will be assembled or updated \n ";
	}

	print " \nProcessing files, please wait . \n ";

	foreach $file (@process_list) {
		my $skip         = 1;
		my $seq          = '';
		my $flag         = 0;
		my $header       = '';
		my $title        = '';
		my $two_seq_flag = $skip;

		my $tmp_cmd = "grep -c \">\" " . $dirs{"clusters"} . "/$file";

		my $wc = `$tmp_cmd`;
		chomp $wc;
		if ( $wc == 2 && $skip == 2 ) { $two_seq_flag = 1; }

		open( CLUSFILE, "<" . $dirs{"clusters"} . "/$file" );
		my $seq_count = 0;
	  CLUSGEN: while ( my $line = <CLUSFILE> ) {
			if ( ( $line =~ /^>/ ) || eof(CLUSFILE) ) {
				if ($seq) {

					&prepare_qual( "$file", $header, $seq, $clone_id,
						$dirs{"seq_qual"}, $dirs{"phrap"}, 0 );

				}
				$seq = '';
			}
			if ( $line =~ /^>(.+)/ ) {
				$seq_count++;
				if ( $seq_count > 400 ) {
					last CLUSGEN;
				} # more than 400 sequences in a cluster can lead to memory problems, why need more than 400 anyway ? (A workaround)
				$header = $1;
			}
			if ( $line !~ /^>/ ) { chomp $line; $seq .= $line; }
		}
		close(CLUSFILE);
	}

	# Run phrap
	chdir( $dirs{"phrap"} );
	opendir( DIR, "./" );
	$size = @process_list;
	foreach $file (@process_list) {
		$size--;
		printf( "\r%5d clusters remaining", $size );
		system(
"$phrap_exec $file -forcelevel $DEFAULT_FORCELEVEL -new_ace  2> /dev/null");
		my $contfile = $file;
		$contfile .= ".contigs";

		unless ( ( -e $contfile ) && ( -s $contfile > 0 ) )
		{    # it didn't assemble so use the largest est as the contig
			push( @pseudo_contig, $file );
			open( FP, "$file " );
			my $length   = 0;
			my $clus_seq = '';
			$seq = '';
			while ( $line = <FP> ) {
				if ( $line =~ /^>/ || eof(FP) ) {
					if ( length($seq) > length($clus_seq) ) {
						$clus_seq = $seq;
					}
					$seq = '';
				}
				else {
					$seq .= $line;
				}
			}
			close(FP);

			open( OUT, ">$file.pseudocontig" );
			print OUT ">$file\n$clus_seq";
			close(OUT);
			open( OUT, ">$file.pseudocontig.qual" ); # create dummy quality file
			print OUT ">$file\n ";
			for ( my $i = 1 ; $i <= length($clus_seq) ; $i++ ) {
				print OUT "$DEFAULT_QUALSCORE ";
				if ( $i % 30 == 0 ) { print OUT "\n"; }
			}
			print OUT "\n";
			close(OUT);
		}
	}
	closedir(DIR);
	chdir( $dirs{"phrap"} );

	# Copy .pseudocontigs to .contig
	foreach $file (@pseudo_contig) {
		if ($file) {
			system( "cp "
				  . $dirs{"phrap"}
				  . "/$file.pseudocontig "
				  . $dirs{"phrap"}
				  . "/$file.contigs" );
			system( "cp "
				  . $dirs{"phrap"}
				  . "/$file.pseudocontig.qual "
				  . $dirs{"phrap"}
				  . "/$file.contigs.qual" );
		}
	}

	&process_clusters( ".contigs", ".fsa", $dirs{"phrap"}, $dirs{"contigs"} );
	&process_clusters( ".contigs.qual", ".qlt", $dirs{"phrap"},
		$dirs{"contigs"} );
	&process_singletons( $singleton_fasta_seq, $clone_id, $cluster_id,
		$dirs{"seq_qual"}, $dirs{"contigs"} );

	# Copy all pseudocontigs across to contigs directory
	my @all_pseudos = glob( $dirs{"phrap"} . "/*pseudocontig" );
	foreach $file (@all_pseudos) {
		my $pattern = $dirs{"phrap"};
		$file =~ s/$pattern\///;
		$file =~ s/\.pseudocontig//;
		my $dummyfile = $file . "_1.fsa";
		system( "cp "
			  . $dirs{"seq_qual"}
			  . "/$file.pseudocontig "
			  . $dirs{"contigs"}
			  . "/$dummyfile" );
	}

	# Create log file
	$timestamp = `date +%D+%H:%M`;
	chomp($timestamp);
	$timestamp =~ s/\//-/g;
	my $file_name = "phrap_$cluster_id" . "_$timestamp.txt";
	my @list      = glob( $dirs{"seq_qual"} . "/*.ace" );
	open( LOGFILE, ">$base_dir/$file_name" )
	  || die " can't open $base_dir/$file_name: $! ";
	print LOGFILE "    ### Cluster\tcontig\tnumber of sequences\n";

	foreach my $file (@list) {
		open( ACEFILE, $file ) || die "Can't find $file\n";
		$file =~ s/\.ace//g;
		$file =~ s/phrap\///g;
		while ( my $line = <ACEFILE> ) {
			if ( $line =~ /^CO\s+(Contig\d+)\s+\d+\s+(\d+)/ ) {
				print LOGFILE "$file\t$1\t$2\n";
			}
		}
	}
	print
"\nA report on the phrap based assembly process has been saved in the file:\n";
	print "$base_dir/$file_name\n";
	close(LOGFILE);
	close(ACEFILE);

	# Create contig output files
	system( "cat "
		  . $dirs{"contigs"}
		  . "/*.fsa > "
		  . $dirs{"out"}
		  . "/contigs.fsa" );
	system( "cp " . $dirs{"out"} . "/contigs.fsa $contig_fasta_seq" );
	system( "cat "
		  . $dirs{"contigs"}
		  . "/*.qlt > "
		  . $dirs{"out"}
		  . "/contigs.qual" );
	system( "cp " . $dirs{"out"} . "/contigs.qual $contig_qual_scores" );

}
#############################  End main() ##########################################

######################## Prepare singleton files ####################################
sub process_singletons() {
	my $sinlgetons_file = $_[0];
	my $clone_id        = $_[1];
	my $cluster_id      = $_[2];
	my $seq_qual_dir    = $_[3];
	my $contigs_dir     = $_[4];    # location of sequence quality files

	# Read in sequences one at a time
	my $line;
	my $file;
	open( SINGFILE, "<$sinlgetons_file" );
	my $flag   = 0;
	my $header = '';
	my $seq    = '';
	while ( $line = <SINGFILE> ) {
		chomp $line;
		if ( $line !~ /^>/ ) { $seq .= $line; next; }
		if ( $line =~ /^>/ && $seq ) {
			&prepare_qual( "$file", $header, $seq, $clone_id, $seq_qual_dir,
				$contigs_dir, 1 );
			$seq = '';
		}
		if ( $line =~ /^>(.+)/ ) {
			$header = $1 . "_1";
			$header =~ /($cluster_id\d\d\d\d\d_1)\s*$/
			  ;    # take cluster id as the file name
			$file = $1;
			next;
		}
		if ( $line !~ /^>/ ) { chomp $line; $seq .= $line; next; }
	}
	close(SINGFILE);
	&prepare_qual( "$file", $header, $seq, $clone_id, $seq_qual_dir,
		$contigs_dir, 1 );    # do the last sequence
}
####################### End process_singletons() ####################################

######################## Prepare quality files for phrap ############################
sub prepare_qual() {
	my $outfile      = $_[0];
	my $header       = $_[1];
	my $seq          = $_[2];
	my $clone_id     = $_[3];
	my $seq_qual_dir = $_[4];
	my $phrap_dir    = $_[5];
	my $singleton    = $_[6];
	my $seqfile      = '';

	my ( $in_seq, $line, $n, $qualseq, @qual_scores, @readseq );
	my ( $end, $start, $startl, $endl );
	my ( $seq_length, $beg_low, $end_low );
	my $trace_problem_flag;

	if ( $header =~ /($clone_id\_\w{2,5}\_[a-z0-9]+)/i ) {  # EGTDC style traces
		$seqfile = $1;
		$seqfile =~ /$clone_id\_([A-Za-z0-9]{2,5})\_[a-z0-9]+/i;
	}

	$outfile =~ /\/(.+)$/;    # take cluster id as the file name
	my $title          = $1;
	my $orig_file_flag = 0;

	if ( $singleton == 1 ) {
		open( QUALFILE, ">>$phrap_dir/$outfile.qlt" )
		  ;    # for singletons this will actually be contigs directory
		open( SEQFILE, ">>$phrap_dir/$outfile.fsa" );
		my $new_header = $outfile;
		$new_header =~ s/protein\///;
		print SEQFILE ">$new_header\n";
		print QUALFILE ">$new_header\n";
	}
	else {
		open( QUALFILE, ">>$phrap_dir/$outfile.qual" );
		open( SEQFILE,  ">>$phrap_dir/$outfile" );
		print SEQFILE ">$header\n";
		print QUALFILE ">$header\n";
	}

	if ( -f "$seq_qual_dir/$seqfile.fsa" ) {
		open( FH, "<$seq_qual_dir/$seqfile.fsa" );
		$in_seq = '';
		while ( $line = <FH> ) {
			if ( $line !~ /^>/ ) { chomp $line; $in_seq .= $line; }
		}
		close(FH);

		my $statfile = '';
		if ( -f "$seq_qual_dir/$seqfile.exp" ) {
			$statfile = "$seq_qual_dir/$seqfile.exp";
		}
		elsif ( -f "$seq_qual_dir/$seqfile.qual" ) {
			$statfile = "$seq_qual_dir/$seqfile.qual";
		}
		if ($statfile) {    # have quality file
			$orig_file_flag = 1;
			open( FHA, "<$statfile" );
			while ( $line = <FHA> ) {
				$line =~ s/^AV\s//;
				if ( $line !~ /^>/ ) { $qualseq .= " " . $line; }
			}
			close(FHA);

			@qual_scores = split( ' ', $qualseq );
			$seq_length  = length($in_seq);
			$n           = 0;
			$startl      = 0;
			$endl        = 0;
			@readseq     = split( '', "$in_seq" );
			for ( my $i = $startl ; $i < $seq_length - $endl ; $i++ ) {
				print SEQFILE "$readseq[$i]";
				print QUALFILE "$qual_scores[$i] ";
				$n++;
				if ( $n % 30 == 0 ) { print QUALFILE "\n"; }
				if ( $n % 60 == 0 ) { print SEQFILE "\n"; }
			}
			if ( $n % 30 != 0 ) { print QUALFILE "\n"; }
			if ( $n % 60 != 0 ) { print SEQFILE "\n"; }
			if ( $seq_length - $endl - $startl < 1 ) {
				$orig_file_flag     = 0;
				$trace_problem_flag = 1;
			}
		}
	}

	#  }
	if ( $orig_file_flag == 0 ) {    # couldn't find original files on system
		$n = 0;
		while ( $seq =~ /(.)/g ) {
			print SEQFILE $1;
			print QUALFILE "$DEFAULT_QUALSCORE ";
			$n++;
			if ( $n % 30 == 0 ) { print QUALFILE "\n"; }
			if ( $n % 60 == 0 ) { print SEQFILE "\n"; }
		}
		if ( $n % 30 != 0 ) { print QUALFILE "\n"; }
		if ( $n % 60 != 0 ) { print SEQFILE "\n"; }
		close(QUALFILE);
		close(SEQFILE);
	}

}
############################ End prepare_qual() #####################################

######## Splits the contig and qual files into individual files for decoder #########
# e.g XHC1000_1, XHC1000_2
sub process_clusters() {
	my $suffix      = $_[0];
	my $newsuffix   = $_[1];
	my $phrap_dir   = $_[2];
	my $contigs_dir = $_[3];

	opendir( DIR, "$phrap_dir" );
	my ( $file, $n, $line, $flag, $in_seq, $newfile, $header );
	while ( defined( $file = readdir(DIR) ) ) {
		if ( $file =~ /$suffix$/ && ( -s "$phrap_dir/$file" > 10 ) ) {
			open( FSAFILE, "<$phrap_dir/$file" );
			while ( $line = <FSAFILE> ) {
				if ( $line =~ /^>(.+)/ ) {
					$header  = $1;
					$newfile = $file;
					$newfile =~ s/$suffix//;
					if ( $header =~ /Contig(\d+)/ ) { $newfile .= "_" . $1 }
					$newfile .= $newsuffix;
					open( OUTFILE, ">$contigs_dir/$newfile" );
					$newfile =~ s/$newsuffix//g;
					print OUTFILE ">$newfile\n";
				}
				else {
					print OUTFILE "$line";
				}
			}
			close(FSAFILE);
		}
	}
}
############################ End process_clusters() ##################################

########### Split file with clusters into individual clusters #######################
############ e.g. split_clusters( "clusters.fasta", "clusters" ) ####################
sub split_clusters() {
	my $file       = shift(@_);
	my $dir        = shift(@_);
	my $cluster_id = shift(@_);
	my $flag       = 0;
	my $in_seq     = "";
	my @heading;
	my $header;
	my $filename;

	open( INPUT, "<$file" ) || die "can't open $file: $!";
	while (<INPUT>) {
		if (/^>/) {
			chop;
			$_ =~ /($cluster_id\d\d\d\d\d)\s*$/;
			my $filename = $1;
			open( OUTFILE, ">>$dir/$filename" )
			  || die "can't open $dir/$filename: $!";
			print OUTFILE "$_\n";
			next;
		}
		print OUTFILE $_ unless (/^>/);
		if (/^>/) { close(OUTFILE); }
	}
	close(OUTFILE);
}
############################ End split_clusters() ###################################
