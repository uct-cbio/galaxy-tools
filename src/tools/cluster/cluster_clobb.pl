#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd qw(realpath);
use Utils::SystemUtils;
use Utils::FastaUtils;
use Utils::FileUtils;

=head1 NAME

cluster.pl - Cluster sequences using CLOBB2.pl and prepare for assembly process. 

=head1 SYNOPSIS

-project <txt> - project directory / project name [input]

-clone_seq <txt> -  clone sequences (fasta format) [input]

-cluster_id <txt> - cluster id [input]

-cluster_seq <txt> - cluster sequences (fasta format) [output]

-singleton_seq <txt> - singleton sequences (fasta format) [output]   

-help - Get more detailed help


=head1 RUNNING

./cluster_clob.pl -project /home/user/project -clone_seq /home/user/clones.fsa -cluster_id XHC -cluster_seq /home/user/cluster.fsa -singleton_seq /home/user/singletons.fsa   

=head1 OUTPUT
 
Directory structure:
    project_dir
        cluster
            sequences
            clusters
            out
                            
=head1 TODO
 
1) If sequences were previoulsy clustered and new sequences are added do a recluster.

2) Create seperate directory for CLOBB processing. Then also link 'sequences' directory to 'seq' directory in base folder.

=cut

########################## main() ##################################################
{
	my $project_dir;
	my $base_dir;
	my %dirs;
	my ( $fasta_file, $cluster_id, $cluster_fasta_seq, $singleton_fasta_seq,
		$help );
	my ( $clobb_exec, $cluster_file, %cluster_content );
	GetOptions(
		"project=s"       => \$project_dir,
		"clone_seq=s"     => \$fasta_file,
		"cluster_id=s"    => \$cluster_id,
		"cluster_seq=s"   => \$cluster_fasta_seq,
		"singleton_seq=s" => \$singleton_fasta_seq,
		"help"            => \$help,
	);
	if (
		($help)
		|| !(
			   ($project_dir)
			&& ($fasta_file)
			&& ($cluster_id)
			&& ($cluster_fasta_seq)
			&& ($singleton_fasta_seq)
		)
	  )
	{
		print
"Cluster sequences using CLOBB2.pl and prepare for assembly process.\n";
		print
"Usage :\n\t cluster.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print "\t -clone_seq <txt> -  clone sequences (fasta format) [input]\n";
		print "\t -cluster_id <txt> - cluster id [input]\n";
		print
		  "\t -cluster_seq <txt> - cluster sequences (fasta format) [output]\n";
		print
"\t -singleton_seq <txt> - singleton sequences (fasta format) [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$project_dir         = realpath($project_dir);
	$fasta_file          = realpath($fasta_file);
	$cluster_fasta_seq   = realpath($cluster_fasta_seq);
	$singleton_fasta_seq = realpath($singleton_fasta_seq);

	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);
	$base_dir = "$project_dir/cluster_$timestamp";
	%dirs     = (
		"seq" => "$base_dir/sequences"
		, # CLOBB needs to have sequence files in specific directory "sequences"
		"clusters" => "$base_dir/clusters",
		"out"      => "$base_dir/out",
		"tmp"      => "$base_dir/tmp"
	);

	# Create directories
	&Utils::FileUtils::create_project_directories( $project_dir, $base_dir,
		%dirs );
	chdir($base_dir);
	&Utils::FastaUtils::split_fasta( $fasta_file, $dirs{"seq"} );
	$clobb_exec = &Utils::SystemUtils::find_program("CLOBB2.pl")
	  ;    # find CLOBB2.pl executable
	system("$clobb_exec $cluster_id");

	# Process results (most of the code is from Partigene_predb.pl,
	# need to check licensing)
	print "\n\nClustering done - now splitting cluster file into individual\n";
	my $file;
	$cluster_file = $base_dir . "/" . $cluster_id . "EST";
	open( INPUT, "<$cluster_file" ) || die "can't open $cluster_file: $!";
	while (<INPUT>) {
		if (/^>/) {
			chop;
			$_ =~ /($cluster_id\d\d\d\d\d)\s*$/;
			$file = $1;
			open( OUTFILE, ">>" . $dirs{"tmp"} . "/$file.fsa" )
			  || die "can't open " . $dirs{"tmp"} . "/$file.fsa $!";
			print OUTFILE "$_\n";
			next;
		}
		print OUTFILE $_ unless (/^>/);
		if (/^>/) { close(OUTFILE); }
	}
	close(OUTFILE);
	close(INPUT);

	opendir( DIR, $dirs{"tmp"} );
	system "touch "
	  . $dirs{"clusters"}
	  . "/singletons.fsa"; # Create a singletons file even if there will be none
	while ( defined( $file = readdir(DIR) ) ) {
		my $no_seqs = 0;
		if ( $file =~ /^$cluster_id/ ) {
			open( INPUT, "<" . $dirs{"tmp"} . "/$file" )
			  || die "can't open " . $dirs{"tmp"} . "/$file $!";
			while (<INPUT>) {
				if (/^>/) { $no_seqs++ }
			}
			$cluster_content{$file} = $no_seqs;
			close INPUT;
			if ( $no_seqs == 1 ) {
				system( "cat "
					  . $dirs{"tmp"}
					  . "/$file >> "
					  . $dirs{"clusters"}
					  . "/singletons.fsa" );
				system( "rm " . $dirs{"tmp"} . "/$file" );
			}
			else {
				system( "mv "
					  . $dirs{"tmp"}
					  . "/$file  "
					  . $dirs{"clusters"}
					  . "/$file" );
			}
		}
	}
	system( "rmdir " . $dirs{"tmp"} );
	my $tmp_cmd = "grep \'^>\' " . $dirs{"clusters"} . "/singletons.fsa\|wc -l";
	my $sing_num = `$tmp_cmd`;
	chomp $sing_num;
	my @clus_num  = glob( $dirs{"clusters"} . "/$cluster_id*" );
	my $clus_num  = scalar @clus_num;
	my $clus_ests = 0;

	foreach my $dummy (@clus_num) {
		$clus_ests += `grep -c \'^>\' $dummy`;
	}
	my $seqs     = $sing_num + $clus_ests;
	my $clus_sum = $sing_num + $clus_num;
	print "\n\nSUMMARY OF CLUSTERING FOR $cluster_id\n";
	print "=============================================================\n";
	printf( "Number of sequences\t\t\t= %9d\n",          $seqs );
	printf( "Total number of clusters\t\t= %9d\n",       $clus_sum );
	printf( "Number of clusters with 1 member\t= %9d\n", $sing_num );
	print "Number of clusters with >1 member\n";
	printf( "(derived from $clus_ests sequences)\t\t= %9d\n", $clus_num );
	print "=============================================================\n";
#### create file with some info
	my $key;
	my $info_file = "CLOBB_$cluster_id\_stats.txt";
	open( OUTFILE, ">$base_dir/$info_file" )
	  || die "can't open $base_dir/$info_file $!";
	print OUTFILE "SUMMARY OF CLUSTERING FOR $cluster_id\n";
	print OUTFILE "=======================================================\n";
	print OUTFILE "Number of sequences \t\t\t=\t$seqs\n";
	print OUTFILE "Total number of clusters\t\t=\t$clus_sum\n";
	print OUTFILE "Number of clusters with 1 member\t= $sing_num\n";
	print OUTFILE "Number of clusters with > 1 member\n";
	print OUTFILE "derived from $clus_ests sequences\t\t= $clus_num\n";
	print OUTFILE "=======================================================\n";

	foreach $key ( sort keys(%cluster_content) ) {
		print OUTFILE "$key $cluster_content{$key}\n";
	}
	close OUTFILE;

	# Create cluster and singleton output files
	system( "mv "
		  . $dirs{"clusters"}
		  . "/singletons.fsa "
		  . $dirs{"out"}
		  . "/singletons.fsa" );
	system( "cp " . $dirs{"out"} . "/singletons.fsa " . $singleton_fasta_seq );
	if ( $clus_num > 0 )
	{ # If clusters exists a cat can be done. Otherwise create an empty file.
		system( "cat "
			  . $dirs{"clusters"} . "/* > "
			  . $dirs{"out"}
			  . "/clusters.fsa" );
		system( "cp " . $dirs{"out"} . "/clusters.fsa " . $cluster_fasta_seq );
	}
	else {
		system("touch  $cluster_fasta_seq");
	}
}
#############################  End main() ##########################################
