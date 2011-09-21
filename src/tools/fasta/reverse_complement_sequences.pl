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

reverse_complement_sequences.pl - Reverse complement sequences. Needs revcomp to run.

=head1 SYNOPSIS

-project <txt> - project directory / project name [input]

-seq <txt> - sequences (fasta format) [input]

-revcomp_seq <txt> - reverse complemented sequences (fasta format) [output]

-help - Get more detailed help

=head1 RUNNING

./reverse_complement_sequences.pl -project /home/user/project -seq /home/user/all.fsa  -revcompseq /home/user/all_reverse_complemented.fsa

=head1 OUTPUT

Directory structure:
    project_dir
        reverse_complement_sequences
            seq
            reverse_complemented_seq
            out
                            
=head1 TODO
    
=cut

########################## main() ##################################################
{
	my $project_dir;
	my $base_dir;
	my %dirs;
	my $log_file;
	my ( $seq, $rev_comp_seq, $help );
	my ($revcomp_exec);
	
	GetOptions(
		"project=s"     => \$project_dir,
		"seq=s"         => \$seq,
		"revcomp_seq=s" => \$rev_comp_seq,
		"help"          => \$help,
	);

	if ( ($help)
		|| !( ($project_dir) && ($seq) && ($rev_comp_seq) ) )
	{
		print "Reverse complemented sequences.\n";
		print
"Usage :\n\t basecall.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print "\t -seq <txt> - sequences (fasta format) [input]\n";
		print
"\t -revcomp_seq <txt> - reverse complemented sequences (fasta format) [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}
	$project_dir  = realpath($project_dir);
	$seq          = realpath($seq);
	$rev_comp_seq = realpath($rev_comp_seq);
	
	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);

	$base_dir = "$project_dir/reverse_complement_sequences_$timestamp";
	%dirs     = (
		"seq"     => "$base_dir/seq",
		"reverse_complemented_seq" => "$base_dir/reverse_complemented_seq",
		"out"     => "$base_dir/out"
	);

	# Create directories
	&Utils::FileUtils::create_project_directories( $project_dir, $base_dir,
		%dirs );

	# Create reverse complement directories
	&Utils::FastaUtils::split_fasta( $seq, $dirs{"seq"} );

	# Run revcomp
	$revcomp_exec = &Utils::SystemUtils::find_program("revcomp");
	system("$revcomp_exec -i " . $dirs{"seq"} . "/* -o " . $dirs{"reverse_complemented_seq"});
	system("cat " . $dirs{"reverse_complemented_seq"} . "/* > " . $dirs{"out"} . "/reverse_complemented.fsa");
	system("cp " . $dirs{"out"} . "/reverse_complemented.fsa " . $rev_comp_seq);
}
#############################  End main() ##########################################
