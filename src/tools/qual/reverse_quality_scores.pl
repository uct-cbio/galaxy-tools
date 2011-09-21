#! /usr/bin/perl
use strict;
use warnings;
use Cwd qw(realpath);
use Getopt::Long;
use Utils::SystemUtils;
use Utils::QualScoreUtils;
use Utils::FileUtils;

=head1 NAME

reverse_quality_scores.pl - Reverse quality scores.

=head1 SYNOPSIS

-project <txt> - project directory / project name [input]

-qual <txt> - quality scores file [input]

-rev_qual <txt> - reversed quality scores file [output]

-help - Get more detailed help

=head1 RUNNING

./reverse_quality_scores.pl -project /home/user/project -qual /home/user/all.qual  -revqual /home/user/all_reverse.qual

=head1 OUTPUT  

Directory structure:
    project_dir
        reverse_quality_scores
            qual            
            reversed_qual
            outp
            
=head1 TODO

=cut

########################## main() ##################################################
{
	my $project_dir;
	my $base_dir;
	my %dirs;
	my $log_file;
	my ( $qual_scores, $reversed_qual_scores, $help );

	GetOptions(
		"project=s"  => \$project_dir,
		"qual=s"     => \$qual_scores,
		"rev_qual=s" => \$reversed_qual_scores,
		"help=s"     => \$help,
	);

	if ( ($help)
		|| !( ($qual_scores) && ($reversed_qual_scores) ) )
	{
		print
"Usage :\n\t reverse_quality_scores.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print "\t -qual <txt> - quality scores file [input]\n";
		print "\t -rev_qual <txt> - reversed quality scores file [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$project_dir          = realpath($project_dir);
	$reversed_qual_scores = realpath($reversed_qual_scores);

	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);

	$base_dir = "$project_dir/reverse_quality_scores_$timestamp";
	%dirs     = (
		"qual"                      => "$base_dir/qual",
		"reversed_qual" => "$base_dir/reversed_qual",
		"out"                      => "$base_dir/out"
	);

	# Create directories
	&Utils::FileUtils::create_project_directories( $project_dir, $base_dir,
		%dirs );


	# Create qual_reversed directories and split qual scores
	&Utils::QualScoreUtils::split_qual( $qual_scores, $dirs{"qual"});

	# Start reversing quality scores
	chdir(  $dirs{"qual"} );
	my @dir;
	@dir = glob("*.qual");
	foreach my $file_name (@dir) {
#		print "Reversing $file_name\n";
		my $reversed_file_name =  $dirs{"reversed_qual"} .  "/$file_name";
		&reverse_scores( $file_name, $reversed_file_name );
	}

	system("cat ". $dirs{"reversed_qual"} ."/* > " . $dirs{"out"} . "/reverse.qual");
	system("cp " . $dirs{"out"} . "/reverse.qual " . $reversed_qual_scores);
}
#############################  End main() ##########################################

########################## Reverse quality scores  ##################################
sub reverse_scores {
	my ( $file_name, $reversed_file_name ) = @_;
	my ( @lines, $header );
	my @scores;

	# read quality score in current order
	open( FILE, "<$file_name" );
	@lines = <FILE>;

	chomp( $lines[0] );
	$header = $lines[0];

	for my $i ( 1 .. @lines - 1 ) {
		chomp( $lines[$i] );
		my $line = $lines[$i];

		#        print "$line\n";
		push( @scores, split( / /, $line ) );
	}
	close(FILE);

	# write reversed scores to file
	open( FILE, ">$reversed_file_name" );
	print FILE "$header\n";

	my $got_newline = 0;

	for my $i ( 0 .. @scores - 1 ) {
		if ( ( $i + 1 ) % 30 ) {
			print FILE "$scores[@scores - $i - 1] ";
			$got_newline = 0;
		}
		else {
			print FILE "$scores[@scores - $i - 1]\n";
			$got_newline = 1;
		}
	}
	if ( !$got_newline ) {
		print FILE "\n";
	}

	close(FILE);
}
############################ End reverse_scores() ###################################

