#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd qw(realpath);
use Utils::FileUtils;

=head1 NAME

concat_files.pl - Concatenate 2 text files.

=head1  SYNOPSIS

-project <txt> - project directory / project name [input]

-file_1 <txt> - file1 [input]\n

-file_2 <txt> - file2 [input]\n

-concat_file <txt> - concatenated file [output]

-help\t - Get more detailed help

=head1 RUNNING

./concat_files.pl -project /home/user/project -file1 /home/user/file1 -file2 /home/user/file2 -concat_file /home/user/concatfile


=head1 OUTPUT

Directory structure:
    project_dir
        concat_files
            output
                        
=head1 TODO                    
                    
=cut

########################## main() ##################################################
{
	my $project_dir;
	my $base_dir;
	my %dirs;
	my $log_file;
	my ( $file_1, $file_2, $concat_file, $help );
	my @concat_files_dirs = ( "concat_files", "output" );
	GetOptions(
		"project=s"     => \$project_dir,
		"file_1=s"      => \$file_1,
		"file_2=s"      => \$file_2,
		"concat_file=s" => \$concat_file,
		"help"          => \$help,
	);

	if ( ($help)
		|| !( ($project_dir) && ($file_1) && ($file_2) && ($concat_file) ) )
	{
		print "Concatenate 2 files.\n";
		print "Usage :\n\t concat_files.pl\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print "\t -file_1 <txt> - file1 [input]\n";
		print "\t -file_2 <txt> - file2 [input]\n";
		print "\t -concat_file <txt> - concatenated file [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$project_dir = realpath($project_dir);
	$file_1      = realpath($file_1);
	$file_2      = realpath($file_2);
	$concat_file = realpath($concat_file);

	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);

	$base_dir = "$project_dir/concat_files_$timestamp";
	%dirs = ( "out" => "$base_dir/seq", );

	# Create blast directories
	&Utils::FileUtils::create_project_directories( $project_dir, $base_dir,
		%dirs );

	# Concatenate files
	system( "cat $file_1 $file_2 > " . $dirs{"out"} . "/concat_file" );
	system( "cp " . $dirs{"out"} . "/concat_file $concat_file" );
}
#############################  End main() ##########################################
