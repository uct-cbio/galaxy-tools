#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd qw(realpath);
use Utils::FileUtils;

=head1 NAME

diff.pl - Find differences between 2 files. Calls the linux diff command. If there is no difference between the 2 files an empty file is returned.

=head1  SYNOPSIS

-file_1 <txt> - file1 [input]\n

-file_2 <txt> - file2 [input]\n

-diff_file <txt> - file showing differences [output]

-help\t - Get more detailed help

=head1 RUNNING

./diff.pl -file_1 /home/user/file1 -file_2 /home/user/file_12 -diff_file /home/user/diff_file

                    
=cut

########################## main() ##################################################
{
	my ( $file_1, $file_2, $diff_file, $help );
	GetOptions(
		"file_1=s"    => \$file_1,
		"file_2=s"    => \$file_2,
		"diff_file=s" => \$diff_file,
		"help"        => \$help,
	);

	if ( ($help)
		|| !( ($file_1) && ($file_2) && ($diff_file) ) )
	{
		print
"Find differences between 2 files. Calls the linux diff command. If there is no difference between the 2 files an empty file is returned.\n";
		print "Usage :\n\t diff.pl <list of arguments>\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print "\t -file_1 <txt> - file1 [input]\n";
		print "\t -file_2 <txt> - file2 [input]\n";
		print "\t -diff_file <txt> - concatenated file [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$file_1    = realpath($file_1);
	$file_2    = realpath($file_2);
	$diff_file = realpath($diff_file);

	system("diff $file_1 $file_2 > $diff_file");
}
#############################  End main() ##########################################
