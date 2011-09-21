#! /usr/bin/perl -w
#######################################################################################
#Name: dbest_submit.pl
#Author: Gerrit Botha, CBIO, University of Capetown
#gerrit@cbio.uct.ac.za
#
#
#Date created: 18/08/09
#Date last modified:
#
#Usage: ./dbest_submit.pl -project /home/user/project -est_seq /home/user/est.fsa -submit false -dbest_file /home/user/dbest_file.txt
#
#Copyright (C) 2009 Gerrit Botha
#
#This program is free software; you can redistribute it and/or
#modify it under the terms of the GNU General Public License
#as published by the Free Software Foundation; either version 3
#of the License, or (at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#######################################################################################

=head1 NAME

dbest_submit.pl - Prepare dbest submission files. These files can be submitted if specied with -submit true flag.

=head1 SYNOPSIS

-project <txt> - project directory / project name [input]

-est <txt> - est nucleotide sequences (fasta format) [input]

-status <txt> - new/update [input]

-cont_name <txt> - contact name must be identical to contact entry [input]

-citation <txt> - mcitation name must be identical to publication title [input]



-submit - submit to dbEST (default off) [input]

-dbest <txt> - dbEST submission file (dbEST format) [output]

-help - Get more detailed help

=head1 RUNNING

./dbest_submit.pl -project /home/user/project -est /home/user/est.fsa -submit false -dbest /home/user/dbest.txt

=head1 OUTPUT

Directory structure:
    project_dir
        dbest_submit
            seq
            submission_files
            log
            output
 
=head1 TODO

 - 
 
=cut

# How to submit info: http://www.ncbi.nlm.nih.gov/dbEST/how_to_submit.html

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd qw(realpath);
use File::Copy;
use Utils::FileUtils;
use Utils::FastaUtils;
use Utils::LogUtils;

########################## main() ##################################################
{

	my @dbest_dirs =
	  ( "dbest_submit", "est", "submission_files", "log", "output" );

	my $project_dir;
	my $base_dir;
	my $log_file;
	my %dirs;
	my ( $est_file, $submit_flag, $dbest_file, $help );

	GetOptions(
		"project=s" => \$project_dir,
		"est=s"     => \$est_file,
		"submit"    => \$submit_flag,
		"dbest=s"   => \$dbest_file,
		"help"      => \$help
	);

	if ( ($help)
		|| !( ($project_dir) && ($est_file) && ($dbest_file) ) )
	{
		print "Predict peptide sequence from nucleotide sequences.\n";
		print
"Usage :\n\t predict_peptide_seq.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print "\t -est <txt> - est sequences [input]\n";
		print "\t -submit - submit to dbEST (default off) [input]\n";
		print "\t -dbest <txt> - dbEST submission file [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$project_dir = realpath($project_dir);
	$est_file    = realpath($est_file);
	$dbest_file  = realpath($dbest_file);

	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);

	$base_dir = "$project_dir/dbest_submit_$timestamp";
	%dirs     = (
		"est"              => "$base_dir/est",
		"submission_files" => "$base_dir/submission_files",
		"out"              => "$base_dir/out"
	);

	# Create process directory structure
	&Utils::FileUtils::create_project_directories( $project_dir, $base_dir,
		%dirs );

	# Create log file
	$log_file = $base_dir . "/$timestamp.log";
	open( my $log_fd, ">$log_file" );

	# Split clone sequences and quality scores into individual files
	&Utils::FastaUtils::split_fasta( $est_file, $dirs{"est"} );

	print "Create submission files...\n";
	&Utils::LogUtils::log_info( $log_fd, "Create submission files..." );
	&create_submission_files( $dirs{"est"}, $dirs{"submission_files"}, "Y" );

	if ($submit_flag) {
		print "Sumbit to dbEST...\n";
		&Utils::LogUtils::log_info( $log_fd, "Sumbit to dbEST..." );
		&submit();
	}
	else {
		print "Do not sumbit to dbEST...\n";
		&Utils::LogUtils::log_info( $log_fd, "Do not sumbit to dbEST..." );
	}

	my $dbest_submit_file = $dirs{"out"} . "/dbest_submit.txt";

	system( "cat " . $dirs{"submission_files"} . "/* > " . $dbest_submit_file );
	system( "cp $dbest_submit_file " . $dbest_file );
	close $log_fd;
}
#############################  End main() ##########################################

#############################  create_submission_files() ###########################
sub create_submission_files() {

	my ( $est_dir, $dbest_dir, $lib ) = @_;
	my ( $estid, $cloneid, $plateid, $row, $column, $primer, $p_end, $homolog );

	my @dir = glob("$est_dir/*.fsa");

	foreach my $fasta_file (@dir) {
		my $sequence = "";
		open( FASTA, "$fasta_file" );
		$fasta_file =~ /.*\/(.*$)/;
		my $dbest_file = $1;
		$dbest_file =~ s/fsa/txt/;
		open( DBEST, ">$dbest_dir/$dbest_file" );
		my @lines = <FASTA>;
		for ( my $i = 1 ; $i < @lines ; $i++ ) {
			$sequence = $sequence . $lines[$i];
			chomp $sequence;
		}

		$estid   = "Y";
		$cloneid = "Y";
		$plateid = "Y";
		$row     = "Y";
		$column  = "Y";
		$primer  = "Y";
		$p_end   = "Y";
		$homolog = "Y";

		my $dbest = "TYPE: EST
STATUS: New
CONT_NAME: X
CITATION:
X
LIBRARY: $lib
EST#: $estid
CLONE: $cloneid
SOURCE: X
ERROR: 0.0001
PLATE: $plateid
ROW: $row
COLUMN: $column
SEQ_PRIMER: $primer
P_END: $p_end
DNA_TYPE: cDNA
PUBLIC:$homolog
COMMENT:
Institutes contributing to this submission are: University of
Washington (UW), Seattle, USA; South African National Bioinformatics
Institute (SANBI), University of the Western Cape, South Africa and
J. Craig Venter Institute (JCVI), Rockville, USA. Tissue samples were
provided by Southwest Foundation for Biomedical Research (SFBR), San
Antonio, USA.
Additional details are available at ftp://baboon.washington.edu/.
SEQUENCE:
";
		print DBEST $dbest;

		my @seq = split( '', $sequence );
		my $n = 0;

		for ( my $i = 0 ; $i < @seq ; $i++ ) {
			print DBEST "$seq[$i]";
			$n++;
			if ( ( $n % 60 ) == 0 ) {
				print DBEST "\n";
			}
		}
		if ( ( $n % 60 ) != 0 ) {
			print DBEST "\n";
		}

		print DBEST "||\n";
		close DBEST;

	}
}
#############################  end create_submission_files() #######################

#############################  submit() ############################################
sub submit() {

}
#############################  end submit() ########################################

