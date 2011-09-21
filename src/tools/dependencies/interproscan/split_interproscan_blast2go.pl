#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Utils::InterProScanUtils;
use Cwd;
use Cwd qw(realpath);
use threads;
use Thread::Queue;

=head1 NAME

split_interproscan_blast2go.pl - Split interproscan xml files into individual files. Files are formatted to be used as input for blast2go. 

=head1 SYNOPSIS

-interproscan_xml <txt> - interproscan xml file OR  [input]

-interproscan_dir <txt> - directory with interproscan xml files [input]

-blast2go_xml <txt> - blast2go xml output directory [output]

-help - Get more detailed help

=head1 RUNNING

./split_interproscan_blast2go.pl -interproscan_xml /home/user/interproscan.xml -blast2go_xml /home/user/blast2go_xml

=cut

########################## main() ##################################################
{
	my ( $interproscan_xml, $interproscan_dir, $blast2go_xml_dir, $help );

	GetOptions(
		"interproscan_xml=s" => \$interproscan_xml,
		"interproscan_dir=s" => \$interproscan_dir,
		"blast2go_xml=s"     => \$blast2go_xml_dir,
		"help"               => \$help,
	);

	if ( ($help)
		|| !(  ( $interproscan_xml || $interproscan_dir )
			&& ($blast2go_xml_dir) ) )
	{
		print
"Split interproscan xml files into individual files. Files are formatted to be used as input for blast2go.\n";
		print
"Usage :\n\t split_interproscan_blast2go.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print
		  "\t -interproscan_xml <txt> - interproscan xml file [input] OR \n";
		print
"\t -interproscan_dir <txt> - directory with interproscan xml files [input]\n";
		print
		  "\t -blast2go_xml <txt> - blast2go xml output directory [output]\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$blast2go_xml_dir = realpath($blast2go_xml_dir);

	if ( -e "$blast2go_xml_dir" ) {
		system("rm -rf $blast2go_xml_dir");
		system("mkdir $blast2go_xml_dir");

	}
	else {
		system("mkdir $blast2go_xml_dir");
	}

	if ( defined($interproscan_xml) ) {
		$interproscan_xml = realpath($interproscan_xml);
		&Utils::InterProScanUtils::split_interproscan( $interproscan_xml,
			$blast2go_xml_dir );
	}

	if ( defined($interproscan_dir) ) {
		$interproscan_dir = realpath($interproscan_dir);
		chdir $interproscan_dir;
		my @interproscan_xmls = glob("*");
		foreach $interproscan_xml (@interproscan_xmls) {
			$interproscan_xml = getcwd() . "/" . $interproscan_xml;
			&Utils::InterProScanUtils::split_interproscan( $interproscan_xml,
				$blast2go_xml_dir );
		}
	}

}
#############################  End main() ###########################################
