#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Utils::FastaUtils;
use Utils::SystemUtils;
use Utils::FileUtils;
use Utils::LogUtils;
use Utils::BlastUtils;
use Cwd;
use Cwd qw(realpath);
use threads;
use Thread::Queue;
use HTML::Template;

=head1 NAME

convert_blast_report.pl - Convert blast report from xml or pairwise to other formats.

=head1 SYNOPSIS

-blast_format <txt> - blast format  (blast default pairwise output (bls) or XML blast output (xml)) [input]

-blast_report <txt> - blast reports  (a concatenated bls or xml file) [input] 

-new_blast_format <txt> - HitTableWriter (htw), HTMLResultWriter (html), GbrowseGFF (gff), HSPTableWriter (hsp), BSMLResultWriter (bsml), ResultTableWriter (rtw), BestHitTable (bht) and pairwise blast report (bls)  (1 file) [input]

-new_blast_report <txt> - output report [ouput]

-prog_work <txt> - program work directory. Can be removed by user after results have been retrieved with the -d flag. [output]

-d - remove program work directory with this flag [input]

-help - Get more detailed help

=head1 RUNNING

./convert_blast_report.pl -blast_format xml -blast_report reports.xml -new_format html -new_report new.html -prog_work out
 
=head1 OUTPUT

Directory structure:
    project_dir
        convert_blast_reports
            reports
            new
            out 
                            
=head1 TODO
 
=cut

########################## main() ##################################################
{
	my $prog_work_dir;
	my $base_dir;
	my %dirs;
	my $log_file;
	my ( $fasta_sequences, $blast_format, $blast_report, $new_format,
		$new_report, $help );
	my $delete_program_work;

	GetOptions(
		"blast_format=s"     => \$blast_format,
		"blast_report=s"     => \$blast_report,
		"new_blast_format=s" => \$new_format,
		"new_blast_report=s" => \$new_report,
		"prog_work=s"        => \$prog_work_dir,
		"d"                  => \$delete_program_work,
		"help"               => \$help,
	);

	if (
		($help)
		|| !(
			   ($prog_work_dir)
			&& ($blast_format)
			&& ($blast_report)
			&& ($new_format)
			&& ($new_report)
		)
	  )
	{
		print
"Convert blast report from xml or pairwise to other formats.\n";
		print
"Usage :\n\t convert_blast_report.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print
"\t -blast_format <txt> - blast format  (blast default pairwise output (bls) or XML blast output (xml)) [input]\n";
		print
"\t -blast_report <txt> -  blast reports  (a concatenated bls or xml file) [input]\n";
		print
"\t -new_blast_format <txt> - HitTableWriter (htw), HTMLResultWriter (html), GbrowseGFF (gff), HSPTableWriter (hsp), BSMLResultWriter (bsml), ResultTableWriter (rtw), BestHitTable (bht)  and pairwise blast report (bls)  (1 file) [input]\n";
		print "\t -new_blast_report <txt> - output report (1 file) [ouput]\n";
		print
"\t -prog_work <txt> - program work direcotry. Can be removed by user after results have been retrieved with the -d flag.\n";
		print "\t -d - remove program work directory\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}

	$prog_work_dir = realpath($prog_work_dir);
	$blast_report  = realpath($blast_report);
	$new_report    = realpath($new_report);
	
	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);

	$base_dir = "$prog_work_dir/convert_blast_report_$timestamp";

	%dirs = (
		"reports"   => "$base_dir/reports",
		"converted" => "$base_dir/converted"
	);

	# Create project work directory
	&Utils::FileUtils::create_project_directories( $prog_work_dir, $base_dir,
		%dirs );

	$log_file = $base_dir . "/$timestamp.log";
	open( my $log_fd, ">$log_file" );

	my @files;
	chdir( $dirs{"reports"} );

	# Check if correct blast format is specified and split blast reports
	if ( $blast_format eq "bls" ) {
		&Utils::BlastUtils::split_pairwise_blast_report( $blast_report,
			$dirs{"reports"} );
		@files = glob("*.bls");

	}
	elsif ( $blast_format eq "xml" ) {
		&Utils::BlastUtils::split_xml_blast_report( $blast_report,
			$dirs{"reports"} );
		@files = glob("*.xml");
	}
	else {
		print
"Input format not supported, currently only blast pairwise (bls)) and blast xml (xml) are supported.\n";
		&Utils::LogUtils::log_error( $log_fd,
"Input format not supported, currently only blast default output (bls) and blast xml (xml) are supported."
		);
		exit();
	}

	my @loop_data = ();

	for my $file (@files) {
		my $file_out = $file;
		$file_out =~ s/$blast_format/$new_format/;
		$file_out = join "/", $dirs{"converted"}, $file_out;

		#		print "Generate new $new_format report for $file...\n";
		&Utils::LogUtils::log_info( $log_fd,
			"Generate new $new_format report for $file..." );

		if ( $blast_format eq "bls" ) {
			if ( $new_format eq "htw" ) {
				my @columns =
				  qw(query_name query_length hit_name hit_length expect hit_description);
				&Utils::BlastUtils::write_hit_table( $file, 'blast', $file_out,
					@columns );

			}
			elsif ( $new_format eq "bht" ) {
				&Utils::BlastUtils::write_best_hit_table( $file, 'blast',
					$file_out );
			}
			elsif ( $new_format eq "hsp" ) {
				&Utils::BlastUtils::write_hsp_table( $file, 'blast',
					$file_out );
			}
			elsif ( $new_format eq "html" ) {
				&Utils::BlastUtils::write_html_blast_report( $file, 'blast',
					$file_out );
			}
			elsif ( $new_format eq "gff" ) {
				&Utils::BlastUtils::write_gff( $file, 'blast', $file_out );
			}
			elsif ( $new_format eq "bsml" ) {
				&Utils::BlastUtils::write_bsml_result( $file, 'blast',
					$file_out );
			}
			elsif ( $new_format eq "rtw" ) {
				&Utils::BlastUtils::write_result_table( $file, 'blast',
					$file_out );
			}
			elsif ( $new_format eq "bls" ) {
				&Utils::BlastUtils::write_pairwise_blast_report( $file, 'blast',
					$file_out );
			}
			else {
				print
"Output format not supported, currently only HitTableWriter (htw), HTMLResultWriter (html), GbrowseGFF (gff), BSMLResultWriter (bsml), ResultTableWriter (rtw), BestHitTable (bht) and pairwise blast report (bls) are supported.\n";
				&Utils::LogUtils::log_error( $log_fd,
"Output format not supported, currently only HitTableWriter (htw), HTMLResultWriter (html), GbrowseGFF (gff), BSMLResultWriter (bsml), ResultTableWriter (rtw), BestHitTable (bht) and pairwise blast report (bls) are supported.\n"
				);
			}
		}
		elsif ( $blast_format eq "xml" ) {
			if ( $new_format eq "htw" ) {
				my @columns =
				  qw(query_name query_length hit_name hit_length expect hit_description);
				&Utils::BlastUtils::write_hit_table( $file, 'blastxml',
					$file_out, @columns );
			}
			elsif ( $new_format eq "bht" ) {
				&Utils::BlastUtils::write_best_hit_table( $file, 'blastxml',
					$file_out );
			}
			elsif ( $new_format eq "hsp" ) {
				&Utils::BlastUtils::write_hsp_table( $file, 'blastxml',
					$file_out );
			}
			elsif ( $new_format eq "html" ) {
				&Utils::BlastUtils::write_html_blast_report( $file, 'blastxml',
					$file_out );
			}
			elsif ( $new_format eq "gff" ) {
				&Utils::BlastUtils::write_gff( $file, 'blastxml', $file_out );
			}
			elsif ( $new_format eq "bsml" ) {
				&Utils::BlastUtils::write_bsml_result( $file, 'blastxml',
					$file_out );
			}
			elsif ( $new_format eq "rtw" ) {
				&Utils::BlastUtils::write_result_table( $file, 'blastxml',
					$file_out );
			}
			elsif ( $new_format eq "bls" ) {
				&Utils::BlastUtils::write_pairwise_blast_report( $file,
					'blastxml', $file_out );
			}
			else {
				print
"Output format not supported, currently only HitTableWriter (htw), HTMLResultWriter (html), GbrowseGFF (gff), BSMLResultWriter (bsml), ResultTableWriter (rtw) and pairwise blast report (bls) are supported.\n";
				&Utils::LogUtils::log_error( $log_fd,
"Output format not supported, currently only HitTableWriter (htw), HTMLResultWriter (html), GbrowseGFF (gff), BSMLResultWriter (bsml), ResultTableWriter (rtw) and pairwise blast report (bls) are supported.\n"
				);
			}

		}
	}

	system( "cat " . $dirs{"converted"} . "/*.$new_format > $new_report" );
	
	close $log_fd;
	
	if ($delete_program_work) {
		system "rm -rf 	$base_dir";
	}
}

#############################  End main() ###########################################

