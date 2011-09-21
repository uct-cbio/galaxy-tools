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

generate_html_blast_reports.pl - Generate HTML BLAST reports from xml or pairwise BLAST reports. An index file is generated that links to the reports.

=head1 SYNOPSIS

-prog_work <txt> - program working directory [output]

-d - is flag is specified the working directory will be deleted when program has completed

-blast_format <txt> - blast format  (blast default pairwise output (bls) or XML blast output (xml)) [input]

-blast_report <txt> - blast reports  (in 1 txt file) [input] 

-html_blast_report <txt> - html reports  (1 html file) [ouput]

-out <txt> - directory for storing single HTML BLAST reports [ouput]

-help - Get more detailed help

=head1 RUNNING

./generate_html_blast_reports.pl -prog_work /home/user/project -blast_format xml -blast_report /home/user/reports.txt -html_blast_report /home/user/reports.html -d
 
=head1 OUTPUT

Directory structure:
    prog_work_dir
        generate_html_blast_reports
            reports
=head1 TODO

=head1 Other

- The HTML pages are generated with the Perl HTML::Template library.

- The HTML BLAST page will contain links to the report file names. Not to the absolute path! This is specifically done so that the output page and links 
 will work within galaxy. Check $galaxy_path variable.
 
=cut

########################## main() ##################################################
{
	my $out_dir;
	my $prog_work_dir;
	my $base_dir;
	my $delete_program_work = '';
	my %dirs;
	my $log_file;
	my ( $fasta_sequences, $blast_format, $blast_report, $html_blast_report, $help );
	my $blast_db_location;    # get from environmental variable
	my $blast_exec;

	GetOptions(
		"out=s"      => \$out_dir,
		"prog_work=s"      => \$prog_work_dir,
		"blast_format=s" => \$blast_format,
		"blast_report=s" => \$blast_report,
		"html_blast_report=s"  => \$html_blast_report,
		"d" => \$delete_program_work,
		"help"           => \$help,
	);

	if (
		($help)
		|| !(
			   ($out_dir)
			&& ($prog_work_dir)
			&& ($blast_format)
			&& ($blast_report)
			&& ($html_blast_report)
		)
	  )
	{
		print
		  "Generate HTML blast reports from xml or pairwise blast reports. An index.html file is generated that link to the reports.\n";
		print
"Usage :\n\t generate_html_blast_reports.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print "\t -blast_report <txt> -  blast reports  (in 1 txt file) [input]\n";
		print "\t -blast_format <txt> - blast format  (blast default pairwise output (bls) or XML blast output (xml)) [input]\n";
	    print "\t -out <txt> - directory for storing single HTML BLAST reports [ouput]\n";
		print "\t -html_blast_report <txt> - html report pages, contains links to single reports (1 html file) [ouput]\n";
		print "\t -prog_work <txt> - program workink directory, where intermediate results and reports will be stored [ouput]\n";
		print "\t -d <flag> - if specified, delete program working directory when done\n";
		print "\t -help\t - Get more detailed help\n";
		exit();
	}
	
	$out_dir  = realpath($out_dir);
	$prog_work_dir  = realpath($prog_work_dir);
	$blast_report  = realpath($blast_report);
	$html_blast_report  = realpath($html_blast_report);
	
	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);
	$base_dir = "$prog_work_dir/generate_blast_html_reports_$timestamp";
	%dirs     = (
		"reports" => "$base_dir/reports",
	);
	
    if ( -e "$out_dir" ) {
        system("rm -rf $out_dir");
        system("mkdir $out_dir");
    } else {
        system("mkdir $out_dir");
    }

	# Create generate_html_blast_reports directories
	&Utils::FileUtils::create_project_directories( $prog_work_dir, $base_dir,
		%dirs );

	# Create log file
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

	open( my $html_blast_report_fd, ">$html_blast_report" );
    print $html_blast_report_fd "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\"\>\n";

	my $html_blast_template_str = &get_html_blast_template();
	my $html_blast_template =
	  HTML::Template->new( scalarref => \$html_blast_template_str);
	$html_blast_template->param( TITLE => 'HTML BLAST Reports' );

	my @loop_data = ();

	for my $file (@files) {
		my $file_out = $file;
		$file_out =~ s/$blast_format/html/;
		$file_out = join "/", $out_dir, $file_out;
		&Utils::LogUtils::log_info( $log_fd,
			"Generate html report for $file..." );
		if ( $blast_format eq "bls" ) {
			&Utils::BlastUtils::write_html_blast_report( $file, 'blast',
				$file_out );
		}
		elsif ( $blast_format eq "xml" ) {
			&Utils::BlastUtils::write_html_blast_report( $file, 'blastxml',
				$file_out );
		}

		my ( %row_data, $real_path, $sequence_name,$relative_path );
		$real_path = $file_out;
		$file_out =~ /\/.*\/(.*)\.html$/;
		$sequence_name = $1;
		my $galaxy_path =  join ".", $sequence_name, 'html';
#		$row_data{SEQUENCE_PATH} = $real_path;
		$row_data{SEQUENCE_PATH} = $galaxy_path;
		$row_data{SEQUENCE_NAME} = $sequence_name;
		push( @loop_data, \%row_data );
	}

	$html_blast_template->param( BLAST_REPORTS => \@loop_data );
	print $html_blast_report_fd  $html_blast_template->output, "\n";
	close $html_blast_report_fd;

	close $log_fd;
	
	if($delete_program_work){
		print "Delete program working directory...\n";	
		system("rm -rf " . $base_dir)
	}
	
	
}

#############################  End main() ###########################################

######################## Get HTML blast template ####################################
sub get_html_blast_template(){
	my $html_blast_template = '
<html>
<head><title><TMPL_VAR NAME=TITLE></title>
<body>
	<h1><TMPL_VAR NAME=TITLE></h1>
	<table cellpadding="3" border="1" cellspacing="1" width="100%">
	<thead>
         <tr>
             <th>BLAST Reports</th>
         </tr>
    </thead>
    <tbody>
	<TMPL_LOOP NAME=BLAST_REPORTS>
		<tr>
		<td><a href="<TMPL_VAR NAME=SEQUENCE_PATH>" name="<TMPL_VAR NAME=SEQUENCE_NAME>"><TMPL_VAR NAME=SEQUENCE_NAME></a></td>
		</tr>
	</TMPL_LOOP>
	</tbody>
</body>
</html>';
	return $html_blast_template;
}
#####################################################################################
