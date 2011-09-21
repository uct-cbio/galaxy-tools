#! /usr/bin/perl -w
# Some parts are similar to the trimming part in trace2dbest.pl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd qw(realpath);
use Utils::FileUtils;
use Utils::SystemUtils;
use Utils::LogUtils;

=head1 NAME

basecall.pl - Unzip compressed file containing ab1 binaries and then do base call using phred. Returns sequences with corresponding quality scores.

=head1 SYNOPSIS

-project <txt> - project directory / project name [input]

-ab1 <txt> - ab1 sequences (zip archive) [input]

-seq <txt> - sequence files (fasta fornat) [output]

-qual <txt> - quality scores (qual format) [output]

-scf <txt> - scf file (scf format) [output]

-phred <txt> - phred output file (txt format) [output]
   
-help - Get more detailed help

=head1 RUNNING

./basecall.pl -project /home/user/project -ab1 /home/user/ab1_sequences.zip -seq /home/user/all.fsa -qual /home/user/all.qual -phred /home/user/all_phred.txt -scf /home/user/all.scf


=head1 OUTPUT

Directory structure:
    project_dir
        basecall
            ab1
            phred
            seq
            qual
            scf
            out
                            
=head1 TODO                            
 
 1) Assume current chromatograph files is in ab1 format. Other formats need to be accommodated.

=cut

########################## main() ##################################################
{
    my $project_dir;
    my $base_dir;
    my %dirs;
    my $log_file;
    my ( $ab1_zipped_input, $phred_output, $fasta_output, $qual_scores_output, $scf_output, $help );
    my ( $phredoptions, $phred_exec );

    GetOptions(
                "project=s" => \$project_dir,
                "ab1seq=s"  => \$ab1_zipped_input,
                "seq=s"     => \$fasta_output,
                "qual=s"    => \$qual_scores_output,
                "phred=s"   => \$phred_output,
                "scf=s"     => \$scf_output,
                "help"      => \$help,
    );

    if (
         ($help)
         || !(
                  ($project_dir)
               && ($phred_output)
               && ($ab1_zipped_input)
               && ($fasta_output)
               && ($qual_scores_output)
               && ($scf_output)
         )
      )
    {
        print
"Unzip compressed file containing ab1 binaries and then do base call using phred. Returns sequences (fasta file) with corresponding quality scores (qual file).\n";
        print "Usage :\n\t basecall.pl <list of arguments> (specify the full path of input and output directories)\n\n";
        print "\t -project <txt> - project directory / project name [input]\n";
        print "\t -ab1 <txt> - ab1 sequences (zip archive) [input]\n";
        print "\t -seq <txt> - sequence files (fasta fornat) [output]\n";
        print "\t -qual <txt> - quality scores (qual format) [output]\n";
        print "\t -phred <txt> - phred output file (txt format) [output]\n";
        print "\t -scf <txt> - scf file (scf format) [output]\n";
        print "\t -help\t - Get more detailed help\n";
        exit();
    }

	$project_dir  = realpath($project_dir);
	$ab1_zipped_input = realpath($ab1_zipped_input);
	$fasta_output = realpath($fasta_output);
	$qual_scores_output = realpath($qual_scores_output);
	$phred_output = realpath($phred_output);
	$scf_output = realpath($scf_output);

    my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
    chomp($timestamp);
    $base_dir = "$project_dir/basecall_$timestamp";
    %dirs = (
              "ab1"  => "$base_dir/ab1",
              "seq"  => "$base_dir/seq",
              "qual" => "$base_dir/qual",
              "phd"  => "$base_dir/phd",
              "scf"  => "$base_dir/scf",
              "out"  => "$base_dir/out"
    );

    # Create basecall directories
    &Utils::FileUtils::create_project_directories( $project_dir, $base_dir, %dirs );

    # Create log file
    $log_file = $base_dir . "/$timestamp.log";
    open( my $log_fd, ">$log_file" );

    # Unzip ab1 binary
    chdir chdir $dirs{"ab1"};
    print "Unzipping...\n";
    &Utils::LogUtils::log_info( $log_fd, "Unzipping..." );
    system "unzip $ab1_zipped_input";
    print "Done.\n";
    &Utils::LogUtils::log_info( $log_fd, "Done" );

    # Set phred options and get executable
    my $phred_log_file = $base_dir . "/phred.log";
    $phredoptions = "-sd "
      . $dirs{"seq"} . " -qd "
      . $dirs{"qual"} . " -pd "
      . $dirs{"phd"} . " -cd "
      . $dirs{"scf"}
      . " -exit_nomatch -trim_alt  \\\"\\\" 2> $phred_log_file 1> /dev/null";
    $phred_exec = &Utils::SystemUtils::find_program("phred");

    # Do base calling using phred
    chdir $dirs{"phd"};
    print "Base calling using phred...\n";
    &Utils::LogUtils::log_info( $log_fd, "Base calling using phred..." );
    system( "$phred_exec -id " . $dirs{"ab1"} . " $phredoptions" );
    print "Done.\n";
    &Utils::LogUtils::log_info( $log_fd, "Done." );

	print "Prepare output files...\n";
	&Utils::LogUtils::log_info( $log_fd, "Prepare output files..." );
    my $seq_file   = $dirs{"out"} . "/seq.fsa";
    my $qual_file  = $dirs{"out"} . "/seq.qual";
    my $phred_file = $dirs{"out"} . "/phred.txt";
    my $scf_file   = $dirs{"out"} . "/seq.scf";

    # Create fasta output file
    system( "cat " . $dirs{"seq"} . "/* > $seq_file" );
    system("sed \"s/\.ab1//g\" -i $seq_file");
    system("cp $seq_file $fasta_output");

    # Create quality scores file
    system( "cat " . $dirs{"qual"} . "/* > $qual_file" );
    system("sed \"s/\.ab1//g\" -i $qual_file");
    system("cp $qual_file $qual_scores_output");

    # Create phred output file
    system( "cat " . $dirs{"phd"} . "/* > $phred_file" );
    system("sed \"s/\.ab1//g\" -i $phred_file");
    system("cp $phred_file $phred_output");

    # Create scf output file
    system( "cat " . $dirs{"scf"} . "/* > $scf_file" );
    system("sed \"s/\.ab1//g\" -i $scf_file");
    system("cp $scf_file $scf_output");
    print "Done.\n";
    &Utils::LogUtils::log_info( $log_fd, "Done." );
    
    close $log_fd;
}
#############################  End main() ##########################################
