#! /usr/bin/perl -w
#######################################################################################
#Name: predict_peptide_seq.pl
#Author: Gerrit Botha, CBIO, University of Cape Town
#gerrit@cbio.uct.ac.za
#
#
#Date created: 14/07/09
#Date last modified:
#
#Usage: ./predict_peptide_seq.pl -project /home/user/project -seqn /home/user/contig.fsa -qualn /home/user/contig.qual -bias "Oryza sativa" -reports /home/user/reports.txt -seqp /home/user/contig_peptide.fsa
#
#Copyright (C) 2009 Gerrit Botha
#Code from prot4EST_2.2.pl created by James Wasmuth, ICAPB, University of Edinburgh was copied and modified. Copyright (C) 2004,2005 James Wasmuth.
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

use strict;
use warnings;
use Getopt::Long;
use Cwd;
use Cwd qw(realpath);
use File::Copy;
use Utils::FileUtils;
use Utils::SystemUtils;
use Utils::FastaUtils;
use Utils::QualScoreUtils;
use Utils::BlastUtils;
use Utils::LogUtils;
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;
use Bio::SeqIO;

use p4e2::tile_path;
use p4e2::getDNA;
use p4e2::estscan_plugin;

=head1 NAME

predict_peptide_seq.pl - Predict peptide sequence from nucleotide sequences. Needs quality score and protein blast reports.

=head1 SYNOPSIS

-project <txt> - project directory / project name [input]

-seqn <txt> - nucleotide sequences (fasta format) [input]

-nuc <txt> - nuclear genetic code [input] 

-mito <txt> - mitochondrial genetic code [input]

-rrna_similarities <txt> - BLAST against ribosomal RNA genes and search for similarities (yes/no) [input]

-mito_similarities <txt> - BLAST against mithocondrial genes and search for similarities (yes/no) [input]

-nuclear_similarities <txt> - Provide protein BLASTx reports and search for nuclear similarities (yes/no) [input] 

-blast_report <txt> - blast reports  (combined reports in 1 file) (none = protein BLASTx reports not provided)[input] 

-blast_format <txt> - blast format  (blast default pairwise output (bls) or XML blast output (xml)) (none = protein BLASTx reports not provided) [input] 

-estscan <txt> - make use of ESTScan (yes/no) [input]

-estscan_matrix <txt> ESTScan score matrix (hs = Homo sapiens, ht = Arabidopsis thaliana, os = Oryza sativa, dm = Drosophila melanogaster, mm = Mus musculus, rn = Rattus norvegicus, zm = Zea mays, none = if ESTScan is not used) [input]

-decoder <txt> - make use of decoder (yes/no) [input]

-qualn <txt> - nucleotide sequence qualtity scores (qual format)  (none = if decoder is not used) [input]

-bias <txt> - select bias codon table of specific organism e.g. (Oryza sativa). The codon table will be fetched from http://www.kazusa.or.jp. (none = if decoder is not used). [input]

-seqp <txt> - predicted peptide sequence [output]

-seqr <txt> - predicted ribosomal RNA genes [output]

-help - Get more detailed help


NCBI Genetic code table (for nuclear and mitochondrial select option)

1 =  Standard
 
2 =  Vertebrate Mitochondrial

3 =  Yeast Mitochondrial

4 =  Mold, Protozoan,  Coelenterate Mitochondrial and Mycoplasma/Spiroplasma

5 =  Invertebrate Mitochondrial

6 =  Ciliate Macronuclear and Dasycladacean

9 =  Echinoderm Mitochondrial

10 = Euplotid Nuclear

11 = Bacterial

12 = Alternative Yeast Nuclear

13 = Ascidian Mitochondrial

14 = Flatworm Mitochondrial

15 = Blepharisma Macronuclear

16 = Chlorophycean Mitochondrial

21 = Trematode Mitochondrial

22 = Scenedesmus obliquus

23 = Thraustochytrium Mitochondrial


=head1 RUNNING

./predict_peptide_seq.pl -project /home/user/project -nseq /home/user/contig.fsa -nuc 1 -mito 16 -rrna_similarities "yes"  -mito_similarities "yes" -nuclear_similarities "yes"  -blast_report /home/user/reports.txt -blast_format xml  -estscan yes -estscan_matrix hs -decoder yes -nqual /home/user/contig.qual -bias "Oryza sativa" -seqp /home/user/contig_peptide.fsa -seqr /home/user/rRNA.fsa

=head1 OUTPUT

Directory structure:
    project_dir
        predict_peptide_seq
            nseq_nqual
            reports
            estscan
            decoder
            out
 
=head1 TODO

 - Check decoder quality files.
 
 - Validate sequences (p4e2::p4e_checker::sequence($mitofile, 'protein'))
 
 - Not sure why this flag.  -alternative        boolean    [N] The default definition of frame '-1' is the reverse-complement of the set of codons
  used in frame 1. (Frame -2 is the set of codons used by frame 2, similarly frames -3 and 3). This is a common standard, used by
  the Staden package and other programs. If you prefer to define frame '-1' as using the set of codons starting with the last codon
  of the sequence, then set this to be true.
 
 - Check if blast report file is correct 
 
 - blastn and blastx evalue are hard coded

 - get latest tables and necessary files from NCBI
 
 - made a change in tile_path.pm. just initialized variables $ext_3 and $ext_5 after defined. Otherwise there is warnings.
 
 - in function setup_decoder set initial value of i = 0 and changed reference to $$seq now: for ( my $i = 0 ; $i < $$seqs ; $i++ ) was: for ( my $i; $i < $#{$seqs} ; $i++ )) 
 
 - took out $er_msg (my $er_msg="\nPlease report this error to me at nematodes.bioinf\@ed.ac.uk\n
   Please send me the .log file and the message given on the terminal screen.\nI'll get back to you as quickly as possible\n\n";
   
 - in FILENAME_CHANGE, change  my @dir        = glob("./") to my @dir        = glob("*"); Otherwise you couldnot specify any other extensions
 then .seq and .fsa
 
 - at the moment copy peptide sequence with extensions to output sequences. This sequence include "X"s instead of "*"'s. Or give the user the option? 
 
 - place environment variables in xml file?
 
 - generate own estscan matrices
 
=cut

########################## main() ##################################################
{
	my $project_dir;
	my $base_dir;
	my %dirs;
	my $log_file;
	my (
		$seqn_file,      $qualn_file,   $nuc_code,       $mito_code,
		$use_rrna_sim,   $use_mito_sim, $use_nuc_sim,    $blastx_report,
		$blast_format,   $use_estscan,  $estscan_matrix, $use_decoder,
		$bias_orgasnism, $seqp_file, $seqr_file, $help
	);
	my ( $rrna_db, $mito_db, $gen_code_file, $estscan_matrices_dir );
	my $BLASTN_EVAL = 1e-65;
	my $BLASTX_EVAL = 1e-8;
	my $DECODER_BIASTABLE =
	  "human_high.cod"; # change the file name depending on what decoder expects
	my $SEQSUFF         = "fsa";
	my $QUALSUFF        = "qual";
	my $hsp_ref_counter = 0;

	GetOptions(
		"project=s"              => \$project_dir,
		"seqn=s"                 => \$seqn_file,
		"nuc=s"                  => \$nuc_code,
		"mito=s"                 => \$mito_code,
		"rrna_similarities=s"    => \$use_rrna_sim,
		"mito_similarities=s"    => \$use_mito_sim,
		"nuclear_similarities=s" => \$use_nuc_sim,
		"blast_report=s"         => \$blastx_report,
		"blast_format=s"         => \$blast_format,
		"estscan=s"              => \$use_estscan,
		"estscan_matrix=s"       => \$estscan_matrix,
		"decoder=s"              => \$use_decoder,
		"bias=s"                 => \$bias_orgasnism,
		"qualn=s"                => \$qualn_file,
		"seqp=s"                 => \$seqp_file,
		"seqr=s"                 => \$seqr_file,
		"help"                   => \$help,
	);

	if (
		($help)
		|| !(
			   ($project_dir)
			&& ($seqn_file)
			&& ($nuc_code)
			&& ($mito_code)
			&& ($use_rrna_sim)
			&& ($use_mito_sim)
			&& ($use_nuc_sim)
			&& ($blastx_report)
			&& ($blast_format)
			&& ($use_estscan)
			&& ($estscan_matrix)
			&& ($use_decoder)
			&& ($bias_orgasnism)
			&& ($qualn_file)
			&& ($seqp_file)
			&& ($seqr_file)
		)
	  )
	{
		print "Predict peptide sequence from nucleotide sequences.\n";
		print
"Usage :\n\t predict_peptide_seq.pl <list of arguments> (specify the full path of input and output directories)\n\n";
		print "\t -project <txt> - project directory / project name [input]\n";
		print "\t -seqn <txt> - nucleotide sequences (fasta format) [input]\n";
		print "\t -nuc  <txt> - nuclear genetic code [input]\n";
		print "\t -mito <txt> - mitochondrial genetic code [input]\n";
		print
"\t -rrna_similarities <txt> - BLAST against ribosomal RNA genes and search for similarities (yes/no) [input]\n";
		print
"\t -mito_similarities <txt> - BLAST against mithocondrial genes and search for similarities (yes/no) [input]\n";
		print
"\t -nuclear_similarities <txt> - Provide protein BLASTx reports and search for nuclear similarities (yes/no) [input] \n";
		print
"\t -blast_report <txt> - blastx report  (combined reports in 1 txt file) [input]\n";
		print
"\t -blast_format <txt> - blast format  (blast default pairwise output (bls) or xml blast output (xml)) [input] \n";
		print "\t -estscan <txt> - make use of estscan (yes/no) [input]\n";
		print
"\t -estscan_matrix <txt> - -estscan_matrix <txt> ESTScan score matrix (hs = Homo sapiens, ht = Arabidopsis thaliana, os = Oryza sativa, dm = Drosophila melanogaster, mm = Mus musculus, rn = Rattus norvegicus, zm = Zea mays, none = if ESTScan is not used) [input]\n";
		print "\t -decoder <txt> - make use of decoder (yes/no) [input]\n";
		print
"\t -qualn <txt> - nucleotide sequence qualtity scores (qual format) (none = if decoder is not used)  [input]\n";
		print
"\t -bias <txt> - select bias codon table of specific organism e.g. (Oryza sativa). The codon table will be fetched from http://www.kazusa.or.jp. (none = if decoder is not used). [input]\n";
		print "\t -seqp <txt> - predicted peptide sequence [output]\n";
		print "\t -seqr <txt> - predicted ribosomal RNA genes  [output]\n";
		print "\t -help\t - Get more detailed help\n\n";
		print
		  "NCBI Genetic code table (for nuclear and mitochondrial select option)
\t 1 =  Standard 
\t 2 =  Vertebrate Mitochondrial
\t 3 =  Yeast Mitochondrial
\t 4 =  Mold, Protozoan,  Coelenterate Mitochondrial and Mycoplasma/Spiroplasma
\t 5 =  Invertebrate Mitochondrial
\t 6 =  Ciliate Macronuclear and Dasycladacean
\t 9 =  Echinoderm Mitochondrial
\t 10 = Euplotid Nuclear
\t 11 = Bacterial
\t 12 = Alternative Yeast Nuclear
\t 13 = Ascidian Mitochondrial
\t 14 = Flatworm Mitochondrial
\t 15 = Blepharisma Macronuclear
\t 16 = Chlorophycean Mitochondrial
\t 21 = Trematode Mitochondrial
\t 22 = Scenedesmus obliquus
\t 23 = Thraustochytrium Mitochondrial\n\n";
		exit();
	}

	$project_dir   = realpath($project_dir);
	$seqn_file     = realpath($seqn_file);
	$blastx_report = realpath($blastx_report);
	$qualn_file    = realpath($qualn_file);
	$seqp_file     = realpath($seqp_file);
	$seqr_file     = realpath($seqr_file);

	my $timestamp = `date +%Y-%m-%d_%H_%M_%S_%N`;
	chomp($timestamp);
	$base_dir = "$project_dir/predict_peptide_seq_$timestamp";
	%dirs     = (
		"seqn_qualn" => "$base_dir/seqn_qualn",
		"reports"    => "$base_dir/reports",
		"decoder"    => "$base_dir/decoder",
		"estscan"    => "$base_dir/estscan",
		"process"    => "$base_dir/process",
		"out"        => "$base_dir/out"
	);

	# Get enviromental variables
	$rrna_db              = $ENV{"RRNA_DB"};
	$mito_db              = $ENV{"MITO_DB"};
	$gen_code_file        = $ENV{"GEN_CODE_FILE"};
	$estscan_matrices_dir = $ENV{"ESTSCAN_MATRICES_DIR"};

	# Create process directory structure
	&Utils::FileUtils::create_project_directories( $project_dir, $base_dir,
		%dirs );

	# Create log file
	$log_file = $base_dir . "/$timestamp.log";
	open( my $log_fd, ">$log_file" );

	# Check if rRNA blast similarities should be searched
	print "Check if rRNA blast similarities should be searched... ";
	&Utils::LogUtils::log_info( $log_fd,
		"Check if rRNA blast similarities should be searched..." );
	if ( $use_rrna_sim eq 'yes' ) {
		print "yes\n";
		&Utils::LogUtils::log_info( $log_fd, "yes" );
	}
	elsif ( $use_rrna_sim eq 'no' ) {
		print "no\n";
		&Utils::LogUtils::log_info( $log_fd, "no" );
		$use_rrna_sim = undef;
	}
	else {
		print "unkown option: $use_rrna_sim. Please specify (yes/no).\n";
		&Utils::LogUtils::log_info( $log_fd,
			"unkown option: $use_rrna_sim. Please specify (yes/no).\n" );
		exit();
	}

	# Check if mitochodrial blast similarities should be searched
	print "Check if mitochodrial blast similarities should be searched... ";
	&Utils::LogUtils::log_info( $log_fd,
		"CCheck if mitochodrial blast similarities should be searched..." );
	if ( $use_mito_sim eq 'yes' ) {
		print "yes\n";
		&Utils::LogUtils::log_info( $log_fd, "yes" );
	}
	elsif ( $use_mito_sim eq 'no' ) {
		print "no\n";
		&Utils::LogUtils::log_info( $log_fd, "no" );
		$use_mito_sim = undef;
	}
	else {
		print "unkown option: $use_mito_sim. Please specify (yes/no).\n";
		&Utils::LogUtils::log_info( $log_fd,
			"unkown option: $use_mito_sim. Please specify (yes/no).\n" );
		exit();
	}

	# Check if protein BLASTx reports are provided
	print "Check if protein BLASTx reports are provided... ";
	&Utils::LogUtils::log_info( $log_fd,
		"Check if protein BLASTx reports are provided..." );
	if ( $use_nuc_sim eq 'yes' ) {
		print "yes\n";
		&Utils::LogUtils::log_info( $log_fd, "yes" );

		# Check if correct blast format is specified and split blast reports
		print "Check if correct blast format is specified... ";
		&Utils::LogUtils::log_info( $log_fd,
			"Check if correct blast format is specified..." );
		if ( $blast_format eq "bls" ) {
			print "yes: normal blast pairwise\n";
			&Utils::LogUtils::log_info( $log_fd, "yes: normal blast pairwise" );
			print "Split blast reports into individual files...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"Split blast reports into individual files..." );
			&Utils::BlastUtils::split_pairwise_blast_report( $blastx_report,
				$dirs{"reports"} );
		}
		elsif ( $blast_format eq "xml" ) {
			print "yes: blast xml\n";
			&Utils::LogUtils::log_info( $log_fd, "yes: blast xml" );
			print "Split blast reports into individual files...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"Split blast reports into individual files......" );
			&Utils::BlastUtils::split_xml_blast_report( $blastx_report,
				$dirs{"reports"} );
		}
		else {
			print "no\n";
			&Utils::LogUtils::log_info( $log_fd, "no" );
			print
"Output format not supported, currently only blast pairwise (bls)) and blast xml (xml) are supported.\n";
			&Utils::LogUtils::log_error( $log_fd,
"Output format not supported, currently only blast default output (bls) and blast xml (xml) are supported."
			);
			exit();
		}
	}
	elsif ( $use_nuc_sim eq 'no' ) {
		print "no\n";
		&Utils::LogUtils::log_info( $log_fd, "no" );
		$use_nuc_sim = undef;
	}
	else {
		print "unkown option: $use_nuc_sim. Please specify (yes/no).\n";
		&Utils::LogUtils::log_info( $log_fd,
			"unkown option: $use_nuc_sim. Please specify (yes/no)." );
		exit();
	}

	# Check if estscan is selected
	print "Check if ESTScan is selected... ";
	&Utils::LogUtils::log_info( $log_fd, "Check if ESTScan is selected..." );
	if ( $use_estscan eq 'yes' ) {
		print "yes\n";
		&Utils::LogUtils::log_info( $log_fd, "yes" );

		# Check if estscan matrix exists
		if ( $estscan_matrix eq "hs" ) {
			$estscan_matrix = join "/", $estscan_matrices_dir,
			  "Homo_sapiens.smat";
		}
		elsif ( $estscan_matrix eq "at" ) {
			$estscan_matrix = join "/", $estscan_matrices_dir,
			  "Arabidopsis_thaliana.smat";
		}
		elsif ( $estscan_matrix eq "os" ) {
			$estscan_matrix = join "/", $estscan_matrices_dir,
			  "Oryza_sativa.smat";
		}
		elsif ( $estscan_matrix eq "dm" ) {
			$estscan_matrix = join "/", $estscan_matrices_dir,
			  "Drosophila_melanogaster.smat";
		}
		elsif ( $estscan_matrix eq "mm" ) {
			$estscan_matrix = join "/", $estscan_matrices_dir,
			  "Mus_musculus.smat";
		}
		elsif ( $estscan_matrix eq "rn" ) {
			$estscan_matrix = join "/", $estscan_matrices_dir,
			  "Rattus_norvegicus.smat";
		}
		elsif ( $estscan_matrix eq "zm" ) {
			$estscan_matrix = join "/", $estscan_matrices_dir, "Zea_mays.smat";
		}
		elsif ( $estscan_matrix eq "none" ) {
			$estscan_matrix = undef;
		}
		else {
			print
"ESTScan scoring matrix do not exist, currently only Homo sapiens (hs), Arabidopsis thaliana (at), Oryza sativa (os), Drosophila melanogaster (dm), Mus musculus (mm), Rattus norvegicus (rn) and Zea mays (zm) exsist.\n";
			&Utils::LogUtils::log_error( $log_fd,
"ESTScan scoring matrix do not exist, currently only Homo sapiens (hs), Arabidopsis thaliana (at), Oryza sativa (os), Drosophila melanogaster (dm), Mus musculus (mm), Rattus norvegicus (rn) and Zea mays (zm) exsist."
			);
			exit();
		}
	}
	elsif ( $use_estscan eq 'no' ) {
		print "no\n";
		&Utils::LogUtils::log_info( $log_fd, "no" );
		$use_estscan    = undef;
		$estscan_matrix = undef;
	}
	else {
		print "unkown option: $use_estscan. Please specify (yes/no).\n";
		&Utils::LogUtils::log_info( $log_fd,
			"unkown option: $use_estscan. Please specify (yes/no).\n" );
		exit();
	}

	# Check if decoder is selected
	print "Check if DECODER is selected... ";
	&Utils::LogUtils::log_info( $log_fd, "Check if DECODER is selected..." );
	if ( $use_decoder eq 'yes' ) {
		print "yes\n";
		&Utils::LogUtils::log_info( $log_fd, "Yes" );

		# Split quality scores into individual files
		&Utils::QualScoreUtils::split_qual( $qualn_file, $dirs{"seqn_qualn"} );
	}
	elsif ( $use_decoder eq 'no' ) {
		print "no\n";
		&Utils::LogUtils::log_info( $log_fd, "No" );
		$use_decoder    = undef;
		$qualn_file     = undef;
		$bias_orgasnism = undef;
	}
	else {
		print "unkown option: $use_decoder. Please specify (yes/no).\n";
		&Utils::LogUtils::log_info( $log_fd,
			"unkown option: $use_decoder. Please specify (yes/no).\n" );
		exit();
	}

	# Split consensus sequences into individual files
	print "Split nucleotide sequences into individual files...\n";
	&Utils::LogUtils::log_info( $log_fd,
		"Split nucleotide sequences into individual files..." );
	&Utils::FastaUtils::split_fasta( $seqn_file, $dirs{"seqn_qualn"} );

	# Start predicting protein sequencess
	print "Now predict protein sequences...\n";
	&Utils::LogUtils::log_info( $log_fd, "Now predict protein sequences..." );

	# Check executables
	my $blastall_exec = &Utils::SystemUtils::find_program("blastall");
	my $transeq_exec  = &Utils::SystemUtils::find_program("transeq");
	my $formatdb_exec = &Utils::SystemUtils::find_program("formatdb");
	my $decoder_exec  = &Utils::SystemUtils::find_program("decoder");
	my $estscan_exec  = &Utils::SystemUtils::find_program("estscan");

	# Sequences that failed to compare codon with NCBI standard genetic code
	my $failed_codon_comp = "";

	# Extract the individual seed sequences from the input file
	my (%seeds);
	my $total_count = 0;
	my $inSeqIO =
	  Bio::SeqIO->new( -file => "$seqn_file", '-format' => 'Fasta' );

	while ( my $seq = $inSeqIO->next_seq() ) {
		$total_count++;
		my $id = $seq->display_id;
		$seeds{$id} = $seq;
	}

	print keys(%seeds) . " sequences will be processed...\n";
	&Utils::LogUtils::log_info( $log_fd,
		keys(%seeds) . " sequences will be processed..." );

	my $table1 = ( $nuc_code == 1  ? 0 : $nuc_code );
	my $table2 = ( $mito_code == 1 ? 0 : $mito_code );

# Write 6-frame nuclear protein translations, not sure about the alternative flag (needed for nuclear similarity and longest ORF prediction)
	chdir( $dirs{"process"} );
	print
"Write 6-frame nuclear protein translations (needed for nuclear similarity and longest ORF prediction)...\n";
	&Utils::LogUtils::log_info( $log_fd,
"Write 6-frame nuclear protein translations (needed for nuclear similarity and longest ORF prediction)..."
	);
	system(
"$transeq_exec $seqn_file nucl_6pep.fsa -frame=6 -table=$table1 -alternative 2> /dev/null"
	);

# Write 6-frame mitochondrial protein translations (only if mitochodrial similarity is selected)
	if ($use_mito_sim) {
		print "Write 6-frame mitochondrial protein translations...\n";
		&Utils::LogUtils::log_info( $log_fd,
			"Write 6-frame mitochondrial protein translations..." );
		system(
"$transeq_exec $seqn_file mito_6pep.fsa -frame=6 -table=$table2 -alternative 2> /dev/null"
		);
	}

	# Download bias codon table if prediction with decoder is selected
	if ($use_decoder) {
		print
"Download bias codon table if prediction with decoder is selected...\n";
		&Utils::LogUtils::log_info( $log_fd,
"Download bias codon table if prediction with decoder is selected..."
		);
		chdir( $dirs{"decoder"} );
		system("cp $decoder_exec .")
		  ; # make a copy of the executable in pwd - requirement of decoder fortran code.
		chmod 0755, "decoder";
		&getcodon( $bias_orgasnism, $DECODER_BIASTABLE );
		chdir( $dirs{"process"} );
	}

	# First look for for sequence similarities
	my @seeds;    # blast needs an array of objects NOT a hash.
	for my $key ( keys %seeds ) {
		push @seeds, $seeds{$key};
	}

	$BLASTN_EVAL = &exp2lin($BLASTN_EVAL);
	$BLASTX_EVAL = &exp2lin($BLASTX_EVAL);

	my ( $rRNA_blast, $blast_ref_rRNA );
	if ($use_rrna_sim) {
		print("\nRun blast on ribosomal RNA genes and split reports...\n");
		&Utils::LogUtils::log_info( $log_fd,
			"Run blast on ribosomal RNA genes and split reports..." );
		####################################
	  # Do a system blast rather than a bioperl blast. Can now catch the STDERR
	  # such as "Could not calculate ungapped Karlin-Altschul parameters due..".
	  # Galaxy cannot process output if something is written to STDERR.
		&run_system_blast( $seqn_file, $rrna_db, $blastall_exec, 'blastn',
			$BLASTN_EVAL, 'rna_blast.xml', undef, '1' );
		( $blast_ref_rRNA, undef ) =
		  &split_blast( 'rna_blast.xml', $BLASTX_EVAL, $total_count, 'blastxml',
			\%seeds );
		############ original ################
		#		$rRNA_blast =
		#		  &run_blast( \@seeds, $rrna_db, 'blastn', $BLASTN_EVAL, 'rna.bls',
		#			undef, undef );
		#		( $blast_ref_rRNA, undef, undef ) =
		#		  &split_blast( $rRNA_blast, $BLASTN_EVAL, $total_count, undef,
		#			\%seeds );
		####################################

		# Clean up memory
		undef $rRNA_blast;
	}

	my ( $mito_blast, $blast_ref_mito );
	if ($use_mito_sim) {
		print("Run blast on mitochondrial genes and split reports...\n");
		&Utils::LogUtils::log_info( $log_fd,
			"Run blast on mitochondrial genes and split reports..." );
		####################################
		# Do a system blast rather than a bioperl blast.
		&run_system_blast( $seqn_file, $mito_db, $blastall_exec, 'blastx',
			$BLASTX_EVAL, 'mito_blast.xml', undef, $mito_code );
		( $blast_ref_mito, undef, undef ) =
		  &split_blast( 'mito_blast.xml', $BLASTX_EVAL, $total_count,
			'blastxml', \%seeds );
		########### original ################
		#		$mito_blast =
		#		  &run_blast( \@seeds, $mito_db, 'blastx', $BLASTX_EVAL, 'mito.bls',
		#			undef, $mito_code );
		#		( $blast_ref_mito, undef, undef ) =
		#		  &split_blast( $mito_blast, $BLASTX_EVAL, $total_count, undef,
		#			\%seeds );
		####################################

		# Clean up memory
		undef $mito_blast;
	}

	my ( $blast_ref_nr, $check_ref );
	if ($use_nuc_sim) {
		print("Split results from previous protein blast...\n");
		&Utils::LogUtils::log_info( $log_fd,
			"Split results from previous protein blast (blastx reports)..." );

		if ( $blast_format eq "bls" ) {
			( $blast_ref_nr, $check_ref, undef ) =
			  &split_blast( $blastx_report, $BLASTX_EVAL, $total_count, 'blast',
				\%seeds );
		}
		elsif ( $blast_format eq "xml" ) {
			( $blast_ref_nr, $check_ref, undef ) = &split_blast(
				$blastx_report, $BLASTX_EVAL, $total_count, 'blastxml',
				\%seeds
			);
		}

		# Non fatal issue
		print(
"\nCheck if sequence file has a blast report files (non fatal issue if not)...\n"
		);
		&Utils::LogUtils::log_info( $log_fd,
"Check if sequence file has a blast report files (non fatal issue if not)..."
		);
		foreach my $seed (@seeds) {
			my $id = $seed->display_id;
			unless ( exists $check_ref->{$id} ) {
				print(
"Sequence $id was not found in the blast report file.  Non-fatal.  Process continuing...\n"
				);
				&Utils::LogUtils::log_info( $log_fd,
"Sequence $id was not found in the blast report file.  Non-fatal.  Process continuing..."
				);
				next;
			}
		}
	}

	my $further;

# Hope this part is not a problem. No need to loop if there is no type of blast similarity to be done
#	if ( $use_rrna_sim || $use_mito_sim || $use_nuc_sim ) {
	foreach my $seq (@seeds)
	{    # cycle through each sequence and remove redundancy in the references.
		my $id = $seq->display_id; # not the most efficient way but is clearest.
		if ( exists $$blast_ref_rRNA{$id} ) {
			delete $$blast_ref_mito{$id};
			delete $$blast_ref_nr{$id};
			next;
		}
		elsif ( exists $$blast_ref_mito{$id} ) {
			delete $$blast_ref_nr{$id};
			next;
		}
		elsif ( exists $$blast_ref_nr{$id} ) {
			next;
		}
		else {
			push @$further, $seq;
		}
	}

	#	}
	my $size_rRNA = keys(%$blast_ref_rRNA);
	my $size_mito = keys(%$blast_ref_mito);
	my $size_nr   = keys(%$blast_ref_nr);

	#	my $size_further = scalar(@$further) if $further || 0;
	my $size_further = ( $further ? scalar(@$further) : 0 );
	print
"\nSUMMARY:\tnumber of rRNA: $size_rRNA\n\t\tnumber of mito: $size_mito\n\t\tnumber of nr hits: $size_nr\n\t\tnumber sent for further processing: $size_further\n\n";
	&Utils::LogUtils::log_info( $log_fd,
"SUMMARY:\tnumber of rRNA: $size_rRNA\n\t\tnumber of mito: $size_mito\n\t\tnumber of nr hits: $size_nr\n\t\tnumber sent for further processing: $size_further"
	);

	if ($use_rrna_sim) {
		print("\nCheck rRNA hits and process if there is a hit...\n");
		&Utils::LogUtils::log_info( $log_fd,
			"Check rRNA hits and process if there is a hit..." );

		# Get rRNA blast results and write similarity matches
		my $count_rRNA = 0;
		if ($size_rRNA) {
			open( RNA, ">rRNAhits.txt" )
			  || die "\nERROR: Cannot open rRNAhits.txt\n";
			print "Processing putative rRNA genes...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"Processing putative rRNA genes..." );
			for my $seed_id ( keys %$blast_ref_rRNA ) {
				progress( ++$count_rRNA, $size_rRNA );
				pop @{ $blast_ref_rRNA->{$seed_id}
				  };    # remove the final 'hsp' as it's the seq length
				print RNA ">$seed_id\tputative rRNA\t",
				  $blast_ref_rRNA->{$seed_id}[0]->evalue,
				  "\n";    # display: name, desc, evalue
				print RNA $blast_ref_rRNA->{$seed_id}[0]->query_string,
				  "\n";    # sequence
			}
			close RNA;
		}
		else {
			system ("touch rRNAhits.txt");
			print "No rRNA genes to process.  Proceeding...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"No rRNA genes to process.  Proceeding..." );
		}
	}
	else{
		system ("touch rRNAhits.txt");
	}

	my ( $sixframe_file, @sixframe_seqs, $transln_mito, $transln_nr );
	if ($use_mito_sim) {

		# Get mitochondrial blast results
		print("\nCheck mitochondrial hits and traslate if there is a hit...\n");
		&Utils::LogUtils::log_info( $log_fd,
			"Check mitochondrial hits and if there is a hit..." );
		if ($size_mito) {
			print "\n\nTranslating putative mitochondrial genes...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"Translating putative mitochondrial genes..." );
			$sixframe_file = Bio::SeqIO->new(
				'-file'   => "mito_6pep.fsa",
				'-format' => 'Fasta'
			);
			my $seq_stats_mito;
			while ( my $seq = $sixframe_file->next_seq() ) {
				push @sixframe_seqs, $seq;
			}
			$seq_stats_mito =
			  p4e2::tile_path::tile_path( $blast_ref_mito, \@sixframe_seqs,
				'10', \%seeds );
			print "\nExtending mitochondial protein hits...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"Extending mitochondial protein hits..." );
			$transln_mito =
			  p4e2::tile_path::extend( $seq_stats_mito, \@sixframe_seqs,
				'blast' );
			$transln_mito =
			  p4e2::getDNA::fromBlast( $transln_mito, \%seeds, $mito_code,
				$gen_code_file, \$failed_codon_comp );
			if ( $failed_codon_comp ne "" ) {
				print
"Mitochondrial peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp\n";
				&Utils::LogUtils::log_error( $log_fd,
"Mitochondrial peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp"
				);
				$failed_codon_comp = "";
			}

		}
		else {
			print "No mitochondrial genes to translate.  Proceeding...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"No mitochondrial genes to translate.  Proceeding..." );
		}
	}

	undef @sixframe_seqs;
	$sixframe_file =
	  Bio::SeqIO->new( '-file' => "nucl_6pep.fsa", '-format' => 'Fasta' );
	while ( my $seq = $sixframe_file->next_seq() ) {
		push @sixframe_seqs, $seq;
	}

	if ($use_nuc_sim) {

		# Get precomputed protein blast results
		print(
"\nCheck precomputed protein blast results hits and traslate if there is a hit...\n"
		);
		&Utils::LogUtils::log_info( $log_fd,
"Check precomputed protein blast results hits and traslate if there is a hit..."
		);
		if ($size_nr) {
			print "\n\nTranslating nuclear encoded genes...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"Translating nuclear encoded genes..." );
			my $seq_stats_nr;
			$seq_stats_nr =
			  p4e2::tile_path::tile_path( $blast_ref_nr, \@sixframe_seqs, '10',
				\%seeds );
			print "\nExtending nuclear protein hits...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"Extending nuclear protein hits..." );
			$transln_nr =
			  p4e2::tile_path::extend( $seq_stats_nr, \@sixframe_seqs,
				'blast' );
			$transln_nr =
			  p4e2::getDNA::fromBlast( $transln_nr, \%seeds, $nuc_code,
				$gen_code_file, \$failed_codon_comp );

			if ( $failed_codon_comp ne "" ) {
				print
"Nuclear protein peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp\n";
				&Utils::LogUtils::log_error( $log_fd,
"Nuclear protein peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp"
				);
				$failed_codon_comp = "";
			}
		}
		else {
			print "No nuclear genes to translate.  Proceeding...\n";
			&Utils::LogUtils::log_info( $log_fd,
				"No nuclear genes to translate.  Proceeding..." );
		}
	}

#----------------------------------------------------------------------------------#
	my $count = 1;
	my $size2print;
	$size2print = $size_mito + $size_nr;

# Write peptide sequences from similarity matches
#----------------------------------------------------------------------------------#
# Mithocondrial matches
	if ($use_mito_sim) {
		print
		  "\nWriting translations from mithocondrial similarity matches...\n";
		&Utils::LogUtils::log_info( $log_fd,
			"Writing translations from mithocondrial similarity matches..." );
		&outfiles('open');
		$count = 1;
		my $mito_write = 0;

		for my $id ( keys %$transln_mito ) {

			progress( $count++, $size2print );
			( my $id_p = $id ) =~ s/^(\w\w)C/$1P/;
			my ( $main, $hsps, $prop ) =
			  &get_stats_4_psql( $id_p, $transln_mito->{$id},
				$hsp_ref_counter );

			$main->[0] = $id_p;
			$main->[1] = $id;
			$main->[2] = 'similarity';
			$main->[3] = 'mitochondrially-encoded';
			my $desc =
"putative mitochondrial protein  Method: similarity and extension";
			my $xtn = Bio::SeqIO->new(
				-file   => ">>translations_xtn.fsa",
				-format => 'fasta'
			);
			my $noxtn = Bio::SeqIO->new(
				-file   => ">>translations_noxtn.fsa",
				-format => 'fasta'
			);
			my $coding =
			  Bio::SeqIO->new( -file => ">>nt_coding.fsa", -format => 'fasta' );
			my $seqobj = Bio::Seq->new(
				-display_id => $id_p,
				-desc       => $desc,
				-seq        => $transln_mito->{$id}->[3]
			);
			$xtn->write_seq($seqobj);
			$seqobj->seq( $transln_mito->{$id}->[7] );
			$coding->write_seq($seqobj);
			my $d = $seqobj->desc;
			$d =~ s/\sand\sextension//g;
			$seqobj->desc($d);
			$seqobj->seq( $transln_mito->{$id}->[2] );
			$noxtn->write_seq($seqobj);
			$mito_write++;
			print DB_MAIN join ",", @$main;
			print DB_MAIN "\n";

			foreach my $hsp_stats (@$hsps) {
				print DB_HSPS join ",", @$hsp_stats;
				print DB_HSPS "\n";
				$hsp_ref_counter++;
			}
		}
		&outfiles('close');
	}

#----------------------------------------------------------------------------------#

	my $size = $further ? scalar @$further : 0;

#----------------------------------------------------------------------------------#
# Protein database matches
	if ($use_nuc_sim) {
		print
		  "\nWriting translations from nuclear encoded similarity matches...\n";
		&Utils::LogUtils::log_info( $log_fd,
			"Writing translations from nuclear encoded similarity matches..." );
		&outfiles('open');
		my $nuc_write = 0;
		for my $id ( keys %$transln_nr ) {

			progress( $count++, $size2print );
			( my $id_p = $id ) =~ s/^(\w\w)C/$1P/;

			my ( $main, $hsps, $prop ) =
			  &get_stats_4_psql( $id_p, $transln_nr->{$id}, $hsp_ref_counter );

			my $desc =
"putative nuclear encoded protein  Method: similarity and extension";
			$main->[0] = $id_p;
			$main->[1] = $id;
			$main->[2] = 'similarity';
			$main->[3] = 'nuclear-encoded';

			my $xtn = Bio::SeqIO->new(
				-file   => ">>translations_xtn.fsa",
				-format => 'fasta'
			);
			my $noxtn = Bio::SeqIO->new(
				-file   => ">>translations_noxtn.fsa",
				-format => 'fasta'
			);
			my $coding =
			  Bio::SeqIO->new( -file => ">>nt_coding.fsa", -format => 'fasta' );
			my $seqobj = Bio::Seq->new(
				-display_id => $id_p,
				-desc       => $desc,
				-seq        => $transln_nr->{$id}->[3]
			);
			$xtn->write_seq($seqobj);

			$seqobj->seq( $transln_nr->{$id}->[7] );
			$coding->write_seq($seqobj);
			my $d = $seqobj->desc;
			$d =~ s/\sand\sextension//g;
			$seqobj->desc($d);
			$seqobj->seq( $transln_nr->{$id}->[2] );
			$noxtn->write_seq($seqobj);
			$nuc_write++;
			print DB_MAIN join ",", @$main;
			print DB_MAIN "\n";

			foreach my $hsp_stats (@$hsps) {
				print DB_HSPS join ",", @$hsp_stats;
				print DB_HSPS "\n";
				$hsp_ref_counter++;
			}
		}
		&outfiles('close');
		$size = $further ? scalar @$further : 0;
		print
"\n\n$size sequences about to be passed through ESTScan and DECODER if selected\n";
		&Utils::LogUtils::log_info( $log_fd,
"$size sequences about to be passed through ESTScan and DECODER if selected"
		);
	}

#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
# Sequences with no similarity matches are send to ESTScan for further processing and possible prediction
# For now use already generated matrices
#----------------------------------------------------------------------------------#

	if ( $use_estscan && $estscan_matrix && $size > 0 ) {
		my $estscan_write = 0;
		my $estscan_res;

		print "ESTScan selected: start predicting with ESTScan...\n";
		&Utils::LogUtils::log_info( $log_fd,
			"ESTScan selected: start predicting with ESTScan..." );

		chdir( $dirs{"estscan"} );

		open( TEMP, ">estscan.fsa" )
		  || die "\nERROR: Cannot open estscan.fsa\n";

		foreach my $seqobj (@$further) {
			print TEMP ">", $seqobj->display_id, "\n", $seqobj->seq, "\n";
		}
		close TEMP;

		p4e2::estscan_plugin::run_estscan( $estscan_exec, 'estscan.fsa',
			$estscan_matrix, 'estscan_nt.fsa', 'estscan_pep.fsa', '90' );
		die "Unable to find estscan_nt.fsa\n" if ( !-f 'estscan_nt.fsa' );

		#		system "rm -f estscan.fsa";

		#now have a fasta file with all the estscan predictions.
		#it needs to be parsed similar to the DECoder results.

		my $transln_estscan =
		  p4e2::estscan_plugin::parse_results_mod( $further, 'estscan_nt.fsa',
			'estscan_pep.fsa', \%seeds );
		my $size_estscan = keys %$estscan_res;

		my $count = 1;

		my ($estscan_seq);

 #go through all the seqobjs in further and see which gave 'robust' translations
		my $toDECODER;
	  PARSE: foreach my $seqobj (@$further) {
			my $id = $seqobj->display_id;
			unless ( exists $transln_estscan->{$id} ) {
				push @$toDECODER, $seqobj;
				next PARSE;
			}
		}
		$further      = $toDECODER;
		$size_estscan = keys(%$transln_estscan);

		$transln_estscan =
		  p4e2::getDNA::fromESTScan( $transln_estscan, 'estscan_nt.fsa',
			\%seeds, $nuc_code, $gen_code_file, \$failed_codon_comp );
		if ( $failed_codon_comp ne "" ) {
			print
"ESTScan peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp\n";
			&Utils::LogUtils::log_error( $log_fd,
"ESTScan peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp"
			);
			$failed_codon_comp = "";
		}

		print("\nWriting ESTScan translations...\n");
		&Utils::LogUtils::log_info( $log_fd,
			"Writing ESTScan translations..." );

		chdir( $dirs{"process"} );
		&outfiles('open');

		for my $id ( keys %$transln_estscan ) {
			progress( $count++, $size_estscan );

			( my $id_p = $id ) =~ s/^(\w\w)C/$1P/;
			my $dub = 0;
			$dub = 1
			  if scalar @{ $transln_estscan->{$id} } == 2; #there's two of them.

			foreach my $transln ( @{ $transln_estscan->{$id} } ) {
				$estscan_write++;
				my ( $main, $hsp, $prop ) =
				  &get_stats_4_psql( $id_p, $transln, $hsp_ref_counter );
				my $desc = "putative nuclear encoded protein  Method: ESTScan";
				if ($dub) {
					$id_p =~ s/(_p)?$/_n/ if $main->[6] == '-1';
					$id_p =~ s/(_n)?$/_p/ if $main->[6] == '1';
				}
				$main->[0] = $id_p;
				$main->[1] = $id;
				$main->[2] = 'estscan';
				$main->[3] = 'nuclear-encoded';

				my $xtn = Bio::SeqIO->new(
					-file   => ">>translations_xtn.fsa",
					-format => 'fasta'
				);
				my $noxtn = Bio::SeqIO->new(
					-file   => ">>translations_noxtn.fsa",
					-format => 'fasta'
				);
				my $coding = Bio::SeqIO->new(
					-file   => ">>nt_coding.fsa",
					-format => 'fasta'
				);
				my $seqobj = Bio::Seq->new(
					-display_id => $id_p,
					-desc       => $desc,
					-seq        => $transln->[3]
				);
				$xtn->write_seq($seqobj);

				#change the object to write coding region
				$seqobj->seq( $transln->[7] );
				$coding->write_seq($seqobj);
				my $d = $seqobj->desc;
				$d =~ s/\sand\sextension//g;
				$seqobj->desc($d);
				$seqobj->seq( $transln->[3] );
				$noxtn->write_seq($seqobj);

				print DB_MAIN join ",", @$main;
				print DB_MAIN "\n";

			}
		}
		&outfiles('close');

		#		if ($further) {
		#			$size_further = scalar @$further;
		#		}
		$size_further = ( $further ? scalar(@$further) : 0 );
		print
"\n\nSUMMARY:\tESTScan completed\n\t\tNumber of sequences translated: $size_estscan\n\t\tNumber of sequences sent for further processing: $size_further\n\n";
		&Utils::LogUtils::log_info( $log_fd,
"\n\nSUMMARY:\tESTScan completed\n\t\tNumber of sequences translated: $size_estscan\n\t\tNumber of sequences sent for further processing: $size_further"
		);

	}

#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
# Sequences with no ESTScan matches are send to DECODER for further processing and possible prediction
	if ( $further && $use_decoder && $qualn_file && $bias_orgasnism ) {
		my $dec_write = 0;
		my $pre_DEC   = "";
		print "DECODER selected: start predicting with DECODER...\n";
		&Utils::LogUtils::log_info( $log_fd,
			"DECODER selected: start predicting with DECODER..." );
		my $off2further;
		unless ($pre_DEC) {
			$off2further =
			  &setup_decoder( $dirs{"seqn_qualn"}, $QUALSUFF, $SEQSUFF,
				$DECODER_BIASTABLE, $dirs{"decoder"}, $further );
			&run_decoder( $dirs{"decoder"}, '5',
				$dirs{"decoder"} . "/decoder" );
			&run_decoder( $dirs{"decoder"}, '3',
				$dirs{"decoder"} . "/decoder" );
			$pre_DEC = $dirs{"decoder"};
		}

		my $count = 1;
		my ( $dec_seq, $failed );

		$size_further = @$further;
		foreach my $seqobj (@$further) {
			my $dec_pred;
			my $id = $seqobj->display_id;
			$dec_pred->{$id} = &decoder_parse( $id, $pre_DEC )
			  ;    # check the length of the putative translation

			progress( $count++, $size_further );
			if ( $dec_pred->{$id}->[5] == 0 ) {    # if poor DECoder prediction
				push @$off2further, $seqobj;
				next;
			}
			else {
				my $start = $dec_pred->{$id}->[0];    # this needs tidying up.
				my $end =
				  $dec_pred->{$id}->[1]
				  ;    # messy because of a new location method
				my $frame    = $dec_pred->{$id}->[2];
				my $ntlen    = $dec_pred->{$id}->[3];
				my $pred_seq = $dec_pred->{$id}->[4];
				$dec_seq->{$id}->[0][1] = $start;
				$dec_seq->{$id}->[0][0] = $frame;
				$dec_seq->{$id}->[1][1] = $end;
				$dec_seq->{$id}->[1][0] = $frame;
				$dec_seq->{$id}->[6]    = $ntlen;
				$dec_seq->{$id}->[2]    = $pred_seq;
			}

		}

		my $size_dec = keys(%$dec_seq);

  # below line is an artifact from above statement and will get sorted out ?????
		my $transln_dec = $dec_seq;
		print "\n\nTRANS $transln_dec\n\n";
		$transln_dec =
		  p4e2::getDNA::fromDECODER( $transln_dec, $pre_DEC, $nuc_code,
			$gen_code_file, \$failed_codon_comp );
		if ( $failed_codon_comp ne "" ) {
			print
"DECODER peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp\n";
			&Utils::LogUtils::log_error( $log_fd,
"DECODER peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp"
			);
			$failed_codon_comp = "";
		}

		# Write peptide sequences from DECODER predicitions
		print("\nWriting DECODER translations...\n");
		&Utils::LogUtils::log_info( $log_fd,
			"Writing DECODER translations..." );
		$count = 1;
		chdir( $dirs{"process"} );
		&outfiles('open');

		for my $id ( keys %$transln_dec ) {
			progress( $count++, $size_dec );
			( my $id_p = $id ) =~ s/^(\w\w)C/$1P/;

			my ( $main, $hsp, $prop ) =
			  &get_stats_4_psql( $id_p, $transln_dec->{$id}, $hsp_ref_counter );

			my $desc = "putative nuclear encoded protein  Method: decoder";
			$main->[0] = $id_p;
			$main->[1] = $id;
			$main->[2] = 'decoder';
			$main->[3] = 'nuclear-encoded';

			my $xtn = Bio::SeqIO->new(
				-file   => ">>translations_xtn.fsa",
				-format => 'fasta'
			);
			my $noxtn = Bio::SeqIO->new(
				-file   => ">>translations_noxtn.fsa",
				-format => 'fasta'
			);
			my $coding = Bio::SeqIO->new(
				-file   => ">>nt_coding.fsa",
				-format => 'fasta'
			);
			my $seqobj = Bio::Seq->new(
				-display_id => $id_p,
				-desc       => $desc,
				-seq        => $transln_dec->{$id}->[2]
			);
			$xtn->write_seq($seqobj);

			#change the object to write coding region
			$seqobj->seq( $transln_dec->{$id}->[7] );
			$coding->write_seq($seqobj);
			my $d = $seqobj->desc;
			$d =~ s/\sand\sextension//g;
			$seqobj->seq( $transln_dec->{$id}->[2] );
			$seqobj->desc($d);
			$noxtn->write_seq($seqobj);
			$dec_write++;
			print DB_MAIN join ",", @$main;
			print DB_MAIN "\n";

		}
		&outfiles('close');

		#		$size_further = scalar @$off2further if ( defined @$off2further );
		$size_further = ( $off2further ? scalar(@$off2further) : 0 );
		print
"\n\nSUMMARY:\tDECODER process completed\n\t\tNumber of sequences translated: $size_dec\n\t\tNumber of sequences sent for further processing: $size_further\n\n";
		&Utils::LogUtils::log_info( $log_fd,
"SUMMARY:\tDECODER process completed\n\t\tNumber of sequences translated: $size_dec\n\t\tNumber of sequences sent for further processing: $size_further"
		);
		$further = $off2further;
	}

#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
# If rRNA, mitochondrial and nuclear similarity matches, DECODER or ESTSCAN prediction was unsuccesful/(not selected) the longest Open Reading Frame are found and written
	my $long_write = 0;
	my $j          = 0;
	if ( $size_further > 0 ) {
		print
"$size_further sequences could not be translated using blast similarity matches, ESTScan or DECODER. For these sequences the longest ORF is calculated...\n";
		&Utils::LogUtils::log_info( $log_fd,
"$size_further sequences could not be translated using blast similarity matches, ESTScan or DECODER. For these sequences the longest ORF  is calculated..."
		);
		my $size2sixf = scalar(@$further);
		my $transln_sixf;
		my $count = 1;

#to speed things up need to create a hash with each sequence's sixframe translations.
		my $seqIO = new Bio::SeqIO( -file => 'nucl_6pep.fsa' );
		my %sixframe;
		while ( my $seq = $seqIO->next_seq ) {
			my $id = $seq->display_id;
			$id =~ s/_\d+$//;
			push @{ $sixframe{$id} }, $seq->seq;
		}
		foreach my $seq (@$further) {
			#######################
	  # Hack to discard sequences for which all nucleic acids are flagged with N
			my $seq_str = $seq->seq();
			$seq_str =~ s/[nN]//g;
			unless ( ( length $seq_str ) <= 6 ) { # peptide should be at least 2 amino acids
				#######################
				my $id = $seq->display_id;
				my $seq_fasta =
				  &parse_sixframe( $id, $size2sixf, $count++, $seq->length,
					$sixframe{$id} );
				$transln_sixf->{$id} = $seq_fasta;
			}
			else {
				print $seq->display_id
				  . " cannot be translated, all nucleic acids are flagged with N\n";
				&Utils::LogUtils::log_error( $log_fd,
					$seq->display_id
					  . " cannot be translated, all nucleic acids are flagged with N"
				);
			}
		}
		p4e2::getDNA::fromORF( $transln_sixf, \%seeds, $nuc_code,
			$gen_code_file, \$failed_codon_comp );
		if ( $failed_codon_comp ne "" ) {
			print
"Longest ORF peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp\n";
			&Utils::LogUtils::log_error( $log_fd,
"Longest ORF peptide sequences which failed to match codon table of NCBI standard genetic code\n$failed_codon_comp"
			);
			$failed_codon_comp = "";
		}

		print("\nWriting longest ORF translations...\n");
		&Utils::LogUtils::log_info( $log_fd,
			"Writing longest ORF translations..." );

		$count = 1;
		&outfiles('open');

		for my $id ( keys %$transln_sixf ) {
			progress( $count++, $size2sixf );
			( my $id_p = $id ) =~ s/^(\w\w)C/$1P/;

			my ( $main, undef, $prop ) =
			  &get_stats_4_psql( $id_p, $transln_sixf->{$id} );
			my $desc = "putative nuclear encoded protein  Method: Longest ORF";
			$main->[0] = $id_p;
			$main->[1] = $id;
			$main->[2] = 'longest_orf';
			$main->[3] = 'nuclear-encoded';

			my $xtn = Bio::SeqIO->new(
				-file   => ">>translations_xtn.fsa",
				-format => 'fasta'
			);
			my $noxtn = Bio::SeqIO->new(
				-file   => ">>translations_noxtn.fsa",
				-format => 'fasta'
			);
			my $coding = Bio::SeqIO->new(
				-file   => ">>nt_coding.fsa",
				-format => 'fasta'
			);
			my $seqobj = Bio::Seq->new(
				-display_id => $id_p,
				-desc       => $desc,
				-seq        => $transln_sixf->{$id}->[2]
			);

			$xtn->write_seq($seqobj);

			#change the object to write coding region
			$seqobj->seq( $transln_sixf->{$id}->[7] );
			$coding->write_seq($seqobj);
			my $d = $seqobj->desc;
			$d =~ s/\sand\sextension//g;
			$seqobj->desc($d);
			$seqobj->seq( $transln_sixf->{$id}->[2] );
			$noxtn->write_seq($seqobj);
			$long_write++;
			print DB_MAIN join ",", @$main;
			print DB_MAIN "\n";
		}
		&outfiles('close');
	}

	print "\n\nAll sequences have now been processed.\n\n";
	&Utils::LogUtils::log_info( $log_fd,
		"All sequences have now been processed.." );

	my $predicted_seqp_file = $dirs{"out"} . "/seqp.fsa";
	my $predicted_seqr_file = $dirs{"out"} . "/seqr.fsa";

	# Prepare peptide output file output files
	system( "cp "
		  . $dirs{"process"}
		  . "/translations_xtn.fsa $predicted_seqp_file" );
	system("cp  $predicted_seqp_file $seqp_file");
		system( "cp "
		  . $dirs{"process"}
		  . "/rRNAhits.txt $predicted_seqr_file" );
	system("cp $predicted_seqr_file $seqr_file");
	close $log_fd;
}
#############################  End main() ##########################################

########### Get codon table of a specific organism e.g. Oryza sativa ###############
###### For now choose largest codon table from the related organism ################
################ Exits if no codon table is found for the organism #################
sub getcodon {

	my $spe_input   = shift;
	my $output_file = shift;

	print(
"Fetching a codon bias table for $spe_input.\nthis may take a few seconds\n"
	);

	my ( $spe_query, $return_data, $codon_table, $spe, $line1, $ch_add,
		$thisURL );
	my $query_return = "qret";
	my $file_return  = "fret";
	my $file_count   = 1;

	my $max_nr_codons = 0;
	my $codons_url    = "";
	my $codons_spe    = "";

	( $spe_query = $spe_input ) =~ s/ /\+/g;    # insert + for url

	print("Searching, please wait...\n");
	system(
"wget -w 10 -q -Oqret 'http://www.kazusa.or.jp/codon/cgi-bin/spsearch.cgi?species=$spe_query&c=i'"
	);
	print(
"wget -w 10 -q -Oqret 'http://www.kazusa.or.jp/codon/cgi-bin/spsearch.cgi?species=$spe_query&c=i'\n"
	);

	open( Q_RET, "<$query_return" );    # open the wget returned web page file

	while ( my $line = <Q_RET> ) {
		if ( $line =~
			/<A HREF\="\/codon\/cgi-bin\/showcodon.cgi\?species=(.*)">/ )
		{

			system(
"wget -qOcret$file_count 'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=$1&aa=1&style=GCG'"
			);
			print(
"wget -qOcret$file_count 'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=$1&aa=1&style=GCG'\n"
			);
			$thisURL = $1;
			open( T_RET, "<cret$file_count" )
			  ;    # for each species open the wget returned web page file
			while ( my $line1 = <T_RET> ) {

				if ( $line1 =~ /<STRONG>Not found/ )
				{    # if there is no codon table
					$thisURL = '';
					next;
				}
				elsif ( $line1 =~ /<STRONG><i>(.*)<\/STRONG>/ )
				{    # if codon table present, parse + keep the details
					$ch_add = $1;
					$ch_add =~ s/<\/i>.*:/-/;
					$ch_add =~ s/<\/STRONG>//;
					if ($thisURL) {
						if ( $ch_add =~ /\((\d+)\scodons\)/ ) {
							my $codons = $1;
							if ( $codons > $max_nr_codons ) {
								$max_nr_codons = $codons;
								$codons_url    = $thisURL;
								$codons_spe    = $ch_add;
							}
						}
					}
				}
			}
			$file_count++;
		}
	}

	if ( $max_nr_codons == 0 ) {
		print("0 matches for $spe_input.\nIs the organism name correct?\n");
		print(
"Choose another, closely related organism.  For full list check out http://www.kazusa.or.jp/codon/\n\n"
		);
		exit;
	}

	system(
"wget -qOfret 'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=$codons_url&aa=1&style=GCG'"
	);
	print(
"wget -qOfret 'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=$codons_url&aa=1&style=GCG'\n"
	);

	open( F_RET, "<$file_return" );
	while ( my $line = <F_RET> ) {    # extract codon table from web page
		$return_data .= $line;
	}

	if ( $return_data =~ /Fraction.*?\.\.\n(.*)<\/PRE>/s ) {
		$codon_table = $1;
	}

	open( OUTPUT, ">$output_file" );
	print OUTPUT $codon_table;
	close OUTPUT;
	print(
		"\nThe codon table for $codons_spe has been saved to $output_file\n\n\n"
	);
	system("rm cret* fret qret");
	return 1;
}
######################### end getcodon() ###########################################

######################### Run blast ################################################
sub run_blast {
	my ( $seqs2blast, $dbase_path, $blast_method, $e, $o, $f, $gcode ) = @_;
	$f = 'T' if ( !$f );
	my @params;

	if ($o) {
		@params = (
			'd'           => "$dbase_path",
			'program'     => "$blast_method",
			'e'           => "$e",
			'o'           => "$o",
			'b'           => '25',
			'v'           => '25',
			'_READMETHOD' => "Blast",
			'F'           => $f,
			'Q'           => $gcode,
		);
	}
	else {
		@params = (
			'd'           => "$dbase_path",
			'program'     => "$blast_method",
			'e'           => "$e",
			'b'           => '25',
			'v'           => '25',
			'_READMETHOD' => "Blast",
			'F'           => $f,
			'Q'           => $gcode,
		);
	}

	my $factory      = Bio::Tools::Run::StandAloneBlast->new(@params);
	my $blast_report = $factory->blastall($seqs2blast);

	return ($blast_report);
}
######################### end run_blast() ###########################################

######################### Run system blast ################################################
sub run_system_blast {
	my ( $blast_file, $dbase_path, $blast_exec, $blast_method, $e, $o, $f,
		$gcode )
	  = @_;
	$f = 'T' if ( !$f );
	system(
"$blast_exec -p $blast_method -d $dbase_path -i $blast_file -v 25 -b 25 -m 7 -e $e -F f -Q $gcode -o $o 2> /dev/null"
	);
}
######################### end run_blast() ###########################################

######################### Split blast file ##########################################
sub split_blast {
	my $in           = shift;
	my $evalue       = shift;
	my $total_count  = shift;
	my $blast_format = shift;
	my $seeds =
	  shift;   # this is a ref to a has which holds the sequence for querying.
	           # its only used when parsing the blast results against swissprot.
	           # if the user selects a previous set of results then they may be
	           # entries for queries which are not in the query set.
	    # this speeds things up greatly and makes downstream events more robust.
	my $check_ref;
	my %hsp_box;
	my $tophit_name;
	my $current_count = 1;

	my $type = ref($in);
	my $searchio;
	my %search_res;
	my $stop;
	my $num_of_seeds = keys %$seeds;
	if ( $type =~ m/blast/i ) {    # the input is already an object
		$searchio = $in;
	}
	else {
		$searchio = new Bio::SearchIO(
			-format => $blast_format,
			-file   => $in
		);
	}

  BLAST: while ( my $results = $searchio->next_result ) {
		my $query_name = $results->query_name;
		$query_name =~ s/,.*//;    # added for hgmp output
		$query_name =~ s/;//;      # stray semi-colons need removing
		my $hsp_count;
		next BLAST unless ( exists $seeds->{$query_name} );
		next BLAST if ( exists $hsp_box{$query_name} );
		if ($total_count) {
			progress( $current_count++, $num_of_seeds );
		}
		$check_ref->{$query_name} = '1';

		my @hsps;
		my $count = 0;
	  HIT: while ( my $hit = $results->next_hit ) {
			if ( $hit->significance < $evalue ) {
				if ( !$count ) {
					$tophit_name = $hit->name;
					$count       = 1
					  ; # have top hit now don't come back in here 'til next result.
				}
				while ( my $hsp = $hit->next_hsp ) {
					if ( $hsp->evalue < $evalue ) {
						push @{ $hsp_box{$query_name} }, $hsp;
						$hsp_count++;
						if ( $hsp_count > 5 ) {
							push @{ $hsp_box{$query_name} },
							  $results->query_length;
							next BLAST;
						}
					}
					elsif ( exists $hsp_box{$query_name}
						&& scalar @{ $hsp_box{$query_name} } == 0 )
					{  # this is similar to $res->num_hits with an Evalue cutoff
						    # different from that of the report used.
						    # If no hit was under Evalue then there will not be
						    # an entry in %hsp_box
						push @{ $hsp_box{$query_name} }, $results->query_length;
						next BLAST;
					}
				}
			}
			else {
				next HIT;
			}
		}

		push @{ $hsp_box{$query_name} }, $results->query_length
		  if ( $hsp_box{$query_name} );
		my $num = 0;
		if ( exists $hsp_box{$query_name} ) {
			foreach my $hsp ( @{ $hsp_box{$query_name} } ) {
				if ( $hsp =~ m/^\d+$/ ) {
					$num++;
				}
				if ( $num > 1 ) {
					print
"Error is parsing the blast lengths:\n$query_name\t$num\t$hsp\n";
					die;
				}
			}
		}
	}
	undef %search_res;
	undef $searchio;
	return ( \%hsp_box, $check_ref, $tophit_name );
}
######################### end split_blast() ########################################

######################### progress() ###############################################
sub progress {
	my $current    = shift;
	my $total      = shift;
	my $percentage = ( $current / $total ) * 100;
	my $progress   = sprintf "%5.2f", $percentage;
	print("\r$progress%");
}
######################### end progress() ###########################################

######################### outfiles() ###############################################
sub outfiles {
	my $action = shift;

	if ( $action eq 'open' ) {
		open DB_MAIN, ">>prot_main.psql" or die;
		open DB_HSPS, ">>prot_hsps.psql" or die;
		open DB_PROP, ">>prot_prop.psql" or die;
	}
	elsif ( $action eq 'close' ) {
		close DB_MAIN;
		close DB_HSPS;
		close DB_PROP;
	}
	else {
		print "Cannot determine file action: $action!\n";
		die;
	}
}
######################### end outfiles() ############################################

######################### exp2lin() #################################################
sub exp2lin {
	my $e = shift;
	if ( $e =~ m/((?:\d+\.?\d*|\.\d+)[eE][+-]?\d+)/ ) {
		my $exp = $1;
		my $limit;
		( $limit = $e ) =~ s/(?:\d+\.?\d*|\.\d+)[eE][+-]?(\d+)/$1/
		  ;          # need to do this as issues if limit set too high
		$limit++;    # compared to that following e/E.
		$limit .= 'f';
		$_ = sprintf "%.$limit", $exp;
		s/0+$//;
		s/\.$//;
		return $_;
	}
	elsif ( $e =~ m/^\s*\d+\.?\d+\s*$/ ) {
		return $e;    # do nothing to e_rRNA as already linear
	}
	else {
		return 0;     # error - can't recognise the number
	}
}
######################### end exp2lin() #############################################

######################### Reverse and complement dna string #########################
sub revDNA {

	my ( $file, $flag, $title1, $seq, $line, $title );
	opendir( DIR, "./Translate.5/" );

	while ( $file = readdir(DIR) ) {
		if ( $file =~ /seq/ ) {
			open( INPUT, "<./Translate.5/$file" )
			  || die "Unable to open Input\n";
			$flag   = 0;
			$title1 = '';
			while ( $line = <INPUT> ) {
				if ( $line =~ /^(>.+)/ ) {
					$title = $1;
					if ( $title1 =~ /./ ) { gorev( $seq, $title1, $file ); }
					$seq = '';
					next;
				}
				else {
					chomp $line;
					$seq .= $line;
					$title1 = $title;
				}
			}
			gorev( $seq, $title1, $file );
			close(INPUT);
		}
		if ( $file =~ /qlt/ ) {
			open( INPUT, "<./Translate.5/$file" )
			  || die "Unable to open Input\n";
			$flag   = 0;
			$title1 = '';
			while ( $line = <INPUT> ) {
				if ( $line =~ /^(>.+)/ ) {
					$title = $1;
					if ( $title1 =~ /./ ) {
						revqual( $seq, $title1, $file );
					}
					$seq = '';
					next;
				}
				else {
					chomp $line;
					$seq .= $line;
					$title1 = $title;
				}
			}
			revqual( $seq, $title1, $file );
			close(INPUT);
		}
	}
}
######################### end revDNA() ##############################################

######################### Reverse dna string ########################################
sub gorev() {
	my $seq    = shift;
	my $title1 = shift;
	my $file   = shift;
	my @tmpseq = split( '', $seq );
	my $l      = length($seq);
	my @revcomp;
	my $i;
	for ( $i = 0 ; $i < $l ; $i++ ) {
		if    ( $tmpseq[$i] eq ( 'A' || 'a' ) ) { $revcomp[$i] = 'T'; }
		elsif ( $tmpseq[$i] eq ( 'C' || 'c' ) ) { $revcomp[$i] = 'G'; }
		elsif ( $tmpseq[$i] eq ( 'G' || 'g' ) ) { $revcomp[$i] = 'C'; }
		elsif ( $tmpseq[$i] eq ( 'T' || 't' ) ) { $revcomp[$i] = 'A'; }
		elsif ( $tmpseq[$i] eq 'N' ) { $revcomp[$i] = 'N'; }
		elsif ( $tmpseq[$i] eq 'X' ) { $revcomp[$i] = 'X'; }
		elsif ( $tmpseq[$i] ne ' ' ) { $revcomp[$i] = $tmpseq[$i]; }
	}
	open( OUTFILE, ">>./Translate.3/$file" )
	  || die "Can't open $file.rev\n";
	print OUTFILE "$title1\n";
	my $n = 0;
	for ( $i = $l - 1 ; $i > -1 ; $i-- ) {
		print OUTFILE "$revcomp[$i]";
		$n++;
		if ( ( $n % 50 ) == 0 ) {
			print OUTFILE "\n";
		}
	}
	if ( ( $n % 50 ) != 0 ) { print OUTFILE "\n"; }
	close OUTFILE;
}
######################### end gorev() ###############################################

######################### Reverse quality scores ####################################
sub revqual() {
	my $seq    = shift;
	my $title1 = shift;
	my $file   = shift;
	my @tmpseq = split( ' ', $seq );
	my $l      = scalar(@tmpseq);
	open( OUTFILE, ">>./Translate.3/$file" )
	  || die "Can't open $file.rev\n";
	print OUTFILE "$title1\n";
	my $n = 0;
	my $i;

	for ( $i = $l - 1 ; $i > -1 ; $i-- ) {
		print OUTFILE "$tmpseq[$i] ";
		$n++;
		if ( ( $n % 50 ) == 0 ) {
			print OUTFILE "\n";
		}
	}
	if ( ( $n % 50 ) != 0 ) { print OUTFILE "\n"; }
	close OUTFILE;
}
######################### end revqual() #############################################

######################### Decoder setup #############################################
sub setup_decoder {
	print "Setting up DECODER\nThis may take a few minutes\n";

	my $path2files = shift;
	my $qualsuff   = shift;
	my $seqsuff    = shift;
	my $cut        = shift;
	my $dir        = shift;
	my $seqs       = shift;

	my $T    = "$dir/Translate";
	my $T5   = "$dir/Translate.5";
	my $T3   = "$dir/Translate.3";
	my $P5   = "$dir/Peptide.5";
	my $P3   = "$dir/Peptide.3";
	my @dirs = ( $T, $T5, $T3, $P5, $P3 );
	my @notThere;
	foreach my $dir (@dirs) { mkdir $dir, 0755; }

	for ( my $i = 0 ; $i < @$seqs ; $i++ ) {
		my $id = $seqs->[$i]->display_id;
		$id =~ s/\.Contig/_/;
		if ( -e "$path2files/$id.$qualsuff" ) {
			copy( "$path2files/$id.$qualsuff", $T5 ) or die;
			copy( "$path2files/$id.$seqsuff",  $T5 ) or die;

		}
		else {
			push @notThere, $seqs->[$i];
			splice @$seqs, $i, 1;    # remove id from the array
		}
	}
	chdir $T5;

	&FILENAME_CHANGE( $qualsuff, 'qlt' );
	&FILENAME_CHANGE( $seqsuff,  'seq' );
	chdir '../';                     # move into decoder directory
	&revDNA;
	chdir '../';                     # move into outdir
	return \@notThere;
}
######################### end setup_decoder() #######################################

######################### Run DECODER ###############################################
sub run_decoder {
	my $dir          = shift;
	my $direction    = shift;
	my $path2decoder = shift;
	chdir $dir or die;
	my @totalfiles    = glob("./Translate.$direction/*.seq");
	my $total_count   = scalar(@totalfiles);
	my $current_count = 0;
	opendir( DIR, "./Translate.$direction" )
	  or die;    #o pen Translate.5 or .3 Dir
	my $file;
	print "\nTranslating $direction' sequences with DECODER\n";

	while ( defined( $file = readdir(DIR) ) )
	{            # while there are files in the Dir
		if ( $file =~ /\.seq/ ) {    # if it's a seq file
			system("cp Translate.$direction/$file input.seq")
			  ;                      # copy the file to input.seq
			$file =~ s/seq/qlt/;     # change from .seq to .qlt
			system("cp Translate.$direction/$file input.qlt")
			  ;                      # copy the file to input.qlt

# system("$path2decoder >& /dev/null");		# run decoder to make result.pep & .rqt & .tfa
			system("$path2decoder 2> /dev/null")
			  ;    # run decoder to make result.pep & .rqt & .tfa

			$file =~ s/qlt/pep/;    # change from .qlt to .pep
			system("mv result.pep ./Peptide.$direction/$file")
			  ;                     # move result.pep to Peptide dir and rename
			$file =~ s/pep/rqt/;    # change from .pep to .rqt
			system("mv result.qlt ./Peptide.$direction/$file")
			  ;                     # move result.rqt to Peptide dir and rename
			$file =~ s/rqt/tfa/;    # change from .rqt to .tfa
			system("mv result.tfa Peptide.$direction/$file")
			  ;                     # move result.tfa to Peptide dir and rename
			$current_count++;
			my $percentage = ( $current_count / $total_count ) * 100;
			my $progress = sprintf "%5.2f", $percentage;
			print "\r$progress%";
		}
	}
	closedir(DIR);                  # close directory
	chdir '../';
}
########################## end run_decoder() ########################################

######################### Parse DECODER results #####################################
sub decoder_parse {
	my ( $input_id, $path ) = @_;
	my ( $off2ESTScan, $decoderfile5, $decoderfile3, $big, $decoder_pred,
		$strand );

	#	print LOG "<DECODER_PARSE>\n";

	$input_id =~ s/\.Contig/_/g;
	$path =~ s/\/$//;    # remove end / if in path.
	if ( $path =~ m/prot4estdecoderfiles/ )
	{    # uses the results created by this process
		$path =~ s/prot4est//;    # remove identifier
		$decoderfile5 = "./$path/Peptide.5/$input_id.pep";
		$decoderfile3 = "./$path/Peptide.3/$input_id.pep";
	}
	else {                        # use the absolte path given in config file
		$decoderfile5 = "$path/Peptide.5/$input_id.pep";
		$decoderfile3 = "$path/Peptide.3/$input_id.pep";
	}

	# both of these could be the correct translation!
	# Pick longest.
	if ( -s "$decoderfile3" >
		-s "$decoderfile5" )  # if file in .3 dir is greater than file in .5 dir
	{
		$big    = $decoderfile3;
		$strand = '-1';
	}
	else {
		$big    = $decoderfile5;
		$strand = '1';
	}

	if ( -f $big ) {
		my $in = Bio::SeqIO->new( -file => "$big", -format => 'Fasta' );
		while ( my $seq = $in->next_seq() ) {
			my $acc  = $seq->display_id;
			my $desc = $seq->desc();
			( my $seqstring = $seq->seq ) =~ s/X+$//g;

			my $seqlen = length($seqstring);

			if ( $seqlen == 0 ) {
				$decoder_pred->[5] = '0';
				next;
			}
			$seq->seq($seqstring);

			$decoder_pred->[3] = $seqlen;
			my ( $bp, $aa );

			if ( $acc && $desc ) {
				$acc =~ s/^(\w+_?\d?).*/$1/i;
				$seq->display_id($acc);
			}
			else {
				die "\n\nERROR: acc=$acc\tdesc=$desc\tlength=$seqlen\n";
			}

			if ( $desc =~ m/\w+\s+(\d+)\sbp\s+(\d+)\saa.+/ ) {
				$bp = $1;
				$aa = $2;
			}
			else {
				print
"\n\n### ERROR ###\n$acc - problem with header in DECODER proediction\n";
				exit(1);
			}

			# determine start and end of peptide w.r.t. the EST
			if ( $desc =~ m/cds\s+(\d+)\s+-\s+(\d+)/ ) {
				my %mod;
				@mod{ 0, 1, 2 } = ( 3, 1, 2 );    # remainder maps to frame.
				my $start = $1;
				my $end   = $2;
				$decoder_pred->[0] = $start;
				$decoder_pred->[1] = $end;
				my $rem;
				if ( $strand eq '1' ) {

			#positive strand.  So modulus so start by three will indicate frame.
					$rem = $start % 3;
				}
				else {

	 #minus strand so need to use diff between start and end to determine frame.
					$rem = ( $seqlen - $start + 1 ) % 3;
				}
				$decoder_pred->[2] = $strand * $mod{$rem};
			}

			my $max_trans =
			  sprintf( "%.0f", $bp / 3 );    #convert length of contig to aa

			my $max_10pc = sprintf( "%.0f", $max_trans * 0.1 )
			  ;    #10% of $max_trans.  Dear user feel free to alter this.
			       #Let me know how you get on.
			if ( $seqlen < 30 ) {    #too short!
				$decoder_pred->[5] = '0';

				#$decoder_pred->[4]=$seq->seq;
			}
			elsif ( $seqlen < $max_10pc ) {    #too short relative to contig
				$decoder_pred->[5] = '0';

				#$decoder_pred->[4]=$seq->seq;
			}
			else {                             #sequence is okay!?!
				$decoder_pred->[5] = '1';
				$decoder_pred->[4] = $seq->seq;
			}
		}
	}
	else {                                     #Non Existant
		$decoder_pred->[5] = '0';
		$decoder_pred->[4] = undef;
	}

	#	print LOG "Leaving <DECODER_PARSE>\n";
	return $decoder_pred;
}
######################### end decoder_parse() #######################################

######################### get_stats_4_psql() ########################################
sub get_stats_4_psql {
	my $id        = shift;
	my $seq_stats = shift;
	my $counter   = shift;

	my ( $xtn, $xtn_start, $xtn_end, $note, $hsps );
	my ( @main, @hsps, @prop );
	$main[5] = $seq_stats->[0][1];    #conf start
	$main[6] = $seq_stats->[0][0];    #start frame
	$main[7] = $seq_stats->[1][1];    #conf end
	$main[9] = $seq_stats->[1][0];    #end frame

	if ( $seq_stats->[0][2] ) {
		$main[4] = $seq_stats->[0][2];    #extn start
	}
	else {
		$main[4] = '';
	}
	if ( $seq_stats->[1][2] ) {
		$main[8] = $seq_stats->[1][2];    #extn end
	}
	else {
		$main[8] = '';
	}

	#remove for release.
	foreach ( $main[4], $main[5], $main[7], $main[8] ) {
		if ( $_ =~ m/-/ ) {
			print "$id caused a crash: $main[4],$main[5],$main[7],$main[8]\n";
		}
	}

	my $count = 0;
	foreach my $hsp ( @{ $seq_stats->[5] } ) {
		my @indiv_hsp;
		$indiv_hsp[0] = $counter++;
		$indiv_hsp[1] = $id;
		for ( my $i = 0 ; $i < 5 ; $i++ ) {
			$indiv_hsp[ $i + 2 ] = $hsp->[$i];
		}
		push @hsps, \@indiv_hsp;
	}

	return ( \@main, \@hsps, \@prop );
}
######################### end get_stats_4_psql() ####################################

######################### Change filename extenstion ################################
sub FILENAME_CHANGE {
	my $old_symbol = $_[0];       #part of file name to change ie contigs.qual
	my $new_symbol = $_[1];       #replacement string
	my @dir        = glob("*");
	my $count      = 0;
	foreach my $file (@dir) {
		print "\r", ++$count;

	#Substitute in names of blast files
	#if the new string is empty or the current file name doesn't contain the new
	#string and doesn't start with a stop
		if (   ( $new_symbol eq '' || $file !~ /$new_symbol/ )
			&& $file !~ /^\./
			&& $file =~ m/$old_symbol/ )
		{
			my $new_file = $file;    #define $new_file as old file name
			$new_file =~
			  s/$old_symbol/$new_symbol/g;    #substitute old string with new
			system "mv '$file' '$new_file'";  #rename old file
		}
	}
}
######################### end FILENAME_CHANGE() #####################################

######################### parse_sixframe() ##########################################
sub parse_sixframe {

	my $seqid        = shift;
	my $size         = shift;
	my $count        = shift;
	my $seed_seq_len = shift;
	my $seqs = shift;    #an aref to all sixframes of the sequence in question

	#	print LOG "<PARSE_SIXFRAME>\n";
	my $i       = 0;
	my $frame   = 0;
	my $longest = 0;
	my ( $start, $end, $longest_frame, $aa_start, $aa_end );
	my @allframes;
	my $longest_seq;

	#	my $infile = "nucl_6pep.fsa";
	#	my $seqIO = new Bio::SeqIO( -file => $infile );

	foreach my $seq (@$seqs) {
		$frame++;
		while ( $seq =~ /\**(\w+)\**?/gs )
		{    #this allows nnnnnnn* , *nnnnnn , *nnnnnn*
			my $current = $1;
			$current =~ s/X+$//i;
			$current =~ s/^X+//i;
			my $seq_len = length($current);
			if ( $seq_len > $longest ) {
				$longest       = $seq_len;
				$longest_seq   = $current;
				$longest_frame = $frame;
			}
		}
	}

	#now need to find its location w.r.t to the EST
	my $est_seq = $seqs->[ $longest_frame - 1 ];

	my $sixf_stats;
	if ( $est_seq =~ m/\Q$longest_seq\E/ ) {
		$aa_start = $-[0] + 1;
		$aa_end   = $+[0] - 1;
	}
	else {
		print "###Error###\nCannot determine location of putative peptide\n";
		print "$longest_frame\n$est_seq\n$longest_seq\n";
		print ":";
		print join "\n:", @$seqs;
		print "\n";
		exit;
	}

	if ( $longest_frame > 3 ) {    #minus strand
		my $frame  = $longest_frame - 3;
		my $offset = ( $aa_start * 3 ) - ( 7 - $longest_frame );
		my $l1     = $seed_seq_len;
		my $nt_end =
		  $l1 -
		  $offset
		  ; #second term is the distance between end of EST and the start of peptide.
		my $l2 = length($longest_seq) * 3;
		my $nt_start = ( ( $nt_end - $l2 ) + 1 );
		if ( $nt_start < 1 ) {
			if ( $longest_seq =~ m/[LVSPTARG]/ )
			{    #examples where the 3rd position in codon is redundant
				$nt_start = 1;    #so transeq uses a two ntide codon.
			}
		}
		die
"DIED $nt_start, $longest_frame, $aa_start, $l1, $l2, $nt_end, $seqid,\n$longest_seq\n"
		  if $nt_start < 1;
		die
		  if $nt_end > $l1;
		$sixf_stats->[0][1] = $nt_start;
		$sixf_stats->[1][1] = $nt_end;
		$longest_frame      = ( $longest_frame - 3 ) * -1;
	}
	else {
		my $offset   = $longest_frame - 1;
		my $nt_start = $aa_start * 3 - 2 + $offset;
		$sixf_stats->[0][1] = $nt_start;
		$sixf_stats->[1][1] = $nt_start + 3 * length($longest_seq) - 1;
	}
	$sixf_stats->[0][0] = $longest_frame;
	$sixf_stats->[1][0] = $longest_frame;
	$sixf_stats->[2]    = $longest_seq;

	my $start_met = index( $longest_seq, "M" ) + 1;
	my $desc;
	if ( $start_met > 0 ) {
		$desc = "Possible METstart at position: $start_met";
	}

	#	print LOG "Leaving <PARSE_SIXFRAME>\n";

	progress( $count, $size );

	return $sixf_stats;
}
######################### end parse_sixframe() ######################################

