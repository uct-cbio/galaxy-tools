#!/usr/bin/perl

package p4e2::estscan_plugin;

$VERSION = 1.4;

#10/08/04 - Have moved some of the subs in here from main code. Also collects and builds pseudo CDS
#12/10/04 - SRS @ ebi changed download procedure - now compatible
#04/01/05 - now allows two "viable" translations from ESTScan per input sequence (one from each strand)
#06/06/05 - fixes bug: copying a a second matrix did not work. Thanks to Cyril Cabane

use strict;
use Bio::SeqIO;
use Bio::Species;
use File::Copy;
use Term::ANSIColor;
use Cwd;

sub get_embl_seqs {
	my $species = shift;
	my $dir     = shift;
	my $flag    = 1;
	my $id;
	my $out;
	my $currentdir = getcwd();

  FLAG: while ($flag) {
		( my $query = $species ) =~
		  s/([^\w()'*~!.-])/sprintf('%%%02x', ord($1))/ge;

		#makes the species name suitable for URL.
		my $lines = '';
		$species =~ s/^(\w)/\u$1/;    #capitalise the first character in genus
		if ( ( $id = $species ) =~ s/^\s*(\w)\w+\s+(\w).*/$1$2/ ) {
			$out = "$id" . "_embl.db";
		}

		else {
			print "### WARNING ###";
			print "\n$species does not appear to be in a parsable format.\n";
			print
"To use ESTScan then please enter an organism name for matrix building OR type 'n' to contiue without ESTScan preditions\n";

			chomp( my $answer = <STDIN> );
			if ( $answer =~ m/^\s*[n,N]\s*$/ ) {
				return 0;
			}
			else { $species = $answer; next FLAG; }
		}
		my $desc1 = 'complete';
		my $desc2 = 'cds';
		my $desc3 = 'RNA | mRNA';

		#my $desc4='mRNA';

		#build query for srs search
		my $srs1 =
"wget http://srs.ebi.ac.uk/srs7bin/cgi-bin/wgetz'?-e+-vn+2+[embl-Organism:";
		my $srs2 = "]%20&%20[embl-Description:";
		my $srs3 = "]%20&%20[embl-Description:";
		my $srs4 = "]%20&%20[embl-Molecule:";
		my $srs5 = "]' -O $out";

		my $srs = join( "",
			$srs1,  $query, $srs2,  $desc1, $srs3,
			$desc2, $srs4,  $desc3, $srs5 );

#my $srs = "wget http://srs.ebi.ac.uk/srs7bin/cgi-bin/wgetz?-e+%5Bembl-Organism:Lumbricus%20rubellus%5%20&%20%5Bembl-Description:complete%5%20&%20%5Bembl-Description:cds%5%20&%20%5Bembl-Molecule:RNA%20%7C%20mRNA%5";
		system "$srs";

		my $count = 0;
		open EMBL, "$out" or die;
		open EMBL_EDIT, ">tmp.embl";
		while ( my $line = <EMBL> ) {
			$lines .= $line;
		}
		my $check;
		my $nohits =
"<TD VALIGN=\"center\"><span class=\"errmsg1\">Error: No entries found<\/span><\/TD>"
		  ;    #need to be explicit to speed up search over very large files
		while ( ( $lines =~ m/$nohits/io ) && ( $check == 0 ) )
		{      #none were found.
			print
"Unable to find any complete CDS entries for $species in the EMBL database!\n"
			  . "Complete CDS are required to construct a matrix which is used by the hidden Markov models in ESTScan.\n"
			  . "If you wish to exploit ESTScan then please enter another organism name - one that is closely related is optimal.\n"
			  . "This will mean more peptides are predicted based upon longest ORF from six-frame translation.\n"
			  . "Please another choose species OR type 'n' to cancel connection to EMBL database.\n";
			my $answer = <STDIN>;
			chomp $answer;
			if ( $answer =~ m/^\s*[nN]\s*$/ ) {
				unlink $out;
				return
				  0;  #this will fail condition to enter ESTScan part of script.
			}
			elsif ( $answer =~ m/(\w+)/ ) {
				$species = $answer;
				$flag    = 1;
				$check   = 1;
				redo FLAG;
			}
			else { $check = 0; }    #go back an ask the question!
		}
		if ( $lines =~ m/<pre>/ ) {    #file contains embl entries
			$flag = 0;
		}
		else { die; }
		$lines =~ s/\<.*?\>//g;                #remove html formatting
		$lines =~ s/(ID.*;?) mRNA;/$1 RNA;/g
		  ;    #editing because of small bug in buildmodel script
		print EMBL_EDIT "$lines";
	}
	close EMBL_EDIT;
	close EMBL;
	move( 'tmp.embl', $out );

	my $matrix;
	( $matrix, $species ) = &build_matrix( $species, $id, $out );

	my @ret = ( $matrix, $out, $species );
	return @ret;
}
################################################################################################################
sub build_matrix {
	my $species  = shift;
	my $id       = shift;
	my $seq_file = shift;

	my $matrix = $id . ".smat";
	my $conf   = $id . "_embl.conf";
	open CONF, ">$conf";
	print CONF
"################################################################################\n"
	  . "#\n#   Parameters for $species, use EMBL data\n#   (use PERL syntax!)\n#\n\n"
	  . "\$organism       = \"$species\";\n"
	  . "\$dbfiles        = \"$seq_file\";\n"
	  . "\$ugdata         = \"\";\n"
	  . "\$estdata        = \"\";\n\n"
	  . "\$datadir        = \".\";\n"
	  . "\$nb_isochores   = 1;\n"
	  . "\$tuplesize      = 6;\n"
	  . "\$minmask        = 30;\n\n"
	  . "#\n#   End of File\n#\n"
	  . "################################################################################";

	system "build_model $conf";
	my @matrixold = glob('Matrices/*.smat');
	move( "$matrixold[0]", "./$matrix" );
	die if ( -f $matrixold[0] );
	my @report = glob('Report/*_process_data.log');
	unless ( $report[0] ) {
		print
"###error: cannot find the process_data log file!\nPlease email me at nematode.bioinf\@ed.ac.uk\n\n";
		die;
	}
	my $num_coding = &count_codingNt( $report[0] );
	my $limit      = 150000;
	my $check      = 0;
	while ( ( $num_coding < $limit ) && ( $check == 0 ) ) {
		print "### WARNING ###\n";
		print "Have found only $num_coding coding nucleotides for $species\n";
		print "We suggest that at least $limit coding nucleotides are used.\n";
		print "Using fewer *may* compromise construction of ESTScan matrix.\n";
		print
"These may be complemented with pseudo entries built by prot4EST from the similarity stages.\n";
		print
"It is a NON-fatal issue.  Do you wish to use these results and utilise ESTScan? [Y/N]\n";
		my $answer = <STDIN>;
		if ( $answer =~ m/^\s*[nN]\s*$/ ) {
			return 0;
		}
		elsif ( $answer =~ m/^\s*[Yy]\s*$/ ) {
			$check = 1;
			next;
		}    #keep running process
		else { $check = 0; }    #go back an ask the question!
	}
	print "\n\n";
	return ( $matrix, $species );
}

sub count_codingNt {
	my $in = shift;
	open PRO, "<$in" or die "Cannot open: $in\n$!\n";
	my $file;
	while (<PRO>) {
		$file .= $_;
	}
	( my $num_coding = $file ) =~ s/.*masked\s*(\d+)\sof\s*(\d+)\s.*/$2-$1/se;
	unless ( $num_coding =~ m/^\d+$/ ) {    #so has to be a number
		print colored (
"### Error! ###\nThere appears to be a problem with the model building for ESTScan!\nPlease contact nematode.bioinf\@ed.ac.uk\n$2 - $1 = $num_coding\n",
			'red bold'
		);
		exit;
	}
	return $num_coding;
}

sub pseudo_collect {
	my $blast_ref = shift;
	my $species   = shift;
	my $embl_db   = shift;
	my $seeds     = shift;
	my $count     = 0;
	( my $name = $species ) =~ s/\s/_/g;
	while ( my ( $id, $hsps ) = each %$blast_ref ) {
		foreach my $hsp (@$hsps) {
			if ( $hsp->start('sbjct') == 1 ) {    #then we can use it
				my ( $q_start, $q_end ) = $hsp->range('query');
				my $seq = $seeds->{$id}->seq;
				write_cds( [ $q_start, $q_end, $seq ], $species, undef, $name );
				$count++;
			}
			else { next; }
		}
	}
	return ( $count, "$name.embl" );
}

sub write_cds {
	my ( $info, $name, $taxid, $out ) = @_;
	if ( ref $info ne 'ARRAY' ) {
		print "$_[0] needs to be an array_ref NOT ", ref $info, " !";
		exit;
	}

	my $io = Bio::SeqIO->new(
		-format => "embl",
		-file   => ">>$out.embl"
	);
	$taxid = 1 unless $taxid;
	my ( $start, $end, $seq ) = @$info;
	my $len = $end - $start;
	my @tags;
	push @tags, new Bio::SeqFeature::Generic(
		-start       => $start,
		-end         => $end,
		-primary_tag => 'source',
		-mol_type    => 'RNA'
	);
	push @tags, new Bio::SeqFeature::Generic(
		-start       => $start,
		-end         => $end,
		-primary_tag => 'CDS'
	);
	my $embl =
	  Bio::Seq::RichSeq->new( -display_id => "undef_cds", -seq => $seq );
	$embl->molecule('RNA');
	$embl->desc("standard; mRNA; UNK, $len BP");
	my $sp = Bio::Species->new( -ncbi_taxid => $taxid );
	my @class = split " ", $name;
	$sp->genus( shift @class );
	$sp->species( join " ", @class );
	$embl->species($sp);
	$embl->add_SeqFeature(@tags);
	$io->write_seq($embl);
}

sub build_matrix_with_pseudo {
	my $new_cds = shift;
	my $old_cds = shift;
	my $species = shift;
	my $id;
	my @counts;
	$counts[0] = '0';
	if ( -e $old_cds ) {    #so we have a file
		my $count = `grep -c "ID " $old_cds`;
		chomp($count);
		$counts[0] = $count;
		system
"rm -rf Evalute Isochores Matrices Report Shuffled Analysis mrna.seq Shuffled"
		  ;                 #need to delete old directories
		system "cat $new_cds $old_cds > cds.tmp";
		move( 'cds.tmp', $old_cds );
		($id) = $old_cds =~ m/([A-Za-z]+)/;
	}
	else {
		( $id = $species ) =~ s/^\s*(\w)\w+\s+(\w).*/$1$2/;
		$old_cds = $new_cds;
	}
	$counts[1] = `grep -c "ID " $old_cds`;
	print
"An additional $counts[0] sequences added to the training set to give $counts[1]\n";
	my $matrix;
	( $matrix, $species ) = &build_matrix( $species, $id, $old_cds );

	return $matrix;
}

sub run_estscan {

	my ( $estscan_exec, $infile, $matrix, $outNT, $outPEP, $len_lim ) = @_;

	unless ( -f $matrix ) {
		print
"Error cannot find matrix file for ESTScan!\nDoes not exist:$matrix\n";
		die;
	}
	unless ( -f $infile ) {
		print
		  "Error cannot find input file for ESTScan!\nDoes not exist:$infile\n";
		die;
	}
	system "$estscan_exec -M $matrix -l $len_lim -o $outNT -t $outPEP $infile";
}

sub parse_results {

	my ( $toESTScan, $ntfile, $pepfile, $seeds ) = @_;
	my $seqInfo;
	my %synom;
	@synom{ 'CT', 'GT', 'TC', 'CC', 'AC', 'GC', 'CG', 'GG' } =
	  ( 'L', 'V', 'S', 'P', 'T', 'A', 'R', 'G' );

	#these are the synomous codons
	#therefore ESTScan awards a amino acid.
	my %term;
	@term{ 'TAA', 'TAG', 'TGA' } = ( 1, 1, 1 );

	my $ntIO  = Bio::SeqIO->new( -file => $ntfile,  -format => 'Fasta' );
	my $pepIO = Bio::SeqIO->new( -file => $pepfile, -format => 'Fasta' );

	my %ntSeq;

	while ( my $seq = $ntIO->next_seq ) {
		my $id = $seq->display_id;
		$id =~ s/;//g;
		if ( exists $ntSeq{$id} ) {    #there are two possible translations
			    #one on each strand. Both have strong support
			$ntSeq{ $id . "_pos" } = $ntSeq{$id};
			$ntSeq{ $id . "_neg" } = $seq->seq;
			delete $ntSeq{$id};
		}
		else {
			$ntSeq{$id} = $seq->seq;
		}
	}

	while ( my $seqO = $pepIO->next_seq ) {
		my $id = $seqO->display_id;
		$id =~ s/;//g;
		my $cont   = 0;
		my $double = 0;
	  SWITCH: {
			if ( exists $ntSeq{$id} ) { $cont = 1; last SWITCH; }
			if ( exists $ntSeq{ $id . "_pos" } ) {
				$cont   = 1;
				$double = 1;
				last SWITCH;
			}
			if ( exists $ntSeq{ $id . "_neg" } ) {
				$cont   = 1;
				$double = 1;
				last SWITCH;
			}
		}
		next unless $cont == 1;
		my $desc = $seqO->desc;
		my $start;
		my $end;
		my $strand;

		if ( $desc =~ m/^-?\d+\s(\d+)\s(\d+)\s+(\w+)?/ ) {
			$start  = $1;
			$end    = $2;
			$strand = $3 ? '-1' : '1';
		}
		else { die; }

		my $id_tmp;
		if ($double) {
			$id_tmp = ( $strand == '-1' ) ? $id . '_neg' : $id . '_pos';
		}
		else {
			$id_tmp = $id;
		}
		my $pep   = $seqO->seq;
		my $len   = $seqO->length;
		my $ntide = $ntSeq{$id_tmp};
		$ntide =~ s/[acgtn]//g;
		my $len_nt  = length($ntide);
		my $seedlen = $seeds->{$id}->length;

		my $rem;
		my $x = () = ( $ntide =~ /X/g );
		my $y = () = ( $ntide =~ /[a-z]/g );

		while (1) {
			my $first_codon = substr $ntide, 0, 3;
			if ( $first_codon =~ m/(?:^[XN]\w\w)|(?:^\w[NX])/i ) {
				my $num = () = ( $first_codon =~ /[^XN]/g );
				$start += $num;
				$ntide =~ s/^$first_codon//g;
				redo;
			}
			elsif ( $first_codon =~ m/^(\w\w)[XN]/ ) {
				unless ( exists $synom{$1} ) {
					$start += length($1);
					$ntide =~ s/^$first_codon//g;
					redo;
				}
				last;
			}
			if ( exists $term{$first_codon} ) {
				$start += 3;
				$ntide =~ s/^$first_codon//g;    #first codon is a stop codon!
				redo;
			}
			last;
		}
		while (1) {
			my $last_codon = substr $ntide, -3;
			if ( $last_codon =~ m/[XN]{2}\w|[XN]\w\w|\w[XN]{2}/i ) {
				( my $nt = $& ) =~ s/X//g;
				$end -= length($nt);
				$ntide =~ s/$last_codon$//g;
				$pep   =~ s/X$//g;
				redo;
			}
			if ( $last_codon =~ m/^(\w\w)[XN]/ ) {
				unless ( exists $synom{$1} ) {
					$end -= length($1);
					$ntide =~ s/$last_codon$//g;
					$pep   =~ s/X$//g;
					redo;
				}
				last;
			}
			last;
		}
		$id =~ s/(_pos|_neg)//;
		push @{ $seqInfo->{$id} },
		  [
			[ $strand, $start ],
			[ $strand, $end ],
			undef, $pep, undef, undef, length($ntide), $ntide
		  ];   #record stats
		       #stored as an array of array as there may be more than one viable
		       #polypeptide produced by ESTScan (i.e. different strands.
		       #now check "quality" of sequence.
		my $max_10pc = sprintf( "%.0f", $len_nt * 0.1 );    #10% of $max_trans
		my $numX = () = ( $pep =~ m/X/g );

		if ( $len < 30 || $len < $max_10pc || ( ($numX) / $len > 0.05 ) ) {
			delete $seqInfo->{$id}[-1]
			  ;    # a bad 'un so deletes the last entry for that id.
		}
		if ( scalar @{ $seqInfo->{$id} } == 0 ) {    #do some tidying up...
			delete $seqInfo->{$id};
		}
	}
	return $seqInfo;
}

sub parse_results_mod {

	my ( $toESTScan, $ntfile, $pepfile, $seeds ) = @_;
	my $seqInfo;
	my %synom;
	@synom{ 'CT', 'GT', 'TC', 'CC', 'AC', 'GC', 'CG', 'GG' } =
	  ( 'L', 'V', 'S', 'P', 'T', 'A', 'R', 'G' );

	#these are the synomous codons
	#therefore ESTScan awards a amino acid.
	my %term;
	@term{ 'TAA', 'TAG', 'TGA' } = ( 1, 1, 1 );

	my $ntIO  = Bio::SeqIO->new( -file => $ntfile,  -format => 'Fasta' );
	my $pepIO = Bio::SeqIO->new( -file => $pepfile, -format => 'Fasta' );

	my %ntSeq;

	while ( my $seq = $ntIO->next_seq ) {
		my $id = $seq->display_id;

		# $id =~ s/;//g;

		# New estscan(3.0.3) seems to create headers differently
		# (or previous version actually did not work).
		# Example:

		#>cn35538; 49 1 164  LEN=164
		#XAAGCAGTGGTATCAACGCAGAGTGGCCATTACGGCCGGGGATCTCCCAACTTGTCCATG
		#GAGGTCTTCACTTCAGAAATCCAAGACTCATATTCATCCAGCTTGGTGTCAAGTGGGCTG
		#TTGCTGCCAGAATTATCTTGTGATTATTTGAGAGATGTATCAGTT
		#>cn35538;a 49 1 164  LEN=164; minus strand
		#XXAACTGATACATCTCTCAAATAATCACAAGATAATTCTGGCAGCAACAGCCCACTTGAC
		#ACCAAGCTGGATGAATATGAGTCTTGGATTTCTGAAGTGAAGACCTCCATGGACAAGTTG
		#GGAGATCCCCGGCCGTAATGGCCACTCTGCGTTGATACCACTGCTTXX
		# The ';a' needs to be removed not only the ';'
		$id =~ s/;.*//g;

		if ( exists $ntSeq{$id} ) {    #there are two possible translations
			    #one on each strand. Both have strong support
			$ntSeq{ $id . "_pos" } = $ntSeq{$id};
			$ntSeq{ $id . "_neg" } = $seq->seq;
			delete $ntSeq{$id};
		}
		else {
			$ntSeq{$id} = $seq->seq;
		}
	}

	while ( my $seqO = $pepIO->next_seq ) {
		my $id = $seqO->display_id;

		# $id =~ s/;//g;
		$id =~ s/;.*//g;    # hack, same reason as for nucleotide sequence above
		my $cont   = 0;
		my $double = 0;
	  SWITCH: {
			if ( exists $ntSeq{$id} ) { $cont = 1; last SWITCH; }
			if ( exists $ntSeq{ $id . "_pos" } ) {
				$cont   = 1;
				$double = 1;
				last SWITCH;
			}
			if ( exists $ntSeq{ $id . "_neg" } ) {
				$cont   = 1;
				$double = 1;
				last SWITCH;
			}
		}
		next unless $cont == 1;
		my $desc = $seqO->desc;
		my $start;
		my $end;
		my $strand;

		#		if ( $desc =~ m/^-?\d+\s(\d+)\s(\d+)\s+(\w+)?/ ) {
		#			$start  = $1;
		#			$end    = $2;
		#			$strand = $3 ? '-1' : '1';
		#		}
		#		else { die; }

		# The estscan (3.0.3) header seems to be different.
		# Method to look for strand direction changed.
		if ( $desc =~ m/^-?\d+\s(\d+)\s(\d+)\s+(.+)?/ ) {
			$start = $1;
			$end   = $2;
			my $other = $3;
			$strand = ( $other =~ m/minus/ ? -1 : 1 );
		}
		else {
			die;
		}

		my $id_tmp;
		if ($double) {
			$id_tmp = ( $strand == '-1' ) ? $id . '_neg' : $id . '_pos';
		}
		else {
			$id_tmp = $id;
		}
		my $pep   = $seqO->seq;
		my $len   = $seqO->length;
		my $ntide = $ntSeq{$id_tmp};
		$ntide =~ s/[acgtn]//g;
		my $len_nt = length($ntide);
		my $seedlen = $seeds->{$id}->length;

		my $rem;

		#		my $x=()=($ntide=~/X/g); # removed this
		#		my $y=()=($ntide=~/[a-z]/g);

		while (1) {
			my $first_codon = substr $ntide, 0, 3;
			if ( $first_codon =~ m/(?:^[XN]\w\w)|(?:^\w[NX])/i ) {
				my $num = () = ( $first_codon =~ /[^XN]/g );
				$start += $num;
				$ntide =~ s/^$first_codon//g;
				redo;
			}
			elsif ( $first_codon =~ m/^(\w\w)[XN]/ ) {
				unless ( exists $synom{$1} ) {
					$start += length($1);
					$ntide =~ s/^$first_codon//g;
					redo;
				}
				last;
			}
			if ( exists $term{$first_codon} ) {
				$start += 3;
				$ntide =~ s/^$first_codon//g;    #first codon is a stop codon!
				redo;
			}
			last;
		}
		while (1) {
			my $last_codon = substr $ntide, -3;
			if ( $last_codon =~ m/[XN]{2}\w|[XN]\w\w|\w[XN]{2}/i ) {
				( my $nt = $& ) =~ s/X//g;
				$end -= length($nt);
				$ntide =~ s/$last_codon$//g;
				$pep   =~ s/X$//g;
				redo;
			}
			if ( $last_codon =~ m/^(\w\w)[XN]/ ) {
				unless ( exists $synom{$1} ) {
					$end -= length($1);
					$ntide =~ s/$last_codon$//g;
					$pep   =~ s/X$//g;
					redo;
				}
				last;
			}
			last;
		}
		$id =~ s/(_pos|_neg)//;

#		push @{$seqInfo->{$id}},[[$strand,$start],[$strand,$end],undef,$pep,undef,undef,length($ntide),$ntide];	#record stats
#																#stored as an array of array as there may be more than one viable
#																#polypeptide produced by ESTScan (i.e. different strands.
# removed and added $ntSeq{$id_tmp} instead of $ntide
		push @{ $seqInfo->{$id} },
		  [
			[ $strand, $start ], [ $strand, $end ],
			undef,          $pep,
			undef,          undef,
			length($ntide), $ntSeq{$id_tmp}
		  ];

		#now check "quality" of sequence.
		my $max_10pc = sprintf( "%.0f", $len_nt * 0.1 );    #10% of $max_trans
		my $numX = () = ( $pep =~ m/X/g );
		if ( $len < 30 || $len < $max_10pc || ( ($numX) / $len > 0.05 ) ) {
			delete $seqInfo->{$id}[-1]
			  ;    # a bad 'un so deletes the last entry for that id.
		}
		if ( scalar @{ $seqInfo->{$id} } == 0 ) {    #do some tidying up...
			delete $seqInfo->{$id};
		}
	}
	return $seqInfo;
}

1;
