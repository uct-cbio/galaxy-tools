#!/usr/bin/perl

package p4e2::getDNA;

$VERSION = 1.1;

#22/11/04
#28/07/05 some changes to cope with ambiguous amino acids

use Term::ANSIColor;
use Bio::SeqIO;
use strict;

my %synom;
@synom{ 'CT', 'GT', 'TC', 'CC', 'AC', 'GC', 'CG', 'GG' } =
  ( 'L', 'V', 'S', 'P', 'T', 'A', 'R', 'G' );

#these are the synomous codons
#therefore ESTScan awards a amino acid.
my %term;
@term{ 'TAA', 'TAG', 'TGA' } = ( 1, 1, 1 );

sub fromBlast {

	#need to cycle through the HSPs collecting the nucleotides
	#from the sixframe. Then look at any extension regions and
	#fetch the DNA from these.

	my $seq_ref       = shift;
	my $nt_seq        = shift;
	my $code          = shift;
	my $gen_code_file = shift;
	my $failed_codon_comp = shift;

	for my $id ( keys %$seq_ref ) {
		my $seq   = $seq_ref->{$id};
		my $hsps  = $seq->[5];
		my $pep   = $seq->[3];             #this is the extended sequence
		my $tmp   = $seq->[3];
		my $frame = $seq->[0][0];
		my $ntide = $nt_seq->{$id}->seq;
		my @coding;
		foreach my $hsp (@$hsps) {
			my $start = $hsp->[0];
			my $end   = $hsp->[1];

			#print "\$start=$start\t$end\n$ntide\n";
			if ( $start > $end ) {
				print "Error: unable to return HSP $id => $start :: $end\n";
				die;
			}
			if ( $frame =~ /-/ ) {
				unshift @coding,
				  substr( $ntide, $start - 1, ( $end - $start + 1 ) );
			}
			else {
				push @coding,
				  substr( $ntide, $start - 1, ( $end - $start + 1 ) );
			}
		}
		if ( my $xtn_start = $seq->[0][2] ) {
			my $conf_start = $seq->[0][1];

			#print "1: $xtn_start\t$conf_start\n";
			if ( $frame =~ /-/ ) {
				push @coding,
				  substr(
					$ntide,
					( $xtn_start - 1 ),
					( $conf_start - $xtn_start )
				  );
			}
			else {
				unshift @coding,
				  substr(
					$ntide,
					( $xtn_start - 1 ),
					( $conf_start - $xtn_start )
				  );
			}
		}
		if ( my $xtn_end = $seq->[1][2] ) {
			my $conf_end = $seq->[1][1];

			#print "2: $conf_end\t$xtn_end\n";
			if ( $frame =~ /-/ ) {
				unshift @coding,
				  substr( $ntide, ($conf_end), ( $xtn_end - $conf_end ) );
			}
			else {
				push @coding,
				  substr( $ntide, ($conf_end), ( $xtn_end - $conf_end ) );
			}
		}
		my $coding_region;
		foreach my $seg (@coding) {
			if ( $frame =~ /-/ ) {
				$seg =~ tr/AGCTagct/TCGAtcga/;
				$seg = reverse($seg);
			}

			#			print "\$seg=$seg\n";
			my $seg_len = length($seg) / 3;
			my $pattern = "^\\w{$seg_len}(B+)";    #print "\$tmpBefore=$tmp\n";
			my $num_n   = 0;
			if ( $tmp =~ m/$pattern/ ) {

				#print "match=",$&,"\n";
				my $b = $1 || '';
				$num_n = 3 * length($b);
			}

			#print "correction factor is $seg_len+$num_n\n";
			$tmp = substr( $tmp, ( $seg_len + $num_n / 3 ) );

			#print "\$tmpAfter=$tmp\n";
			$coding_region .=
			  $seg . ( 'N' x $num_n );    #print "$num_n is a number\n\n\n";
		}
		my $gcode = store_gen_code($gen_code_file);
		$pep =~ s/B/X/g;
		&match_nt_aa( $coding_region, $pep, $gcode->[$code]->[1], $id, $failed_codon_comp );
		$seq->[7] = $coding_region;
	}
	return $seq_ref;
}

sub fromESTScan {

	my $seq_ref       = shift;
	my $estscan_file  = shift;
	my $ntide_seqs    = shift;
	my $code          = shift;
	my $gen_code_file = shift;
	my $failed_codon_comp = shift;

	my @ids = keys %$seq_ref;

	foreach my $id (@ids) {
		my $translns = $seq_ref->{$id};
		foreach my $transln (@$translns) {
			my $coding = $transln->[-1];
			my $pep    = $transln->[3];
			my $gcode  = store_gen_code($gen_code_file);
			&match_nt_aa( $coding, $pep, $gcode->[$code]->[1], $id, $failed_codon_comp );
		}
	}
	return $seq_ref;
}

sub fromDECODER {
	my $seq_ref       = shift;
	my $path2files    = shift;
	my $code          = shift;
	my $gen_code_file = shift;
	my $failed_codon_comp = shift;

	my @ids = keys %$seq_ref;
	foreach my $id (@ids) {
		( my $fileID = $id ) =~ s/\.Contig/_/;
		my $seqInfo = $seq_ref->{$id};
		my $frame   = $seqInfo->[0][0];
		my $dir     = ( $frame =~ m/-/ ) ? 3 : 5;
		my $ntIO    = Bio::SeqIO->new(
			-file   => "$path2files/Peptide.$dir/$fileID.tfa",
			-format => 'fasta'
		);
		my $pepIO = Bio::SeqIO->new(
			-file   => "$path2files/Peptide.$dir/$fileID.pep",
			-format => 'fasta'
		);
		my $start = $seqInfo->[0][1];
		my $end   = $seqInfo->[1][1];

		my $ntide;
		my $pep;
		while ( my $seqO = $ntIO->next_seq ) {
			$ntide = $seqO->seq;
		}
		while ( my $seqO = $pepIO->next_seq ) {
			$pep = $seqO->seq;
		}
		$ntide = substr( $ntide, $start - 1, ( $end - $start ) );
		my $gcode = store_gen_code($gen_code_file);
		&match_nt_aa( $ntide, $pep, $gcode->[$code]->[1], $id, $failed_codon_comp );
		$seqInfo->[7] = $ntide;
	}
	return $seq_ref;
}

sub fromORF {
	my $seq_ref       = shift;
	my $nt_seq        = shift;
	my $code          = shift;
	my $gen_code_file = shift;
	my $failed_codon_comp = shift;

	my @ids = keys %$seq_ref;
	foreach my $id (@ids) {
		my $seqInfo = $seq_ref->{$id};
		my $start   = $seqInfo->[0][1];
		my $end     = $seqInfo->[1][1];
		my $pep     = $seqInfo->[2];
		my $ntide   = $nt_seq->{$id}->seq;
		my $frame   = $seqInfo->[0][0];
		$ntide = substr $ntide, ( $start - 1 ), ( $end - $start + 1 );
		if ( $frame =~ /-/ ) {
			$ntide =~ tr/AGCTagct/TCGAtcga/;
			$ntide = reverse($ntide);
		}
		my $gcode = store_gen_code($gen_code_file);
		&match_nt_aa( $ntide, $pep, $gcode->[$code]->[1], $id, $failed_codon_comp );
		$seqInfo->[7] = $ntide;
	}
	return $seq_ref;
}

sub match_nt_aa {
	my ( $coding, $pep, $gcode, $id, $failed_codon_comp) = @_;
	my $tmp = $pep;
	my $x   = 0;
	while ( $coding =~ m/(\w{3})/g ) {
		my $codon = $1;
		my $pep = substr( $pep, 0, 1, "" );
		$x++;    #print "$pep($codon):";
		next if $codon =~ m/[^ACGT]/i;
		$codon = uc($codon);
		next if exists $term{$codon};

		next if $pep eq 'X';

		unless ( $gcode->{$codon} eq $pep ) {
			print "\n###error: Unable to match codon to amino acid residue\n",
			  "$id:\t$codon : ", $gcode->{$codon}, ":\t$pep\n";
			print "$x\n$tmp\n$coding\n";
			### added 2 lines ###
			$$failed_codon_comp = $$failed_codon_comp . $id . "\n";
			last;
			##################### 
#			die;
		}
	}
}

sub store_gen_code {
	my ( $file, $fetch ) = @_;
	local $/;
	open IN, "<$file" or die "Cannot open file: $file\n$!\n";
	my @gcodes;
	my $all = <IN>;
	while ( $all =~ m/\{\n(.*?)\}/sg ) {
		my $entry = $1;
		( my $name ) = $entry =~ m/name\s"(.*?)" ,/s;
		$name =~ s/\n//g;
		( my $id ) = $entry =~ m/id\s(\d+)\s,/;
		my ( @aa, @b1, @b2, @b3 );

		if ( $entry =~ m/\s+ncbieaa\s+"([\w\*]*)"/ ) {
			@aa = split //, $1;
		}
		if ( $entry =~ m/\s+--\sBase1\s+([AGCTU]*)/ ) {
			@b1 = split //, $1;
		}
		if ( $entry =~ m/\s+--\sBase2\s+([AGCTU]*)/ ) {
			@b2 = split //, $1;
		}
		if ( $entry =~ m/\s+--\sBase3\s+([AGCTU]*)/ ) {
			@b3 = split //, $1;
		}
		my %code;
		for ( my $i = 0 ; $i <= $#aa ; $i++ ) {
			my $triplet = join "", $b1[$i], $b2[$i], $b3[$i];
			$code{$triplet} = $aa[$i];
		}
		$gcodes[$id] = [ $name, \%code ];
	}
	close IN;
	return \@gcodes;
}
return 1;
