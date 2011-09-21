#!/usr/bin/perl

$VERSION = 2.4;

#09/04/04 - extension process now deals with potential conflicts
#			caused through hits to minus strand.
#02/04/04 - extend now gives correct coords if minus strand.
#01/05/04 - replaces ALL X's from the blast search. New implementation.
#16/06/04 - bug fixed in replace X. No longer carries out replace if no X's exist
#05/10/04 - deals more stringently with repetitive sequence.
#05/10/04 - mysterious Js no longer appear in final sequence
#10/10/04 - bug fixed wrt location of polypeptide in EST
#24/11/04 - completely overhauled the tile_path method
#			tiling may give very slightly different results
#04/01/05 - a few minor bug fixes and streamlining of code.
#14/02/05 - tightened up extension process
#25/05/05 - does not try and include HSPs contained completely
#			within the tile_path - bug fix
#27/07/05 - bug fixed in building tile path routine '<' => '>'

package p4e2::tile_path;
use Term::ANSIColor;
use strict;

sub tile_path {
	my $hsp_ref      = shift;
	my $sixframe_ref = shift;
	my $gap_lim      = shift;
	my %sixframe_store;

	if ( ref $sixframe_ref eq 'ARRAY' ) {    #an aref
		foreach my $seqobj (@$sixframe_ref) {
			$sixframe_store{ $seqobj->display_id } = $seqobj;
		}
	}
	elsif ( ref $sixframe_ref eq 'HASH' ) {    #then a href
		%sixframe_store = %$sixframe_ref;
	}
	else {                                     #don't know!
		print colored( "Error: Unable to determine sixframe reference\n",
			'red bold' );
		die;
	}

	my %hsps = %$hsp_ref;                      #holds the hsp information

	my $seq_stats;
	my $total   = keys %hsps;
	my $current = 0;
	for my $id ( keys %hsps ) {                #key=id, query=set of hsps
		                                       #extract info from top hsp
		&progress( ++$current, $total );
		my $seq_len = pop @{ $hsps{$id} };
		unless ( $seq_len =~ m/^\d+$/ ) {
			print "\nERROR: $id - Unable to determine seq_len - $seq_len\n\n";
			die;
		}
		my $top_hsp = $hsps{$id}[0];
		unless ($top_hsp) {
			print "\nERROR: $id - Unable to tile fragment\n\n";
			die;
		}

		my $t_strand    = $top_hsp->query->strand;
		my $t_frame     = ( $top_hsp->query->frame + 1 ) * $t_strand;
		my $t_start_hit = $top_hsp->hit->start;
		my $t_end_hit   = $top_hsp->hit->end;
		my $t_start_sbj = $top_hsp->query->start;
		my $t_end_sbj   = $top_hsp->query->end;
		my $t_e         = $top_hsp->evalue;
		my $t_score     = $top_hsp->score;
		my $t_seq       = $top_hsp->query_string;
		$t_seq =~ s/-//g;    #remove gaps

   #the bioperl assignment of hsp start/end are w.r.t to the hsp.
   #with regard to the translation the hsps must be considered to be in the same
   #orientation.
   #so:  		start---------------------->end
   #					end<--------------------start
   #
   #becomes:	start---------------------->end
   #					start<------------------end
   #
   #the ->query_string is already corrected.

#with filtering turned on it is likely that there are X's in the sequence.  These need to be replaced with true sequence.
		my $tile_path;
		if ( $t_seq =~ m/X/ ) {
			my $locator
			  ;    #this will be the end or the start depending on the strand.
			if ( $t_strand == '-1' ) {
				$locator = $seq_len - $t_end_sbj;
			}
			elsif ( $t_strand == '1' ) {
				$locator = $t_start_sbj;
			}
			else {
				print "Cannot determine the strand!\n";
				die;
			}
			$tile_path = replace_x( $top_hsp->query_string, $id, $t_frame,
				\%sixframe_store, $locator );
		}
		else {
			$tile_path = $t_seq;
		}

		my $complete = 0;

		#hmmmmm can this be made to look tidier?
		$seq_stats->{$id} = [
			[ $t_frame, $t_start_sbj ],
			[ $t_frame, $t_end_sbj ],
			[],
			[],
			[],
			[
				[
					$t_start_sbj, $t_end_sbj,   $t_frame, $t_e,
					$t_score,     $t_start_hit, $t_end_hit
				]
			],
			$seq_len
		];    #create an anonymous hash of a multi-D array
		 #print "$t_start_sbj, $t_end_sbj, $t_frame, $t_e, $t_score,$t_start_hit,$t_end_hit\n";

		#[0] - the start position of the tile path with regard to the EST
		#[1] - the end position of the tile path with regard to the EST
		#[2] - the amino acid sequence from the tile path process
		#[3] - the full length sequence after tile path and extension
		#[4] - notes/comments
		#[5] - anonymous array which stores information relating to each hsp
		#[6] - length of nucleotide sequence
		#[7] - the ntide coding region that is added later.
		my %seen;
	  HSPS: while ( !$complete ) {    #is the whole process complete 1-no, 0-yes
			for ( my $i = 1 ; $i <= $#{ $hsps{$id} } ; $i++ )
			{                         #go thru each hsp for this query (EST)
				my $new_hsp = $hsps{$id}->[$i];
				next if exists $seen{$i};    #seen and used it before
				my $n_strand    = $new_hsp->query->strand;
				my $n_frame     = ( $new_hsp->query->frame + 1 ) * $n_strand;
				my $n_start_hit = $new_hsp->hit->start;
				my $n_end_hit   = $new_hsp->hit->end;
				my $n_start_sbj = $new_hsp->query->start;
				my $n_end_sbj   = $new_hsp->query->end;
				my $n_e         = $new_hsp->evalue;
				my $n_score     = $new_hsp->score;
				my $n_seq;

				if ( $new_hsp->query_string =~ m/X/ ) {
					my $locator
					  ; #this will be the end or the start depending on the strand.
					if ( $n_strand == '-1' ) {
						$locator = $seq_len - $n_end_sbj;
					}
					elsif ( $n_strand == '1' ) {
						$locator = $n_start_sbj;
					}
					else {
						print "Cannot determine the strand!\n";
						die;
					}
					$n_seq = replace_x( $new_hsp->query_string, $id, $n_frame,
						\%sixframe_store, $locator );
				}
				else {
					$n_seq = $new_hsp->query_string;
				}
				$n_seq =~ s/-//g;
				my $n_len_aa = length($n_seq);

				my @t_hsp;    #need to find the 5' most hsp for comparison.
				if ( $t_frame =~ /-/ ) {
					@t_hsp = @{ $seq_stats->{$id}->[5]->[-1] };
				}
				else {
					@t_hsp = @{ $seq_stats->{$id}->[5]->[0] };
				}
				my ( $t_start_sbj, $t_end_sbj, $t_frame, $t_e, $t_score,
					$t_start_hit, $t_end_hit )
				  = @t_hsp;

				if ( $t_strand eq $n_strand ) {    #possible frameshifts
					    #can a nearby hsp extend the profile upstream?
					my $diff5 =
					  $t_start_hit -
					  $n_end_hit
					  ;    #calculate the difference between start & stop points
					my ( $dir5, $mod5 ) = ( $diff5 =~ m/(-?)(\d+)/ )
					  ;    #need the size of difference not direction.
					if ( $mod5 <= $gap_lim )
					{ #so need to compare the new hsp with 5' hsp currently in the tile.
						next if ( $t_frame eq $n_frame );
						next
						  if ( $t_start_hit - $n_start_hit <= 1 )
						  ; #if an hsp is short it may be completely within the tile_path
						my $extension = '';
						my $filler    = '';
						if ( $dir5 eq '-' || $mod5 == 0 )
						{    #overlap so need to crop from one sequence
							my $crop_aa = $mod5 + 1;
							my $crop_nt = $crop_aa * 3;
							if ( $t_e <= $n_e ) {    #crop from incoming hsp
								if ( $n_strand eq '-1' ) {
									$seq_stats = &add_hsp(
										$t_strand,
										[
											$n_start_sbj + $crop_nt,
											$n_end_sbj,
											$n_frame,
											$n_e,
											$n_score,
											$n_start_hit + $crop_aa,
											$n_end_hit
										],
										$seq_stats,
										$id
									);
								}
								else {
									$seq_stats = &add_hsp(
										$t_strand,
										[
											$n_start_sbj,
											$n_end_sbj - $crop_nt,
											$n_frame,
											$n_e,
											$n_score,
											$n_start_hit,
											$n_end_hit - $crop_aa
										],
										$seq_stats,
										$id
									);
								}
								$extension = substr $n_seq, 0,
								  ( $n_len_aa - $crop_aa );
							}
							else {    #remove from previously 5' most hsp

							$tile_path=substr $tile_path, ($crop_aa);
#							$tile_path=substr $tile_path, ($crop_aa); # REMOVED THIS LINE HERE FOR NOW
								$seq_stats->{$id}->[5]->[0]->[0] += $crop_nt
								  ; #make start coord correction in current 5' hsp
								$seq_stats->{$id}->[5]->[0]->[5] += $crop_aa;
								$seq_stats = &add_hsp(
									$t_strand,
									[
										$n_start_sbj, $n_end_sbj,
										$n_frame,     $n_e,
										$n_score,     $n_start_hit,
										$n_end_hit
									],
									$seq_stats,
									$id
								);
								$extension = $n_seq;
							}
						}
						else
						{ #the new hsp and current tile_path are very close but do not overlap
							    #this means filling in is required.
							my $add_aa = $mod5 - 1;
							$filler    = 'B' x $add_aa;
							$seq_stats = &add_hsp(
								$t_strand,
								[
									$n_start_sbj, $n_end_sbj, $n_frame,
									$n_e,         $n_score,   $n_start_hit,
									$n_start_hit, $n_end_hit
								],
								$seq_stats,
								$id
							);
							$extension = $n_seq;
						}

						$tile_path = $extension . $filler . $tile_path;

						if ( $t_strand =~ m/-/ ) {
							$seq_stats->{$id}->[1] = [ $n_frame, $n_end_sbj ];
						}
						else {
							$seq_stats->{$id}->[0] = [ $n_frame, $n_start_sbj ];
						}
						$seen{$i} = '';
						redo HSPS
						  ; #go back and see if a previously considered (yet ignored)
						    #HSP can move the tile_path more 5'wards.
					}
				}

				#now lets look downstream
				#need to find the 3' most hsp for comparison.
				if ( $t_frame =~ /-/ ) {
					@t_hsp = @{ $seq_stats->{$id}->[5]->[0] };
				}
				else {
					@t_hsp = @{ $seq_stats->{$id}->[5]->[-1] };
				}
				(
					$t_start_sbj, $t_end_sbj, $t_frame, $t_e, $t_score,
					$t_start_hit, $t_end_hit
				) = @t_hsp;
				if ( $t_strand eq $n_strand ) {    #possible frameshift;
					next if ( $t_frame eq $n_frame );
					my $diff3 =
					  $n_start_hit -
					  $t_end_hit
					  ;    #calculate the difference between start & stop points
					my ( $dir3, $mod3 ) = ( $diff3 =~ m/(-?)(\d+)/ );
					if ( $mod3 <= $gap_lim ) {    #looks like a frameshift event
						 #need to make sure that the "new" hsp is not contained completely within the tile path.
						 #Can happen for short matches.
						next if ( $t_end_hit - $n_end_hit >= 1 );
						my $extension = '';
						my $filler    = '';

						if ( $dir3 eq '-' || $mod3 == 0 )
						{    #overlap so remove from one of the hsps
							my $crop_aa = $mod3 + 1;
							my $crop_nt = $crop_aa * 3;
							if ( $t_e <= $n_e ) {    #crop of new sequence
								if ( $n_strand =~ /-/ ) {
									$n_end_sbj -= $crop_nt;
								}
								else {
									$n_start_sbj += $crop_nt;
								}
								$seq_stats = &add_hsp(
									( $t_strand * -1 ),
									[
										$n_start_sbj, $n_end_sbj,
										$n_frame,     $n_e,
										$n_score,     $n_start_hit + $crop_aa,
										$n_end_hit
									],
									$seq_stats,
									$id
								);
								$extension = substr $n_seq, $crop_aa;
							}
							else {    #remove overhand from current 3' hsp
								$seq_stats->{$id}->[5]->[-1]->[1] -= $crop_nt;
								$seq_stats->{$id}->[5]->[-1]->[6] -= $crop_aa;
								$seq_stats = &add_hsp(
									( $t_strand * -1 ),
									[
										$n_start_sbj, $n_end_sbj,
										$n_frame,     $n_e,
										$n_score,     $n_start_hit,
										$n_end_hit
									],
									$seq_stats,
									$id
								);
								$tile_path = substr $tile_path, 0,
								  ( length($tile_path) - $crop_aa );
								$extension = $n_seq;
							}

						}
						else
						{ #hsps are close but do not overlap. Gaps needs filling
							my $add_aa = $mod3 - 1;
							$filler    = 'B' x $add_aa;
							$seq_stats = &add_hsp(
								( $t_strand * -1 ),
								[
									$n_start_sbj, $n_end_sbj, $n_frame,
									$n_e,         $n_score,   $n_start_hit,
									$n_end_hit
								],
								$seq_stats,
								$id
							);
							$extension = $n_seq;
						}
						$tile_path .= $filler . $extension;
						if ( $t_strand =~ m/-/ ) {
							$seq_stats->{$id}->[0] = [ $n_frame, $n_start_sbj ];
						}
						else {
							$seq_stats->{$id}->[1] = [ $n_frame, $n_end_sbj ];
						}
						$seen{$i} = '';
						redo HSPS;
					}
				}
			}

			$complete = 1;    #process is complete
		}
		$seq_stats->{$id}->[2] = $tile_path;
	}

	return ($seq_stats);
}
#############
sub add_hsp {
	my $strand    = shift;
	my $data_ref  = shift;
	my $seq_stats = shift;
	my $id        = shift;
	if ( $strand =~ m/-/ ) {
		push @{ $seq_stats->{$id}->[5] }, $data_ref;
	}
	else {
		unshift @{ $seq_stats->{$id}->[5] }, $data_ref;
	}
	my $s = scalar @{ $seq_stats->{$id}->[5] };
	return $seq_stats;
}

#############
sub replace_x {
	my $seq       = shift;
	my $id        = shift;
	my $frame_seq = shift;
	my $sixframe  = shift;
	my $locator   = shift;

	my %seqs_6;
	if ( ref $sixframe eq 'HASH' ) {    #is an aref
		%seqs_6 = %$sixframe;
	}
	elsif ( -f $sixframe ) {            #is a file
	}
	else {                              #don't know!
		print colored (
"\nError: Unable to determine sixframe reference\nExpecting a hash reference or valid file name\n\n",
			'red bold'
		);
		die;
	}

	$seq =~ s/-//g;                     #remove gaps first
	my ( $termN, undef, $termC ) = $seq =~ /^([\w\*]{3})(.*)([\w\*]{3})$/;
	my $len = $+[2] - $-[2];

	if ( !$termN || !$termC ) {
		warn
		  "Error in replacing filtered amino acids: $id-  :$termN:  :$termC:\n";
		exit;
	}
	my ( $strand, $frame ) = ( $frame_seq =~ m/(-?)(\d)/ );
	if ( $strand =~ m/-/ ) {
		$frame = ($frame) + 3;
	}
	my $six_seq = $sixframe->{ "$id" . "_$frame" }->seq;
	unless ($six_seq) {
		warn "Error in replacing filtered amino acids: $id";
		die;
	}
	foreach ( $termN, $termC ) {
		$_ =~ s/X/\./ig;
		$_ =~ s/\*/\\*/g;
	}

	my $i      = 0;
	my $bridge = "." x $len;
	while ( $six_seq =~ m/($termN$bridge$termC)/g ) {
		my $match = $1;
		my $location = $-[0] - ( $locator / 3 );
		$location =~ s/-//;
		next unless ( $location < 4 );
		$i++;
		if ( $i != 1 ) {
			warn
"Error in replacing filtered amino acids: $id $i\nsearch pattern not unique!\n";
			exit;
		}
		unless ( length($seq) == length($match) ) {
			print "\nProblem replacing filtered sequence\n$seq\n\n$match\n";
			exit;
		}
		while ( $seq =~ m/X+/g ) {
			my $l_old = length($seq);
			my $line  = $&;
			my $l     = length($line);
			my $s     = $-[0];
			my $new   = substr( $match, $s, $l );
			$new =~ s/X$/J/;         #terminal amino acid
			$seq =~ s/$line/$new/;
			unless ( length($seq) == $l_old ) {
				warn
"Error replacing filtered amino acids: $id\t new_string and old_strings are different lengths\n",
				  length($seq), "\t$l\n";
				exit;
			}
		}
	}
	$seq =~ s/J/X/g
	  ; #some of the sequence may have to be an X because a nucleotide in the codon may be X.
	return $seq;
}
########
sub extend {
	my $seq_ref      = shift;
	my $sixframe_ref = shift;
	my $size         = keys %$seq_ref;
	my %sixframe_store;
	my $method = shift;

	#the method is necessary.
	#if blast then there's the possiblity that different HSPs
	#in the tiling path have opposing frameshifts.
	#this is fine but we don't want to extend these,
	#as they'll extend into one another.

	if ( ref $sixframe_ref eq 'ARRAY' ) {    #an aref
		foreach my $seqobj (@$sixframe_ref) {
			$sixframe_store{ $seqobj->display_id } = $seqobj;
		}
	}
	elsif ( ref $sixframe_ref eq 'HASH' ) {    #then a href
		%sixframe_store = %$sixframe_ref;
	}
	else {                                     #don't know!
		print colored( "Error: Unable to determine sixframe reference\n",
			'red bold' );
		die;
	}

	my $count = 1;

	#extract information from $seq_ref
	for my $query_id ( keys %$seq_ref ) {
		my $frame5        = $seq_ref->{$query_id}[0][0];
		my $term5_nuc     = $seq_ref->{$query_id}[0][1];
		my $nt_len        = $seq_ref->{$query_id}[6];
		my $frame3        = $seq_ref->{$query_id}[1][0];
		my $term3_nuc     = $seq_ref->{$query_id}[1][1];
		my $seq_sofar     = $seq_ref->{$query_id}[2];
		my $hsp_array_ref = $seq_ref->{$query_id}[5]
		  if $seq_ref->{$query_id}[5];
		my $term5_pro = $seq_ref->{$query_id}[5][0][0];
		my $term3_pro = $seq_ref->{$query_id}[5][-1][1];

		if ( $frame5 =~ m/-/ ) {    #then minus frame so need to mess about
			my $old5 = $term5_pro;
			$term5_pro = $nt_len - $term3_pro;
			$term3_pro += ( $term5_pro - $old5 );
		}

		#lets get some info relating to the 5' end of the sequence
		#remember that this is w.r.t the peptide but the coords are w.r.t.
		#the EST.  These may need to be altered.

		my ( $rel_start, $rel_end );
		my ( $runin5,    $runin3 );
		my ( $ext_5,     $ext_3 );     #these will hold any additional sequence
		$ext_5 = '';
		$ext_3 = '';
		my ($ext_seq);

		my $minus5;
		if ( $frame5 =~ m/-/ ) {       #so this is the reverse strand

#need to check whether extending in the either end direction is likely to write into
#another part of the sequence.
#this is the case for sequences built from more than one HSP
#
#if the start/end of the sequence is from the minus strand and there is a
#hsp from the positive strand anywhere in the prediction there is a
#risk that these will collide
#for now if this is the case then there will be no extension.
#later version will give a more precise solution.

			$runin5 = &check_runins( $hsp_array_ref, 'front' );

			$frame5    = ( $frame5 * -1 ) + 3;
			$minus5    = 1;
			$rel_start = $nt_len - $term5_nuc + 1;
		}
		else { $frame5 =~ s!\+!!; $minus5 = 0; $rel_start = $term5_nuc; }

		#now the 3' end w.r.t to the peptide

		my $minus3;
		if ( $frame3 =~ m/-/ ) {    #so this is the reverse strand

			#look at long comment above.
			$runin3  = &check_runins( $hsp_array_ref, 'front' );
			$frame3  = ( $frame3 * -1 ) + 3;
			$minus3  = 1;
			$rel_end = $nt_len - $term3_nuc + 1;
		}
		else { $frame3 =~ s!\+!!; $minus3 = 0; $rel_end = $term3_nuc; }

		my $tile_len = length($seq_sofar);
		my ( $first_aa, $last_aa ) = ( $seq_sofar =~ m/^([\w*]).*([\w*])$/ )
		  ;    #find the aa residues at both end of the sequence

		#first consider 5' direction...
		unless ($runin5) {    #0 for no conflict
			my $seq5 = $sixframe_store{ $query_id . "_" . $frame5 }->seq;
			my $region5 = substr( $seq_sofar, 0, 30 );
			$region5 = quotemeta($region5);
			my $upstream_seq;
			my $count_x;
			while ( $seq5 =~ m/$region5/g ) {
				my $aa_start = $-[0];
				$term5_pro /= 3;
				$term5_pro -= $aa_start;
				$term5_pro =~ s/-//g;
				next unless $term5_pro < 4;
				if ( $seq5 =~ m/(.*?)$region5.*/ )
				{    # 2(30^20) chance of this sequence appearing twice
					$upstream_seq = $1
					  ; #use this if loop.  If it doesn't match then the 5' end is not a perfect match to
					 #the original EST, so we are unsure as to the actual start.  This is used in 3' later.
				}
				else { die; }
				if ( length($upstream_seq) > 0 )
				{    #so exists and has some length
					    #identify the possible start of the sequence
					my $whereMet = -1;
					my @mets;
					my $new_met;
					while ()
					{    #find location of all METs in the upstream section
						$whereMet = index( $upstream_seq, "M", $whereMet + 1 );
						if ( $whereMet == -1 ) {    #no more METs
							last;
						}
						else { push @mets, $whereMet; }
					}
					my $stop = rindex( $upstream_seq, "*" )
					  ;    #find location of last stop codon

					if ( !@mets ) {    #die("here1");
						if ( ( $term5_nuc <= 90 ) && ( $first_aa =~ m/M/i ) )
						{    #the hsp is near the start of the subject and
							 #the first amino acid in the hsp is a Met.     |---MALTWHISKEY
							$ext_5 = ''
							  ; #so no additional sequence to append to the 5' end.
						}
						else
						{ #even we do have a Met as first_aa it is unlikely to be the true one,
							 #so use all the upstream sequence in this frame.   /-------------MALTWHISKEY
							$ext_5 = $upstream_seq
							  ;    #append remainder of upstream sequence.

							if ( $first_aa eq 'M' ) {
								push @{ $seq_ref->{$query_id}[4] },
								  (
"Possible leading Met found at peptide position 1"
								  );
							}
							else {    #no leading Met found in sequence  :o(
								push @{ $seq_ref->{$query_id}[4] },
								  'No leading Met found'
								  ;    #add this information to seq_ref
							}

						 #remove any sequence upstream of the 3' most stop codon
							$ext_5 =~ s/(?:.*\*)?(.*)$/$1/;
						}
					}

					elsif ( @mets && $stop )
					{ #there is both a Met and stop codon in the upstream sequence.
						 #need to find if there is a Met downstream of the stop codon.
						 #die("here2");
						foreach my $met (@mets) {
							if ( $met < $stop )
							{ #there is another putative Met upstream in this frame but with an interceeding stop codon
								$new_met = 'n'
								  ; #this also covers no Met.		/-----M---*---MALTWHISKEY
								next;
							}
							else
							{ #there is a Met which maps out sequence not interupted by a stop codon
								$new_met = $met;
								last;    #only record most 5 prime one.
							}
						}
						if ( $new_met eq 'n' ) {
							$ext_5 = ''
							  ; #so no additional sequence to append to the 5' end.

							if ( $first_aa ne 'M' ) {
								push @{ $seq_ref->{$query_id}[4] },
								  'No leading Met found'
								  ;    #add this information to seq_ref
							}
						}
						else {         #	|----*--M----MALTWHISKEY
							$ext_5 = substr( $upstream_seq, $new_met )
							  ;        #store from the new Met downwards

							if ( $first_aa eq 'M' ) {
								push @{ $seq_ref->{$query_id}[4] },
"Possible leading Met at peptide position 1 (extension process), peptide position "
								  . ( 1 + length($ext_5) )
								  . "(HSP tile path)";
							}
							else {
								push @{ $seq_ref->{$query_id}[4] },
'Possible leading Met at position 1(extension process)';
							}
						}
					}
					elsif ( @mets && $stop == '-1' )
					{    #upstream Met and no stop codon.
						    #die("here3");
						if ( $mets[0] > 5 ) {    #	/-------M------(M)ALTWHISKEY
							$ext_5 = $upstream_seq
							  ; #add all upstream but make note of this putative Met
							    #no need to alter {$query}[0][1]
							if ( $first_aa eq 'M' ) {
								push @{ $seq_ref->{$query_id}[4] },
								  (
"Probable leading Met found at peptide position "
									  . ( 1 + length($ext_5) )
									  . " Possible leading Met found at ",
									( $mets[0] + 1 )
								  );
							}
							else {
								push @{ $seq_ref->{$query_id}[4] },
								  (
									"Possible leading Met found at ",
									$mets[0] + 1
								  );
							}

						}
						else {    #|--M------WHISKEY
							$ext_5 = substr( $upstream_seq, $mets[0] );
							if ( $first_aa eq 'M' ) {
								push @{ $seq_ref->{$query_id}[4] },
								  ( "Likely Met found at peptide position "
									  . 1 + length($ext_5) );
							}
						}
					}
					$ext_5 =~ s/^X//g;
					my $len5 = length($ext_5);
					if ( $len5 == '0' ) {    #die("here4");
						$seq_ref->{$query_id}[0][2] = undef;
						if ( $seq_sofar =~ m/^M/ ) {    #starting Met
							push @{ $seq_ref->{$query_id}[4] },
							  (
"Probable leading Met found at peptide position 1"
							  );
						}
					}
					else {
						if ($minus5) {                  #die("here5");
							my $coord = ( $len5 * 3 + $term3_nuc );
							$coord = 1
							  unless ( $coord > 0 )
							  ;    #compensates for bug in ESTScan
							$seq_ref->{$query_id}[1][2] = $coord;
						}
						else {
							my $coord = ( $term5_nuc - $len5 * 3 );
							$coord = 1
							  unless ( $coord > 0 )
							  ;    #compensates for bug in ESTScan
							$seq_ref->{$query_id}[0][2] = $coord;
						}
					}
				}
				else {             #no extension in 5' direction
					$seq_ref->{$query_id}[0][2] = '';
					if ( $seq_sofar =~ m/^M/ ) {    #starting Met
						push @{ $seq_ref->{$query_id}[4] },
						  ("Probable leading Met found at peptide position 1");
					}
				}
			}
		}

		unless ($runin3) {

			#Now consider 3' direction

			my $seq3 = $sixframe_store{ $query_id . "_" . $frame3 }->seq;
			my $region3 = substr( $seq_sofar, -30 );
			$region3 = quotemeta($region3);
			my $downstream_seq;
			my $count_x = 0;

			#++$count_x while $seq3=~m/$region3/g;
			while ( $seq3 =~ m/$region3/g ) {
				my $aa_end = $+[0];
				$term3_pro /= 3;
				$term3_pro -= $aa_end;
				$term3_pro =~ s/-//g;
				next unless $term3_pro < 4;
				if ( $seq3 =~ m/.*?$region3(.*)/ ) {
					$downstream_seq = $1;
				}
				if ( length($downstream_seq) > 0 ) {    #so exists and has size
					my $stop_pos = index( $downstream_seq, "\*" )
					  ;    #find the putatuve stop codon for the translation

					if (   ( length($downstream_seq) < 30 )
						&& ( $stop_pos == 0 ) )
					{

#so the blastx search has predicted a stop codon at the end and there are less than 30 amino acids (90 nucleotides)
#downstream.  This is likely to changed in version 2.0
						$ext_3 = '';    #nothing to extend
					}
					elsif ( $stop_pos > 0 )
					{                   #a stop codon downstream  SKEY-----*---|
						$ext_3 = substr( $downstream_seq, 0, $stop_pos )
						  ;    #will add this extra seqeunce to prediction.
					}
					elsif ( $stop_pos == -1 )
					{          #no stop codon anywhere		SKEY---------/
						$ext_3 = $downstream_seq;    #add the whole 3' section.
					}
				}
				$ext_3 =~ s/X$//g;
				my $len3 = length($ext_3);
				if ( $len3 == '0' ) { }
				elsif ($minus3) {
					$seq_ref->{$query_id}[0][2] =
					  ( $term5_nuc - length($ext_3) * 3 + 1 );
				}
				else { $seq_ref->{$query_id}[1][2] = $len3 * 3 + $term3_nuc; }
			}
		}
		$seq_sofar =~ s/\*$/J/;    #temporary change of any final stop codon
		$seq_sofar =~
		  s/\*/X/g;    #replace any stop codons from the original HSP with X's
		$seq_sofar =~ s/J$/\*/;    #now put it back
		                           #so now set the new putative translation
		$ext_seq = $ext_5 . $seq_sofar . $ext_3;
		$seq_ref->{$query_id}->[3] = $ext_seq;
	}
	return $seq_ref;
}

sub check_runins {
	my $hsps_ref = shift;
	my $where    = shift;
	my ( $termHsp, $confl );

	if ( $where eq 'front' ) {     #need the fist hsp
		$termHsp = $hsps_ref->[0];
	}
	elsif ( $where eq 'end' ) {    #need last hsp
		$termHsp = $hsps_ref->[-1];
	}
	else {                         #??
		print "\nERROR!:\nCannot determine nature of HSP order!\n\n";
		exit;
	}

	#what strand is the terminal hsp on?
	if ( $termHsp->[2] =~ m/-/ ) {    #need to check possible conflicts
		foreach my $hsp (@$hsps_ref) {
			if ( $hsp->[2] !~ m/-/ ) {
				$confl = 1;
				last;                 #only need to happen once
			}
		}
	}

	#else its fine to extend.
	return $confl;
}

sub progress {
	my $current    = shift;
	my $total      = shift;
	my $percentage = ( $current / $total ) * 100;
	my $progress   = sprintf "%5.2f", $percentage;
	print "\r$progress%";
}

return (1);
