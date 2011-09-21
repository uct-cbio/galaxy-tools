#!/usr/bin/perl -w

package p4e2::p4e_checker;

$VERSION=1.4;
#changes: 
#10/03/03 - added cut_file subrountine.
#12/03/03 - allowed comments in cut file.
#10/08/03 - has sub to look for blast dbs of files.
#17/05/05 - regex bug fix

use strict;
use Term::ANSIColor;

##########
sub exe	{	#check that all the req/optional progs are in the user's path

	my @progs=@_;
	my @errors;
	my $found;
	
	unless (exists $ENV{ESTSCANDIR})	{
		print "\nTo ensure the smooth running of prot4EST you need to set the environmental variable ESTSCANDIR\n".
		      "This should be placed in the .cshrc / .bashrc file\n".
			  "For further queries please consult the user guide\n";
		exit;
	}
	else	{	  
		$ENV{PERL5LIB} .= ":".$ENV{ESTSCANDIR}."/MkTables";
		$ENV{PATH} .= ":".$ENV{ESTSCANDIR};
		$ENV{PATH} .= ":".$ENV{ESTSCANDIR}."/MkTables";
	}
	
	
	PROG: foreach my $prog (@progs)	{
		foreach my $path (split(":","$ENV{'PATH'}"))	{ 
			if (-x "$path/$prog")	{
				$ENV{PERL5LIB} .= ":$path ";
				$ENV{PATH} .= ":$path";
				next PROG;
			
			}
			elsif (-x "$path/MkTables/$prog")	{
				$ENV{PERL5LIB} .= ":$path/MkTables/";
				$ENV{PATH} .= ":$path/MkTables";
				next PROG;
			}
			else	{
				$found=0;
			}
		}
		push @errors, $prog if ($found==0);	
	} 
	return \@errors;
}

#########
sub sequence	{	#check that the sequences contain only legal characters

	my $file=shift;
	my $type=shift;
	
	my $input;
	my @errors;
	open IN, "<$file" or die;
	while (<IN>)	{
		$input .= $_;
	}
	close IN;
	
	my @seqs=split /^>/, $input;
	my $id = shift @seqs,"\n"; 						#removes the first and empty element from array
	
	foreach (@seqs)	{
		if ($type =~ m/IUPAC/i)	{ 		#IUPAC-approved single-letter base code
			unless (m/^(.*)?\n(?:[AGCTUNWYSDKHMRBV]+\n)+/si)	{
				push @errors, $1;
			}
		}
		elsif ($type =~ m/DNA/i)	{
			unless (m/^(.*)?\n(?:[AGCT]+\n)+/si)	{
				push @errors, $1;
			}		
		}
		elsif ($type =~ m/RNA/i)	{
			unless (m/^(.*)?\n(?:[AGCU]+\n)+/si)	{
				push @errors, $1;
			}		
		}
		elsif ($type =~ m/PROTEIN/i)	{
			if (!m/^(.*)?\n(?:[AC-IK-NP-TVWY]+\n)+/si)	{#;	
				push @errors, $1;			
			}
			
		}	
		else {print "Error: Unable to determine sequence type!\n"; die;}
	}
	return \@errors;
}
#########
sub bls_db	{	#look for blast dbs of files.
	
	my $file=shift;
	my $type=shift;
	my @errors;
	#my @nuc=qw/nhr nin nsq/;
	#my @pro=qw/phr pin psq/;
	if ($type eq 'nuc')	{
		foreach my $suf (qw/nhr nin nsq/)	{
			unless (-e $file.".$suf")	{
				push @errors, $file.".$suf";
			}
		}
	}
	elsif ($type eq 'pro')	{
		foreach my $suf (qw/phr pin psq/)	{
			unless (-e $file.".$suf")	{
				push @errors, $file.".$suf";
			}
		}
	}
	else 	{
		print "$type is not acceptable blast format!\n";
		die;
	}
	return \@errors;	
}
#########
sub species	{	#check that species name is in correct format
	
	my $species=shift;
	my @errors;
	unless ($species=~m/[A-Za-z]+\s+[A-Za-z]+(?:\s+(\w|.)+)?/)	{	#this should allow use of sub-species and varietas
		push @errors, $species;
	}
	return \@errors;
}
#########
sub blastfile	{	#is the file a valid blast report
	
	my $file=shift;
	open IN, "<$file" or die;
	my $count=0;
	my $okay=0;
	while ((my $line=<IN>) && $count++ < 5)	{	#only look at top 5 lines
		if ($line=~m/^(?:T|PSI)?BLAST[NXP]?/)	{
			$okay=1;
			return 'none';
		}
	}
	if (!$okay)	{
		return $file;
	}
}
#########
sub estscanmatrix	{
	
	local $/='FORMAT:';
	my @errors;
	my $file=shift;
	open IN, "<$file" or die;
	while (<IN>)	{
		next if (m/^FORMAT:/);
				
		unless (m/^[\/.:+\w\d\s]+?\n		#first line contains some info
		 (?:						#sets format for rest of string
			(?:-?\d{1,3}			#one or two digit number (may be negative)
				\s*){4}?\n)			#need four of them per line
		/sx)	
		{
			push @errors, $file;
		}
	}
	return \@errors;
}
#########
sub wheredecoder	{
	foreach my $path (split(":","$ENV{'PATH'}"))	{
		if (-x "$path/decoder")	{		#found it
			return "$path/decoder";
		}
		else	{next;}
	}
}
#########
sub decoder	{

	my ($sequaldir, $seqsuff, $qualsuff) = @_;
	
	my $error;
	my $qual_err;
	my %bothfiles;
	my @not4decoder;
	my @files=glob("$sequaldir/*");	#take all the files
	foreach my $file (@files)	{
		if ($file=~m/(.*)?\.$seqsuff/)	{
			$bothfiles{$1}->[0]=1;
			#should check that it is in fasta format
		}
		elsif ($file=~m/(.*)?\.$qualsuff/)	{
			$bothfiles{$1}->[1]=1;
			$qual_err = &check_phrap($file);
		}
		else { push @not4decoder, $file; $error=1;}
	}
#check that all entries have a sequence and quality file
	open DEC, ">decoder_errors.txt" or die "cannot open error report!\n$!\n";
	for my $key (keys %bothfiles)	{
		if ($bothfiles{$key}->[0] != '1')	{
			print DEC "No sequence file for $key\n";
			$error=1;
		}
		elsif ($bothfiles{$key}->[1] != '1')	{
			print DEC "No quality file for $key\n";
			$error=1;
		}
	}
	foreach (@not4decoder)	{
		print DEC "file: $_ invalid for DECODER program.  Please delete\n";
	}
	foreach (@$qual_err)	{
		print DEC "file: $_ is not a valid phrap quality file.  Please check\n";
		$error=1;
	}
	return $error;
}
#########
sub check_phrap	{	#looks at the quality file to make sure it is phrap format
	my $qualfile=shift;
	my @errors;

	open IN, "<$qualfile" or die "Cannot open quality file: $qualfile\n$!\n";
	LINE: while (<IN>)	{
		unless (m/^>/ || m/^\s?(?:\d+\s?)+\n/) 	{
			push @errors, $qualfile;
			last LINE;
		}
	}
	return \@errors;
}
#########
sub cut_file	{
	my $cut_in = shift;
	my $error;
	
	open CUT, "<$cut_in" or die;

	while ((!$error) && (my $line=<CUT>))	{
		if ($line=~m/^\s+$/)	{
			$error=0;
		}
		elsif ($line=~m/^\w{3}\s+\w{3}(?:\s+\d+\.\d\d){3}/)	{
			$error=0;
		}
		elsif ($line=~m/^#/)	{
			$error=0;
		}
		else	{
			$error=$line;
			last;
		}
	}
	return $error;
}
#########

return (1);
