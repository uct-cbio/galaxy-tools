#!/usr/bin/perl

package p4e2::blaststuff;

$VERSION=1.0;
#10/03/04

use Term::ANSIColor;
use strict;

sub finddb	{	#much of this code has been modified from Al Anthony's
	my $dbtype=shift;
	my $ask=shift;
	
	#print "\$ENV{BLASTDB} = $ENV{BLASTDB}\n";
	
	unless (exists $ENV{BLASTDB})	{	#has this been set?
		print "Error: The BLASTDB environmental variable has not been set.\nPlease consult the blastall documentation\n";
		die;
	}	
	my @dbs;
	foreach my $path (split(":","$ENV{'BLASTDB'}"))	{
		if (-d $path && -r $path)	{
			if ($dbtype=~m/prot/)	{	#find protein databases
				push @dbs, glob("$path/*.pin");
			}
			elsif ($dbtype=~m/nuc/)	{	#find nucleotide databases
				push @dbs, glob("$path/*.nin");
			}
			elsif ($dbtype=~/both/)	{	#all of them!
				push @dbs, glob("$path/*.*in");
			}
			else {	#eh?
				print "Error: Unable to determine type of database selected\n";
				die;
			}
		}
	}
	if ($ask)	{ askwhich(\@dbs); }
	else { return \@dbs; }
}
###########
sub askwhich	{
	my $dbs=shift;
	
	my $db_idx;
	unless (ref $dbs eq 'ARRAY')	{
		print "Error: expecting an array refence\n";
		die;
	}
	
	my @dbs=@$dbs;	
	#so now offer the blastdb to the user
	
	if (scalar @dbs == '0')	{
		print "Error: No databases were found.\nPlease set the 'BLASTDB' environmental variable to the location of blast databases.\n";
		print "Consult the programs UserGuide for more information.\n";
		die;
	}
	else	{
		print colored ("The following databases have been found:\n",'bold');
		foreach my $n (0..$#dbs)	{
			(my $db = $dbs[$n]) =~ s/.*\/(.*)\.[pn]in/$1/;
			print $n+1,".  ",$db,"\n";
		}
	}
	my $flag;
	#print colored ("Please select one of the databases...\n",'bold');
	while ($flag==0)	{
		print colored ("Please select one of the databases...\n",'bold');
		chomp (my $choice=<STDIN>);
		unless (($choice =~m/\s*\d+\s*/) && ($choice <= ($#dbs+1)) && $choice > 0)	{	#make sure the option is valid
			print "INVALID CHOICE!\n";
			#print colored ("Please select one of the databases...\n",'bold');
			next;
		}
		else	{
			$flag=1;
			$db_idx=$choice-1;
		}
	}
	(my $database=$dbs[$db_idx])=~s/(.*)\.[pn]in/$1/;
	return $database;	
}

return 1;
