#! /usr/bin/perl -w
package Utils::InterProScanUtils;
use strict;
use warnings;

########### Split iprscan xml file into seperate interproscan xml files #############
################### e.g. split_interproscan( "input.xml", "output" ) ################
sub split_interproscan() {
    my $file   = shift(@_);
    my $dir    = shift(@_);
    my $flag   = 0;
    my $in_seq = "";
    my @heading;
    my $header;
    my $filename;
    open( FILE, "$file" ) || die "Can't open $file\n";
    my $content = "";

    while ( my $line = <FILE> ) {
        
        $content .= $line;

        if ( $line =~ /<protein id="(.*?)"/ ) {
            $filename = $1;
            $content = $line;            
        }

        if (( $line =~ /<\/protein>/ ) && $filename) {
            open( OUTFILE, ">$dir/$filename.xml" )
              || die "Can't open $dir/$filename.xml\n";
            print OUTFILE "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>
<EBIInterProScanResults   >
        <Header>
                <program name=\"InterProScan\" version=\"4.0\" citation=\"PMID:11590104\" />
                <parameters>
                        <sequences total=\"1\" />
                        <databases total=\"13\">
                                <database number=\"1\" name=\"PRODOM\" type=\"sequences\" />
                                <database number=\"2\" name=\"PRINTS\" type=\"matrix\" />
                                <database number=\"3\" name=\"PIR\" type=\"model\" />
                                <database number=\"4\" name=\"PFAM\" type=\"model\" />
                                <database number=\"5\" name=\"SMART\" type=\"model\" />
                                <database number=\"6\" name=\"TIGRFAMs\" type=\"model\" />
                                <database number=\"7\" name=\"PROFILE\" type=\"strings\" />
                                <database number=\"8\" name=\"PROSITE\" type=\"strings\" />
                                <database number=\"9\" name=\"SUPERFAMILY\" type=\"model\" />
                                <database number=\"10\" name=\"GENE3D\" type=\"model\" />
                                <database number=\"11\" name=\"PANTHER\" type=\"model\" />
                                <database number=\"12\" name=\"SIGNALP\" type=\"model\" />
                                <database number=\"13\" name=\"TMHMM\" type=\"model\" />
                        </databases>
                </parameters>
        </Header>";
        
            print OUTFILE "<interpro_matches>\n";
            print OUTFILE $content;
            print OUTFILE "</interpro_matches>\n";
            print OUTFILE "</EBIInterProScanResults>";
            close OUTFILE;            
        }

       
    }
        
}
############################ End split_interproscan() ###############################

1;
