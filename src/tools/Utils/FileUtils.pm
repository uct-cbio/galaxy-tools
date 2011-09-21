#! /usr/bin/perl -w
package Utils::FileUtils;
use strict;
use warnings;

##################### Create process directories ####################################
sub create_project_directories() {
    my ( $project_dir, $base_dir, %dirs ) = @_;

    # Create new project if one does not exist
    unless ( -e "$project_dir" ) {
        system("mkdir $project_dir");
    }

    # If previous process directory exits in project project, remove and create a new directory.
    # If no directory exists create new
    if ( -e "$base_dir" ) {
        system("rm -rf $base_dir");
        system("mkdir $base_dir");
    } else {
        system("mkdir $base_dir");
    }

    # Within process directory create additional directories for storing data and output
    foreach my $key ( keys %dirs ) {
        system("mkdir $dirs{$key}");
    }
}
######################## End find_program() #########################################

1;