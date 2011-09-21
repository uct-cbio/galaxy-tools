#! /usr/bin/perl -w
package Utils::LogUtils;
use strict;
use warnings;

######################## Log info message ###########################################
sub log_info() {
    my ($fd, $message) = @_;
    print $fd "Info: " . $message . "\n";

}
######################## End log_info() #############################################


######################## Log failed message #########################################
sub log_fail() {
    my ($fd, $message) = @_;
    print $fd "Fail: " . $message . "\n";

}
######################## End log_fail() #############################################


######################## Log error message ##########################################
sub log_error() {
    my ($fd, $message) = @_;
#    $message = "Error: " . $message . "\n";
     print $fd "Error: " . $message . "\n";
#    print $fd $message;

}
######################## End log_error() #############################################

1;