#!/usr/bin/perl -w

########################################################
#
# Common functions
# Author: Sim Ngak Leng
# First created: 2011-05-16
# Last modified: 2011-05-16
#
########################################################

use strict;

sub readMeta() {
    my ($metafile) = @_;
    my %results = ();
    open(META, "<$metafile") || die "Unable to open $metafile\n";
    while(my $line = <META>) {
	chomp $line;
	if ($line ne "" && $line !~ /^#/) {
	    my ($key, $value) = split("=", $line);
	    $results{$key} = $value;
	}
    } #end while
    close(META);
    return \%results;
} #end readMeta


sub trim {
    my ($str) = @_;
    if ($str ne "") {
	$str =~ s/^\s+//;
	$str =~ s/\s+$//;
    }
    return $str;
}

sub deldir {

    my ($dirtodel) = @_;
    print "Deleting $dirtodel\n";
    my $sep = '/';
    opendir(DIR, $dirtodel);
    my @files = readdir(DIR);
    closedir(DIR);
    
    @files = grep { !/^\.{1,2}/ } @files;
    @files = map { $_ = "$dirtodel$sep$_"} @files;
    @files = map { (-d $_)?deldir($_):unlink($_) } @files;
    
    rmdir($dirtodel);
}


1;
