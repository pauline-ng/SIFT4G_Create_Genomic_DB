#!/usr/bin/perl -w

####################################################
#
# Script to download required dbSNP bcp.gz files
#
# Author: Sim Ngak Leng
# First created: 2011-07-22
# Last modified: 2011-07-22
#
####################################################

use strict;

require 'common-utils.pl';

if (scalar @ARGV != 1) {
    die "Usage: perl $0 <metafile>\n" .
	"Example: perl $0 human.txt\n";
}
my ($meta_file) = @ARGV;
my $meta_href = readMeta($meta_file);
my %meta_hash = %{$meta_href};


my $organism_download_site = $meta_hash{"DBSNP_ORGANISM_DOWNLOAD_SITE"};
 
my @files_from_organism_site = ();
if (exists ($meta_hash{"DBSNP_VCF_FILE"})) {
	 @files_from_organism_site =  ($meta_hash{"DBSNP_VCF_FILE"});
} 

my $download_destination_dir;
if (exists ($meta_hash{"DBSNP_DIR"})) {
	$download_destination_dir = $meta_hash{"PARENT_DIR"} . "/". $meta_hash{"DBSNP_DIR"};
#	if (!glob ("$download_destination_dir/*.vcf.gz")) {
#		print "No dbSNP VCF found\n";
#		exit (0);
#	}
} else {
	exit (0);
}

&download_dbSNP($organism_download_site, \@files_from_organism_site, $download_destination_dir);
print "Completed.\n";

### SUB ROUTINES
sub download_dbSNP() {
    my ($src_site, $files_aref, $dest_dir) = @_;
    print "Downloading from $src_site\n";
    my @files = @{$files_aref};

    if (@files == 0) {
	print "no files listed, downloading everything\n";
	my $command = "wget -r -p -A gz --no-directories --directory-prefix=$dest_dir $src_site";
	print "$command\n";
        my $result = system($command);
        if ($result == 0) { print "Downloading of $src_site successful.\n"; }
        else { print "Unsuccessful in downloading $src_site\n"; }
    }

    # Delete away if exists
    foreach my $f (@files) {
	my $checkFile = $dest_dir . "/" . $f;
	if (-e $checkFile) {
	    print "$checkFile exists, deleting.\n";
	    unlink($checkFile);
	}
    }

    # Now download
    foreach my $file (@files) {
	my $filepath = $src_site . "/" . $file;
	print "Downloading $file\n";
	my $command = "wget --directory-prefix=$dest_dir $filepath";
	print "$command\n";
	my $result = system($command);
	if ($result == 0) { print "Downloading of $file successful.\n"; }
	else { print "Unsuccessful in downloading $file\n"; }
    } #end foreach

} #end download_dbSNP


__END__
