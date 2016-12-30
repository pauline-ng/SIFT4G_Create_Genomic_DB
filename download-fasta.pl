#!/usr/bin/perl -w

########################################################
#
# Download chromosome fa files
# Author: Sim Ngak Leng
# First created: 2011-05-19
# Last modified: 2011-05-19
#
########################################################


use strict;
require 'common-utils.pl';

if (scalar @ARGV != 1) {
    die "Usage: perl $0 <metafile>\n" .
	"Example: perl $0 mouse.txt\n";
}

my ($metafile) = @ARGV;
my $meta_href = readMeta($metafile);
my %meta_hash = %{$meta_href};

my $download_site = $meta_hash{"CHR_DOWNLOAD_SITE"};
my $download_dest =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"CHR_DOWNLOAD_DEST"};


if (! -d $download_dest) {
    mkdir($download_dest, 0755);
}

    my $command = "wget -r -l1 --no-parent --no-directories -A \"*.dna.chromosome.*.fa.gz\" --directory-prefix=$download_dest $download_site";
    my $result = system($command);
    if ($result == 0) { print "Downloading of genome successful.\n"; }

	 $command = "wget -r -l1 --no-parent --no-directories -A \"*.dna.genome.fa.gz\" --directory-prefix=$download_dest $download_site";
         $result = system($command);

	$command = "wget -r -l1 --no-parent --no-directories -A \"*.dna.nonchromosomal.fa.gz\" --directory-prefix=$download_dest $download_site";
         $result = system($command);

print "Completed downloading of chromosome files to $download_dest\n";

#    my @files = `ls $download_dest/*.gz`;
#    foreach my $file (@files) {
#	chomp ($file);
#	print "file is $file\n";
#	print "unzipping $file\n";
#	`gunzip $file`;
#    }	

exit (0);






