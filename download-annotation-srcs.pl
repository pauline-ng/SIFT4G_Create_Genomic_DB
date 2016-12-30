#!/usr/bin/perl -w

########################################################
#
# Download ensGene/refGene/knownGene/ccdsGene data
# Author: Sim Ngak Leng
# First created: 2011-05-16
# Last modified: 2011-05-28
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

my $download_dest = $meta_hash{"PARENT_DIR"} . "/" .  $meta_hash{"GENE_DOWNLOAD_DEST"};

if (! -d $download_dest) {
    mkdir($download_dest, 0755);
}

my $gene_file = $meta_hash{"GENE_DOWNLOAD_SITE"}; 
my @files_to_download;
push (@files_to_download, $gene_file); 

## If there is ensPep.txt.gz download that as well, this is for checking
my $oth_files_to_download = $meta_hash{"PEP_FILE"};
$oth_files_to_download =~ s/^\s+//;
$oth_files_to_download =~ s/\s+$//;
if ($oth_files_to_download ne "") {
    my @otherFiles = split(",", $oth_files_to_download);
    foreach my $f (@otherFiles) {
	push @files_to_download, $f;
    }
}


foreach my $file (@files_to_download) {
    print "Downloading $file \n";
    my $command = "wget --directory-prefix=$download_dest $file";
    my $result = system($command);
    if ($result == 0) { print "Downloading of $file successful.\n"; }
    else { print "Unsuccessful in downloading $file\n"; }
}
print "Completed downloading of files to $download_dest\n";

__END__
