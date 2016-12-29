#!/usr/bin/perl -w

#############################################################
# 
# Script that takes dbSNP file generated from dbSNP/*.pl
# and split by chromosomes 
# This is because the dbSNP output file is large 
# (+1GB for humans, could be larger for other organisms
#
# The output file will be used in map-dbSNP-to-records.pl
# Author: Sim Ngak Leng
# First created: 2011-07-25
# Last modified: 2011-07-26
#
#############################################################

use strict;
require 'common-utils.pl';

if (scalar @ARGV != 1) {
    die "Usage: perl $0 <metafile>\n" .
	"Example: perl $0 human.txt\n";
}

my $TRUE = 0;
my $FALSE = 1;
my $DEBUG = $TRUE;

my ($meta_file) = @ARGV;
my $meta_href = readMeta($meta_file);
#my %meta_hash = %{$meta_href};

### MAIN ###################################################
&split_dbSNP_output_by_chromosomes($meta_href);

### SUB ROUTINE
sub split_dbSNP_output_by_chromosomes() {
    print "Splitting dbSNP output file into different chromosomes.\n";
    my $meta_href = $_[0];
    my %meta = %{$meta_href};
    #my $chrlist = $meta{"CHR_LIST"};
    #my @chromosomes = split(",", $chrlist);
    my $chr_fasta_dir = $meta_hash{"PARENT_DIR"} . "/". $meta{"CHR_DOWNLOAD_DEST"};
    my @chromosomes = getChr ($meta_hash{"PARENT_DIR"} . "/". $meta{"GENE_DOWNLOAD_DEST"});

    my $dbSNPFile = $meta_hash{"PARENT_DIR"} . "/".  $meta{"DBSNP_FINAL_OUTPUT_DIR"} . "/" . $meta{"DBSNP_VCF_FILE"};
    my $cat_command = "cat";
    if ($dbSNPFile =~ /\.gz$/) {
	$cat_command = "zcat"; 
    } 
    print "DBSNP main output file is : $dbSNPFile\n";
    my $outputDir = $meta_hash{"PARENT_DIR"} . "/". $meta{"DBSNP_FINAL_OUTPUT_DIR"};
    foreach my $chr (@chromosomes) {
#	print "Awk-ing chromosome chr$chr\n";
	#if ($chr eq "M") { $chr = "MT"; } # This is to cater to dbSNP calling Mitochrondria chromosome MT
	my $outfile = $outputDir . "/vcf_chr_" . $chr . ".vcf";
	my $chrField = '$1';
	my $command = "$cat_command $dbSNPFile | awk '{ if ($chrField == \"$chr\") { print }}' | gzip > $outfile.gz";
	print "$command\n";
	system($command);

    } #end foreach

} #end split_dbSNP_output_by_chromosomes

__END__


