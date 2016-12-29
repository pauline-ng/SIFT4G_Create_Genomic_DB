#!/usr/bin/perl -w

########################################################
#
# Batch script to 
# 1) map-scores-back-to-records.pl 
# 2) match-replace-dbSNP-vcf.pl
# 3) add-gene-id-output-tsv.pl
#
#  Author: Pauline Ng
# First created: 2012-12-12
# Last modified: 2011-12-12
#
########################################################


use strict;
require 'readMeta.pl';
require 'common-utils.pl';

if (scalar @ARGV != 1) {
    die "Usage: perl $0 <metafile>\n" .
	"Example: perl $0 mouse.txt\n";
}

my ($metafile) = @ARGV;
my $meta_href = readMeta($metafile);
my %meta_hash = %{$meta_href};

my @chromosomes = getChr ( $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"GENE_DOWNLOAD_DEST"});

my $final_outfolder = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"ORG_VERSION"}  ; 
system ("mkdir $final_outfolder");

foreach my $chr (@chromosomes) {
    #my $command = "perl make-single-records-BIOPERL.pl $metafile " . "chr" . $chr;
     my $chr_of_interest =  "chr" . $chr;
	my $command1 = "perl map-scores-back-to-records.pl $metafile $chr"; 
	 print "Executing $command1\n";
	`$command1`;
     my $command2 = "perl match-replace-dbSNP-vcf.pl $metafile $chr"; 
	print "Executing $command2\n";
        `$command2`;
	my $mapped_score_file =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"} . "/" . $chr .  "_scores.Srecords"; # this is output of command1 
	my $db_file = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"} . "/" . $chr .  "_scores.Srecords.with_dbSNPid";

	my $noncod_file = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_BY_CHR_DIR"} . "/". $chr . ".singleRecords_noncoding.with_dbSNPid";
	print "db file $db_file\nnoncoding $noncod_file";
# checking file
        my $check_outfile = $final_outfolder .  "/" . $chr .  "_SIFTDB_stats.txt";
#       `python check_SIFTDB.py $db_file $check_outfile`;

	# now sort the files
	print "final outfolder is $final_outfolder\n";	
	my $db_sorted_gz =  $final_outfolder . "/"  . $chr . ".gz";
	my $db_sorted = $db_file . ".sorted";
	my $line = "\"#Position\tRef_allele\tNew_allele\tTranscript_id\tGene_id\tGene_name\tRegion\tRef_amino_acid\tNew_amino_acid\tPosition_of_amino_acid_substitution\tSIFT_score\tSIFT_median_sequence_info\tNum_seqs_at_position\tdbSNP_id\""; 
	`echo $line > $db_sorted`; 

	`cat $noncod_file $db_file  | sort -k1,1 -k2,2n -k3,3 -k4,4 -k16,16n -k12,12n | cut -f2-15  >> $db_sorted`;
	my $regions_file = $final_outfolder . "/" . $chr . ".regions";
	my $command4 = "python make_regions_file.py $db_sorted $regions_file"; 	
	print "$command4\n";
	`$command4`;
	`gzip -c $db_sorted > $db_sorted_gz`;
     print "Built single records " . $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_BY_CHR_DIR"} . "\n";
	system ("rm $mapped_score_file"); 
	system ("rm $db_file");
}








