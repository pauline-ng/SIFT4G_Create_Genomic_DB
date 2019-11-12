#!/usr/bin/perl -w

########################################################
# 
# Script that takes SIFTprediction files
# and map back to records
# 
# Author: Sim Ngak Leng
# First created: 2011-07-08
# Last modified: 2011-10-19
# 
########################################################

use strict;
#require 'common-utils.pl';
use File::Basename;
use Cwd qw(abs_path);
my $directory_of_script = dirname(abs_path(__FILE__));
require $directory_of_script . '/common-utils.pl';


if (scalar @ARGV != 2) {

    die "Usage: perl $0 <metafile> <chr>  \n" .
	"Example: perl $0 human.txt chr14\n" ; 
}

my $TRUE = 0;
my $FALSE = 1;
my $DEBUG = $FALSE;

my ($meta_file, $chr_of_interest ) = @ARGV;
my $meta_href = readMeta($meta_file);
my %meta_hash = %{$meta_href};

my $sift_prediction_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SIFT_SCORE_DIR"};
my $single_records_file = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_BY_CHR_DIR"} . "/" 
	 .  $chr_of_interest . ".singleRecords";
my  $output_dir =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"};
my $output_file =  $chr_of_interest  . "_scores.Srecords";
#unless ($output_file =~ /^chr/) {
#	$output_file = "chr" . $output_file;
#}
#print  $output_file . "\n";
 
####################################################################################
#
# ALGORITHM
# Get directory where SIFTprediction files are, eg. hg19_siftscores_chr3
#  parse and get gene_id (72455), aa1, aapos2, aa2, score, median if TOLERATED/DELETERIOUS
#  store in dictionary <gene_id:aa1:aapos2:aa2, score:median>
# Get singleRecord file, for each line, get score/median
# Replace score/median part
#
# For each SIFTprediction file
#   if ! (/Warning/ || /Not Scored/) 
#      if (TOLERATED || DELETERIOUS)   I777K   DELETERIOUS     0.00    4.32    2       10
#        get aa1, aapos2, aa2, score, median, others
#        match chr, aa1, aapos2, aa2 against record in hg19_chr3.singleRecords
#        put in score, median
#
####################################################################################


#if ($DEBUG == $TRUE) { $sift_prediction_dir = "/home/simnl/projects/elvimo_data/human/testSiftScoresUpdate"; }


my $siftfiles_aref = &getSiftPredictionFileNames($sift_prediction_dir, $single_records_file); 
my $siftPredictions_href = &getSiftPredictions($sift_prediction_dir, $siftfiles_aref);


&matchAndInsertScores($siftPredictions_href, $single_records_file, $output_dir . "/" . $output_file);

exit (0);
#print "SIFT Prediction files are in : $sift_prediction_dir\n";
#print "Single Records files are in : $single_records_dir\n";


### SUB ROUTINES #############################################
sub writeOutToFile() {
    print "Writing out to file for $chr_of_interest.\n";
    my $data_aref = $_[0];
    my $outputdir = $_[1];
    my $outfilename = $_[2];
    if (! -d $outputdir) {
	print "$outputdir does not exist, creating.\n"; # DO NOT RE-CREATE BECAUSE WE ARE GOING TO DO THIS CHR BY CHR
	make_dir($outputdir);
    }
    
    my @final_data = @{$data_aref};
    my $num_rec = scalar(@final_data);
    print "There are $num_rec records to write out.\n";

    #my $outputfile = $outputdir . "/" . $chr_of_interest . ".record";
    my $outputfile = $outputdir . "/" . $outfilename;
    print "Writing out to $outputfile\n";
    open(OUT, ">$outputfile") || die "Unable to open $outputfile\n";
    foreach my $data (@final_data) {
	print OUT "$data\n";
    }
    close(OUT);
    print "Completed writing out to $outputfile.\n";

} #end writeOutToFile



sub matchAndInsertScores() {

    my $dum_cnt = 0;

    print "Matching and adding sift scores if available.\n";
    my ($sift_scores_href, $single_records_file, $outfile) = @_;
    #my %scores_hash = %{$sift_scores_href};
#    my @single_rec_array = @{$single_records_aref};

#    my @results = (); # At the end of the day, @results == @single_rec_array in terms of size
    my $remainder = ",,,,,,,,,,,";
    my $cnt = 0;
    my $rowsWithScores = 0;
    my $rowsNoScores = 0;
    open (IN_RECORDS, $single_records_file) || die "can't open $single_records_file";
    open (OUT_RECORDS, ">$outfile") || die "can't open $outfile";
    print "output will be in " . $outfile . "\n";
    my $record;
    while ($record = <IN_RECORDS>) {
#    foreach my $record (@single_rec_array) {
	chomp ($record);
	if ($cnt < 5) {	print"$record\n"; $cnt++; } # for debugging and testing.

	# Check
	my @d = split(",", $record);
	my $num = scalar(@d);


#chr1,1247493,1247494,-1,novel,,60924,,CDS,Synonymous,T,C,1337,1338,CCA,CCg,P,P,445,446,1,,,,,,,,,,1,1,1,1,0,1,618,10,25,0,207,860,T:0.188976377953,C:0.811023622047,T:0.982517482517,C:0.0174825174825,T:0.932432432432,C:0.0675675675676,ref:0.00,alt:0.00,T:0.572314049587,C:0.427685950413,T:0.606946983547,C:0.393053016453
	my ($chr, $coord1, $coord2, $orn, $rsid, $ensg, $uniq_key_id, $ensp,
	    $region, $snp, $nt1, $nt2, $ntpos1, $ntpos2, $codon1, $codon2, 
	    $aa1, $aa2, $aapos1, $aapos2, $aa1_valid, $enst_valid, $dummy_score, $dummy_median, $dummy_seqs_rep, @others) = split(",", $record);


	# Pauline removed chr and $aapos1 to try and make key smaller
	# as chr1 keeps getting killed
	my $aa1_pos_aa2 =  $aapos2 . $aa2;
	my $key = "$uniq_key_id:$aa1_pos_aa2";
#	print "in sift records $key\n";
#	if (defined($scores_hash{$key})) {
#	    my $value = $scores_hash{$key};
	if (defined ($sift_scores_href->{$key})) {
	    my  $value = $sift_scores_href->{$key}; 
	    my ($score, $median, $seqs_rep) = split(/:/, $value);
	    # The 1 is for enst_valid which is never used.

	    if ($aa1_valid && $aa1_valid eq "") { $aa1_valid = 1; } # we just assume valid
	    if (!defined($aa1_valid)) { $aa1_valid = 1; }

	    my @withScores = ($chr, $coord1, $coord2, $orn, $rsid, $ensg, $uniq_key_id, $ensp,
			      $region, $snp, $nt1, $nt2, $ntpos1, $ntpos2, $codon1, $codon2, 
			      $aa1, $aa2, $aapos1, $aapos2, $aa1_valid, 1, $score, $median, $seqs_rep, @others); # <<< Added first num here, to confirm with Pauline.
	    my $scoredRecord = join(",", @withScores);

	    $scoredRecord .= $remainder;
	    #push @results, $scoredRecord;
	    print OUT_RECORDS $scoredRecord . "\n";
	    $rowsWithScores++;
	} else {
	    print OUT_RECORDS $record . "\n";
	#push @results, $record; # no scores, just store as per normal
	    $rowsNoScores++;
	}
    } #end foreach
    
    print "Rows with Scores: $rowsWithScores\nRows without Scores: $rowsNoScores\n";

    close (IN_RECORDS);
    close (OUT_RECORDS);
} #end matchAndInsertScores



sub getSiftPredictions() {
    print "Storing SIFT predictions.\n";
    my $sift_prediction_dir = $_[0];
    my $listOfSiftPredictionFiles = $_[1];
    my @sift_pred_array = @{$listOfSiftPredictionFiles};

    my %results = ();

    my $valid = 0;
    my $not_scored = 0;
    my $invalid = 0;
    foreach my $predFile (@sift_pred_array) {
	# just in case periods are in gene id name
	my @gene_id_fields = split (/\./, $predFile);
	my $gene_id = join (".", @gene_id_fields[0..$#gene_id_fields-1]);
#	my ($gene_id, $d) = split(/\./, $predFile);
#	print "chromosome: $chromo\tGeneID: $gene_id\t\n";

	my $siftfile =  $sift_prediction_dir . "/" . $predFile;
	print $siftfile . "\n";
	open(INFILE, "<$siftfile") || die "Unable to open $siftfile\n";
	while(my $siftline = <INFILE>) {
	    chomp $siftline;
	    if ($siftline =~ /TOLERATED/ || $siftline =~ /DELETERIOUS/) {
		$valid++;
		my ($aa1_pos_aa2, $tolerance, 
		    $score, $median, $num1, $num2) = split("\t", $siftline);

		if (defined($score) && defined($median) && defined($num1)) {
		    my $pos_aa2 = substr ($aa1_pos_aa2, 1,100); 
		   # Pauline remove first amino acid because trying to make the key smaller
		    my $key = "$gene_id:$pos_aa2";
#		    my $key = "$chromo:$gene_id:$aa1_pos_aa2";
		    my $value = "$score:$median:$num1";
#		print "in store predictions $key\t$score\t$median\n";
		    $results{$key} = $value;
		}
	    } elsif ($siftline =~ /NOT SCORED/) {
		$not_scored++;
	    } else {
		$invalid++;
	    }
	} #end while
	close(INFILE);
    } #end foreach

    print "Valid: $valid\tInvalid: $invalid\tNot Scored: $not_scored.\n";

    return \%results; # <k,v> = <gene_id:aa1:aapos2:aa2, score:median>
} #end getSiftPredictions


sub getSiftPredictionFileNames() {
    print "Getting SIFTprediction gene names, might take a while\n";
    print "as the directory contains lots of files.\n";

    my ($pred_dir,  $singleRecordsFile) = @_; 
    my @results = ();
    my @ensts = `cat $singleRecordsFile | cut -f7 -d"," | uniq`;
    foreach my $enst (@ensts) {
	chomp ($enst);
	my $prediction_file =  $pred_dir . "/" . $enst . ".SIFTprediction"; 
	if (-e $prediction_file) {
#	    print "PRED: $prediction_file\n";
	    push (@results, $enst . ".SIFTprediction");
	    if ($DEBUG == $TRUE) { print "$prediction_file\n"; }
	}#end if
    } #end while
    my $num = scalar(@results);
    print "There are $num SIFTprediction files for $chr_of_interest.\n";
    return \@results;
} #end getSiftPredictionFileNames




__END__

