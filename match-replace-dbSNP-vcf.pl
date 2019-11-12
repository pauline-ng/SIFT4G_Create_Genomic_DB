#!/usr/bin/perl -w

##############################################
#
# Script to match dbSNP data
# with singleRecords
#
# Author: Sim Ngak Leng
# First created: 2012-01-31
# Last modified: 2012-02-20
#
# 2012-02-15: Fixed issue with missing rsids
# 2012-02-20: Fixed missing obsolete2 when writing out
##############################################


use strict;
#require 'common-utils.pl';
use File::Basename;
use Cwd qw(abs_path);
my $directory_of_script = dirname(abs_path(__FILE__));
require $directory_of_script . '/common-utils.pl';


if (scalar(@ARGV) != 2) {
    die "Usage: perl $0 <metafile> <chr> \n" .
	"Example: perl $0 human.txt 22 \n";
}

my ($metafile, $chr_of_interest) = @ARGV;
my $meta_href = readMeta($metafile);
my %meta_hash = %{$meta_href};


	my $rootfile = $chr_of_interest . "_scores.Srecords";
#	unless ($rootfile =~ /^chr/) {
#		$rootfile = "chr" . $rootfile;
#	}
	my $oldfile = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"} . "/" . $rootfile; 
	my $noncodfile = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_BY_CHR_DIR"} . "/" . $chr_of_interest . ".singleRecords_noncoding";
	my $new_noncod_file = $noncodfile . ".with_dbSNPid";

	
	my $dbSNPFile = "";
	if (-e $meta_hash{"DBSNP_DIR"}) {
		$dbSNPFile = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"DBSNP_DIR"} . "/" . "vcf_chr_" . $chr_of_interest . ".vcf.gz";
	}
	my $newfile = $oldfile . ".with_dbSNPid";
	my $errfile = $newfile . ".ERR";
	my $logfile = $newfile . ".LOG";

	my $dbsnp_href ;
	if (-e $dbSNPFile && $dbSNPFile ne "") {
		$dbsnp_href = &getDBSNPData($dbSNPFile, $errfile);
	} else {
	 	# if dbSNP file does not exist
		my %hash; # create empty hash
		$dbsnp_href = \%hash;
	} 
		
	&matchAndWriteToFile($dbsnp_href, $oldfile, $newfile, $errfile, $logfile);
	&matchAndWriteToFile($dbsnp_href, $noncodfile, $new_noncod_file, $errfile, $logfile);

exit (0);

sub matchAndWriteToFile() {
    my ($dbsnp_href, $oldfile, $newfile, $errfile, $logfile) = @_;
    my %dbsnp_map = %{$dbsnp_href};
    print "Matching and writing to $newfile\n";    
    open(LOG, ">$logfile") || die "Unable to write to $logfile\n";
    open(ERR, ">>$errfile") || die "Unable to append to $errfile\n";
    open(OUT, ">$newfile") || die "Unable to write to $newfile\n";
    open(IN, "<$oldfile") || die "Unable to read from $oldfile\n";
    while(my $row = <IN>) {
	chomp $row;
	my @fields = split (/,/, $row);
	my $chr = $fields[0]; my $coord2 = $fields[2]; my $nt1 = $fields[10];
	my $nt2 = $fields[11];
#	my ($chr,$coord1,$coord2,$orn,$rsid,$obs1,$uniq,$obs2,$region,$snp,$nt1,$nt2,$nt1pos,$nt2pos,
#	    $codon1,$codon2,$aa1,$aa2,$aa1pos,$aa2pos,$aa1_valid,$enst_valid,$score,$median,$seqs_rep, @others) = split(/,/, $row);
	my $chromosome = $chr;
	$chromosome =~ s/chr//; # Gina's dbSNP files do not have chrY, just Y

	my $key = "$chromosome:$coord2:$nt1:$nt2";

	my $rsid = "";
	if (defined($dbsnp_map{$key})) {
#	    $rsid =  $dbsnp_map{$key} . ":" . $nt2;
	    $rsid = $dbsnp_map{$key}; 
	    print LOG "$key -> $rsid -> $row\n"; 
	} else {
	    if ($nt1 eq $nt2) {
		$rsid = "ref";
	    } else {
		$rsid = "novel";
	    }
	}
	$fields[4] = $rsid;
#	$row = join (",", @fields);
#	my $new_row = "$chr,$coord1,$coord2,$orn,$rsid,,$uniq,,$region,$snp,$nt1,$nt2,$nt1pos,$nt2pos," .
#	    "$codon1,$codon2,$aa1,$aa2,$aa1pos,$aa2pos,$aa1_valid,$enst_valid,$score,$median,$seqs_rep," . join (",", @others);
	my $new_row = format_line (@fields);
	print OUT "$new_row";
	
    } #end while
    close(IN);
    close(OUT);
    close(ERR);
    close(LOG);
}


sub format_line {
	my @fields = @_;
	my $chr = $fields[0];
        my $pos = $fields[2];
        my $transcript_id = $fields[6];
        my $rsid = $fields[4];
   	my $region_type = $fields[8]; 
        my $variant_type = $fields[9];
        my $ref_allele = $fields[10];
        my $new_allele = $fields[11];
        my $aa_pos = $fields[19];
        my $aa_old = $fields[16];
        my $aa_new = $fields[17];
	my $gene_id = $fields[28];
	my $gene_name = $fields[29];
        my $sift_score = "";
        if (exists ($fields[22])) {
                $sift_score = $fields[22];
        }
        my $sift_confidence = "";
        if (exists ($fields[23])) {
                $sift_confidence = $fields[23];
        }
        my $num_seqs = "";
        if (exists ($fields[24])) {
                $num_seqs = $fields[24];
        }

        # assign variant type numbers with 0 being most damaging
        # and higher numbers being tolerated
        my $variant_type_num = "";
        if ($aa_new eq "*" && $aa_old ne "*") { # this is introducing a stop
                $variant_type_num = 0;
        } elsif ($aa_old ne $aa_new) {
                $variant_type_num = 1; #nonsynomyous
        } elsif ($aa_old eq $aa_new && $aa_old ne "" && $aa_new ne "") {
                $variant_type_num = 2; #synonymous
        } elsif ($aa_old eq "" && $aa_new eq "") {
		$variant_type_num = 3; # noncoding
	}
        my $name1 = "";
        my $name2 = "";
#       print "$outfile unique gene id is $uniq_gene_id sift $sift_confidence\n";
	$transcript_id = check($transcript_id);
	$gene_id = check ($gene_id);
	$gene_name = check ($gene_name);
	$region_type = check ($region_type);
	$aa_old = check ($aa_old); 
	$aa_new = check ($aa_new);
	$aa_pos = check ($aa_pos);
	$sift_score = check ($sift_score);
	$sift_confidence = check ($sift_confidence); 
	$num_seqs = check ($num_seqs);
	$rsid = check ($rsid);

        my $newline = "$chr\t$pos\t$ref_allele\t$new_allele\t$transcript_id\t$gene_id\t$gene_name\t$region_type\t$aa_old\t$aa_new\t$aa_pos\t$sift_score\t$sift_confidence\t$num_seqs\t$rsid\t$variant_type_num\n";
	return $newline;
}

sub check 
{
	my ($val) = @_;
	if (!defined $val) {
		return "NA";
	}
	if ($val eq "") {
		return "NA";
	} else {
		return $val;
	}
}

sub getDBSNPData() {
    my ($dbSNPFile, $errfile) = @_;


    my %dbSNP_data_hash = ();

    open(ERRORFILE, ">$errfile") || die "Unable to open $errfile for writing.\n";
    print "opening up $dbSNPFile";
    if ($dbSNPFile =~ /\.gz$/) {
	open (DBSNPFILE, "gunzip -c $dbSNPFile |") || die "can't open pipe to $dbSNPFile";
    } else {
	open(DBSNPFILE, "<$dbSNPFile") || die "Unable to open $dbSNPFile for reading\n";
   }

    while(my $dbSNP_data = <DBSNPFILE>) {
	chomp $dbSNP_data;
#	my ($dbsnp_rsid, $dbsnp_chr, $dbsnp_coord, $dbsnp_snp2chrOrn, 
#	    $dbsnp_alleles, $dbsnp_ancestral, $dbsnp_contigID, 
#	    $dbsnp_contig2chrOrn, $dbsnp_contigAllele) = split(":", $dbSNP_data);
	my ($dbsnp_chr, $dbsnp_coord, $dbsnp_rsid, $ref_allele, $dbsnp_new) = split ("\t", $dbSNP_data); 	
	my @alt_alleles = split (',', $dbsnp_new);
	
	    # Ignore indels and multple alleles (GA/CG)

	    foreach my $alt_allele (@alt_alleles) {
		my $key = "$dbsnp_chr:$dbsnp_coord:$ref_allele:$alt_allele";

#		if (length($alt_allele) > 1 || length($ref_allele) > 1) {
#			    print "allele has more than 1 base: $key\n";
#		}
#		if ($key =~ /::/) {
#		    print "Missing elements: $key\n";
#		}
		$dbSNP_data_hash{$key} = $dbsnp_rsid;			
	    }
    } # end while 
    close(DBSNPFILE);
    close(ERRORFILE);

    return \%dbSNP_data_hash;

}




__END__



