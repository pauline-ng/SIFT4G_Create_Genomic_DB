#!/usr/bin/perl -w

#############################################################
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
my %meta_hash = %{$meta_href};
my @chromosomes = getChr ($meta_hash{"PARENT_DIR"} . "/". $meta_hash{"GENE_DOWNLOAD_DEST"});
my $mo_out_dir = $meta_hash{"PARENT_DIR"} . "/". $meta_hash{"ORG_VERSION"};
chomp ($mo_out_dir);

#print join ("," , @chromosomes);
get_genes_in_db ($mo_out_dir,  @chromosomes);

sub count_genes_in_SIFT_db_file
{
	my ($file) = @_;

	my %genes = ();
	my %genes_with_sift_score = ();

	my $cds_pos = 0;
	my $cds_pos_with_sift_score = 0;
	my $sift_confident_score = 0;

	my $line;
	open (DB_IN, $file) || die "can't open $file";
	while ($line = <DB_IN>) {
		if ($line =~ /CDS/) {
			my @fields = split ('\t', $line);
			my $txid = $fields[3];
			my $siftscore = $fields[10]; 
			my $medianscore = $fields[11];
			$genes{$txid} = 1;
			if ($siftscore ne "NA" && $siftscore >= 0.0 && $siftscore ne "") {
				$genes_with_sift_score{$txid} = 1;
			}
			if ($fields[7] ne $fields[8] && $fields[7] ne "*" 
			&& $fields[8] ne "*") {
				$cds_pos++;
				if ($siftscore ne "NA" && $siftscore >= 0.0 
					&& $siftscore ne "") {

					$cds_pos_with_sift_score++;
					if ($medianscore <= 3.5) {
						$sift_confident_score++;
					}
				}	
			}
		}
	}
	close (DB_IN);
	my $num_genes =  scalar keys %genes;
	my $num_genes_sift = scalar keys %genes_with_sift_score;
	return ($num_genes, $num_genes_sift, $cds_pos, 
				$cds_pos_with_sift_score, $sift_confident_score);
}

sub get_genes_in_db
{
	my ($mo_out_dir, @chromosomes) = @_;
	my %count_hash;
	my $results_file = $mo_out_dir . "/CHECK_GENES.LOG";
	open (OUT_RES, ">$results_file") || die "can't write to $results_file";
	print OUT_RES "Chr\tGenes with SIFT Scores\tPos with SIFT scores\tPos with Confident Scores\n";
	my $tmp_file = $meta_hash{"PARENT_DIR"}  . "/SIFT_tmp233_rm.txt";
	if (-e $tmp_file) {
		system ("rm $tmp_file");
	}

	my ($total_cds_gene_count_with_sift,$total_cds_gene_count, 
		$total_cds_pos_with_sift_scores, 
		$total_cds_pos, $total_sift_confident_score)
				= (0) x 5; 

	for (my $i = 0; $i < @chromosomes; $i++) {
		my $chr = $chromosomes[$i]; 
		chomp ($chr);
		my $file = $mo_out_dir . "/" . $chr . ".gz";
		system ("zcat $file > $tmp_file"); 
		my ($cds_gene_count, $cds_gene_count_with_sift, $cds_pos,
			$cds_pos_with_sift_scores, $sift_confident_score) = 
			count_genes_in_SIFT_db_file ($tmp_file);
		my $cds_perc = 0;
		if ($cds_gene_count > 0) {
			$cds_perc = $cds_gene_count_with_sift/$cds_gene_count * 100;
		}
		my $pos_perc = 0;
		if ($cds_pos > 0) {
			$pos_perc = $cds_pos_with_sift_scores / $cds_pos * 100;
		}
		my $sift_conf_perc = 0;
		if ($cds_pos_with_sift_scores > 0) {
			$sift_conf_perc = $sift_confident_score/$cds_pos_with_sift_scores * 100;
		}	
		print OUT_RES "$chr\t" . round ($cds_perc) . 
			" ($cds_gene_count_with_sift/$cds_gene_count)\t" . round ($pos_perc) . " ($cds_pos_with_sift_scores/$cds_pos)\t" . round ($sift_conf_perc) . "($sift_confident_score/$cds_pos_with_sift_scores)\n";
		system ("rm $tmp_file");

		$total_cds_gene_count_with_sift += $cds_gene_count_with_sift;
		$total_cds_gene_count += $cds_gene_count; 
                $total_cds_pos_with_sift_scores += $cds_pos_with_sift_scores;
                $total_cds_pos += $cds_pos;  
		$total_sift_confident_score += $sift_confident_score;
	}
	print OUT_RES "\n";


# do it genome-wide, get totals
	my $total_cds_perc = 0; 
	 if ($total_cds_gene_count > 0) {
             $total_cds_perc = $total_cds_gene_count_with_sift/$total_cds_gene_count * 100;
         }
         my $total_pos_perc = 0;
         if ($total_cds_pos > 0) {
             $total_pos_perc = $total_cds_pos_with_sift_scores / $total_cds_pos * 100;
         }
         my $total_sift_conf_perc = 0;
         if ($total_cds_pos_with_sift_scores > 0) {
             $total_sift_conf_perc = $total_sift_confident_score/$total_cds_pos_with_sift_scores * 100;
         }    


	print OUT_RES "ALL\t" . round ($total_cds_perc) .
                 " ($total_cds_gene_count_with_sift/$total_cds_gene_count)\t" 
	. round ($total_pos_perc) . 
	" ($total_cds_pos_with_sift_scores/$total_cds_pos)\t" . 
	round ($total_sift_conf_perc) . 
	"($total_sift_confident_score/$total_cds_pos_with_sift_scores)\n";

	close (OUT_RES);
}

sub round
{
	my ($float ) = @_;
	return (int ($float+0.5));
}

sub num_genes_per_chr 
{
	my (%meta_hash) = @_;

	my %chr_counts;
	my $gene_transcripts_file =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"GENE_DOWNLOAD_DEST"} . "/protein_coding_genes.txt";
	my @counts = `cat $gene_transcripts_file | cut -f1,8 -d: | sort | uniq | cut -f1 -d: | uniq -c`;
	for (my $i= 0; $i < @counts; $i++) {
		chomp ($counts[$i]);
		$counts[$i] =~ s/^\s+//;
		my ($num, $chr) = split (/\s+/, $counts[$i]);
#		print "chr $chr num $num\n";
		$chr_counts{$chr} = $num;
	}
	return (%chr_counts);
}

