#!/usr/bin/perl -w

###############################################################
# Script to create the following for a given chromosome
# (a) FASTA files
# (b) SUBST files
#
# Author: Sim Ngak Leng
# First created: 2011-06-16
# Last modified: 2013-10-04
#
# Used BioPerl to extract nt sequences to generate mRNA
# Write directly out instead of storing the data and calling writeOutFasta
#
#  Sept 16, 2016 - added back to check that sequence is divisible
#  by 3. This happens to a small number of sequences (8/28,000)
#  where Ensembl has not put a value for reading frame, and creates
#  the wrong sequence
#
###############################################################

use strict;
use DBI;
use Class::Struct;
use Bio::DB::Fasta;
require 'common-utils.pl';
require 'dna_protein_subs.pl';
#use lib "/mnt1/scripts_ensembl";
#require 'DNA_PROT.pm';
#use DNA_PROT qw (chr_is_mito);

my $TRUE  = 0;
my $FALSE = 1;

if ( scalar @ARGV != 1 ) {
	die "Usage: perl $0 <metafile>\n" . "Example: perl $0 metadocs/mouse.txt\n";
}
my ($metafile) = @ARGV;
my $meta_href  = readMeta($metafile);
my %meta_hash  = %{$meta_href};

# GET DATABASE INFORMATION

my $region_in_gene = "CDS";

my $records_table_name = "RECORDS_TABLE";

my %codon_to_amino = (
	'TTT' => 'F',
	'TCT' => 'S',
	'TAT' => 'Y',
	'TGT' => 'C',
	'TTC' => 'F',
	'TCC' => 'S',
	'TAC' => 'Y',
	'TGC' => 'C',
	'TTA' => 'L',
	'TCA' => 'S',
	'TTG' => 'L',
	'TCG' => 'S',
	'TGG' => 'W',
	'CTT' => 'L',
	'CCT' => 'P',
	'CAT' => 'H',
	'CGT' => 'R',
	'CTC' => 'L',
	'CCC' => 'P',
	'CAC' => 'H',
	'CGC' => 'R',
	'CTA' => 'L',
	'CCA' => 'P',
	'CAA' => 'Q',
	'CGA' => 'R',
	'CTG' => 'L',
	'CCG' => 'P',
	'CAG' => 'Q',
	'CGG' => 'R',
	'ATT' => 'I',
	'ACT' => 'T',
	'AAT' => 'N',
	'AGT' => 'S',
	'ATC' => 'I',
	'ACC' => 'T',
	'AAC' => 'N',
	'AGC' => 'S',
	'ATA' => 'I',
	'ACA' => 'T',
	'AAA' => 'K',
	'AGA' => 'R',
	'ATG' => 'M',
	'ACG' => 'T',
	'AAG' => 'K',
	'AGG' => 'R',
	'GTT' => 'V',
	'GCT' => 'A',
	'GAT' => 'D',
	'GGT' => 'G',
	'GTC' => 'V',
	'GCC' => 'A',
	'GAC' => 'D',
	'GGC' => 'G',
	'GTA' => 'V',
	'GCA' => 'A',
	'GAA' => 'E',
	'GGA' => 'G',
	'GTG' => 'V',
	'GCG' => 'A',
	'GAG' => 'E',
	'GGG' => 'G',
	'TAA' => '*',
	'TAG' => '*',
	'TGA' => '*'
);

struct GeneTranscript => {
	id_num => '$',
	chr    => '$',
	orn    => '$',
	txS    => '$',
	txE    => '$',
	cdsS   => '$',
	cdsE   => '$',
	exonS  => '$',
	exonE  => '$',
};

struct SingleRecord => {
	chr      => '$',
	orn      => '$',
	coord1   => '$',
	coord2   => '$',
	identity => '$',
	nt1      => '$',
	nt2      => '$',
	ntpos1   => '$',
	ntpos2   => '$',
	snp      => '$',
	aa1      => '$',
	aa2      => '$',
	aapos1   => '$',
	aapos2   => '$',
	codon1   => '$',
	codon2   => '$',
	exon_num => '$',    # EXON.1, etc.
};

# Set up logging

my $logfile = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"FASTA_LOG"};
open( LOG, ">$logfile" ) || die "Unable to open $logfile for writing.\n";
my $invalid_log = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"INVALID_LOG"};
open( INVALID, ">$invalid_log" ) || die "Unable to open $invalid_log for writing.\n";

my $nCase_log = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"ZLOGFILE"};
open( ZLOG, ">$nCase_log" ) || die "Unable to open $nCase_log for writing.\n";

###### MAIN PROGRAM

print "Completed collecting fasta sequences, safe to run generate UTR script now.\n";

#DEBUG

# Create output directory if it does not exist
my $fasta_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"FASTA_DIR"};
if ( !-d $fasta_dir ) { mkdir( $fasta_dir, 0755 ); }

my $chr_dir   = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"CHR_DOWNLOAD_DEST"};
my $subst_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SUBST_DIR"};
my $gene_transcripts_file = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"GENE_DOWNLOAD_DEST"} . "/protein_coding_genes.txt";

if ( !-d $subst_dir ) { system("mkdir -p $subst_dir"); }
&generateOutput( $meta_href, $chr_dir, $gene_transcripts_file, $fasta_dir, $subst_dir );

print "Completed creation of FASTA files: $fasta_dir\n";

print "Completed SUBST files: $subst_dir\n";

close(ZLOG);
close(INVALID);
close(LOG);

#close(MTDEBUG);    #DEBUG

#####################################################################################################

sub writeOutFasta() {
	print "Writing out fasta.\n";
	my ( $meta_href, $gene_transcripts_href, $data_href ) = @_;
	my %meta_hash = %{$meta_href};
	my $fasta_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"FASTA_DIR"};
	my $model_org = $meta_hash{"ORG"};

	# Must not delete directory
	#    if (-d $fasta_dir) {
	#	print "Deleting $fasta_dir directory.\n";
	#	#deldir($fasta_dir);
	#    }
	#    print "Creating $fasta_dir directory.\n";
	#
	if ( !-d $fasta_dir ) {
		mkdir( $fasta_dir, 0755 );
	}

	my %gt_hash = %{$gene_transcripts_href};    # <k,v> = chr:orn:coord1:... => uniq_key_id
	my %data = %{$data_href};       #
	my @keys = keys(%data);
	foreach my $key (@keys) {
		# key is chr12:1:124086695:124092191:2:124086695,124092146,:124086803,124092191,
		my $uniq_key_id = $gt_hash{$key};

		my ( $chromosome, @others ) = split( ":", $key );

		if ( !defined($uniq_key_id) ) {
			die "writeOutFasta: Unique Key does not exist for $key.\n";
		}
		my $aa_seq = $data{$key};

		#print "ID: $uniq_key_id, Key for FASTA: $key\n";
		#	    print "BEFORE CHOP: $aa_seq\n";
		if ( ( $aa_seq =~ /.*\*$/ ) || ( $aa_seq =~ /.*X$/ ) ) {
			chop($aa_seq);    # chop off the last * or X
		}    #end if ($fasta_seq =~ /.*\*$/)

		#	    print "AFTER CHOP: $aa_seq\n";

		# Dec 4, 2013 Pauline changed header of fasta file to match filename
		# to be compatible with Robert's code
		my $fasta_header_title =  "ORG_" . $chromosome . "_" . $uniq_key_id . "_" . $model_org;
		my $fasta_fpath = $fasta_dir . "/" . $fasta_header_title . ".fasta";
		open( FASTA, ">$fasta_fpath" ) || die "Unable to open $fasta_fpath for writing.\n";
		#	print FASTA ">$model_org" . "_" . "$uniq_key_id\n";
		print FASTA ">$fasta_header_title\n";
		print FASTA "$aa_seq\n\n";
		#print ">$model_org" . "_" . "$uniq_key_id\n";
		#print "$aa_seq\n";
		close(FASTA);

	}    #end foreach my $key (@keys)
	print "Completed writing out fasta files to $fasta_dir\n";
}    #end writeOutFasta

sub writeOutSubstFiles() {
	print "Writing out subst.\n";
	my ( $meta_href, $gt_href, $data_href ) = @_;
	my %data = %{$data_href};    # hash of key --> \@single_records

	my %gt_hash   = %{$gt_href};
	my %meta_hash = %{$meta_href};

	my $subst_dir = $meta_hash{"PARENT_DIR"} . "/subst";

	if ( !-d $subst_dir )
	{    # do not re-create, since we are doing chromosome by chromosome
		print "Creating $subst_dir directory.\n";
		mkdir( $subst_dir, 0755 );
	}

	my $model_org = $meta_hash{"ORG"};

	my @keys    = keys(%data);
	my $numKeys = scalar(@keys);
	for ( my $k = 0 ; $k < $numKeys ; $k++ ) {    # each key is chr:orn:...
		my $key = $keys[$k];

		my ( $chromosome, @others ) = split( ":", $key );

		my $uniq_key_id = $gt_hash{$key};
		if ( !defined($uniq_key_id) ) {
			die "writeOutSubstFiles: unique Key does not exist for $key.\n";
		}

		my $single_records_aref   = $data{$key};
		my @arr_of_single_rec_arr = @{$single_records_aref
		  };    # array of arrays @(@array_of_single_records_in_same_aapos)

		# TESTING
		#	my $num1 = scalar(@arr_of_single_rec_arr);
		#	print "Number of arrays in arr_of_single_rec_arr: $num1\n";

		# Prepare transposed array
		my @transposed_arr = ();
		for ( my $i = 0 ; $i < 9 ; $i++ ) {
			my @arr = ();
			$transposed_arr[$i] = \@arr;
		}       #end for (my $i = 0; $i < 9; $i++)

		# TESTING
		#	my $num2 = scalar(@transposed_arr);
		#	print "Number of elements in the transposed array: $num2\n";

		my $len_of_peptide = scalar(@arr_of_single_rec_arr);

		for ( my $i = 0 ; $i < $len_of_peptide ; $i++ ) {

			my $single_records_aref = $arr_of_single_rec_arr[$i];
			my @single_records_arr  = @{$single_records_aref};    # This is the 9 or less mutations for a given position

			my $n = scalar(@single_records_arr);
			for ( my $j = 0 ; $j < $n ; $j++ ) {
				my $single_record = $single_records_arr[$j];
				if ( defined($single_record) ) {
					my $id     = $single_record->identity;
					my $aa1    = $single_record->aa1;
					my $aa2    = $single_record->aa2;
					my $aapos2 = $single_record->aapos2;

					if ( !defined($aa1) || !defined($aa2) || $aa1 eq "Z" || $aa2 eq "Z" ){
						my $znt1    = $single_record->nt1;
						my $znt2    = $single_record->nt2;
						my $zcodon1 = $single_record->codon1;
						my $zcodon2 = $single_record->codon2;
						my $zid     = $single_record->identity;
						my $zcoord1 = $single_record->coord1;
						my $zcoord2 = $single_record->coord2;
						my $zntpos1 = $single_record->ntpos1;
						my $zntpos2 = $single_record->ntpos2;
						print ZLOG "$zid:$zcoord1:$zcoord2:$znt1:$znt2:$zntpos1:$zntpos2:$zcodon1:$zcodon2\n";
					}

					# For substitution files, if it's U or X, we cannot have it as well as they represent stops.
					if ( $aa1 eq "U" || $aa1 eq "X" ) { $aa1 = "*"; }
					if ( $aa2 eq "U" || $aa2 eq "X" ) { $aa2 = "*"; }

   					# TEST
  					#		    if ($aa1 eq "*") { print "AA1: $aapos2 <-> $len_of_peptide: $aa1\n"; }
   					#		    if ($aa2 eq "*") { print "AA2: $aapos2 <-> $len_of_peptide: $aa2\n"; }
   					# END TEST

					# Write out if meets requirements
					# 1. Cannot have the last stop codon, but stop codons within sequences should be tried out
					# 2. Cannot have aa1 == aa2 (Synonymous SNPs, makes no sense to have sift scores for that)
					if ($aa1 ne "" && $aa2 ne "" && $aa1 ne "*" && $aa2 ne "*" && $aa1 ne "Z" && $aa2 ne "Z"){
					    #$aa1 ne $aa2 && We want to compute synonymous SNP scores.
						my $substitution = $aa1 . $aapos2 . $aa2;
						my $t_aref       = $transposed_arr[$j];
						my @t_arr        = @{$t_aref};
						push @t_arr, $substitution;
						$transposed_arr[$j] = \@t_arr;
					}    # if meets requirements
				}    #end if (defined($single_record))
			}    #end for (my $j = 0; $j < $n; $j++)
		}    #end for (my $i = 0; $i < scalar(@arr_of_single_rec_arr); $i++)

		# Now to write it out to a file

		for ( my $z = 0 ; $z < scalar(@transposed_arr) ; $z++ ) {
			my $subst_fpath = $subst_dir . "/ORG_". $chromosome . "_" . $uniq_key_id . "_". $model_org . "_". $z  . "_.subst";
			open( SUBST, ">$subst_fpath" )|| die "Unable to open $subst_fpath for writing.\n";
			my $aref = $transposed_arr[$z];
			my @a    = @{$aref};
			foreach my $subst (@a) {
				print SUBST "$subst\n";
			}
			close(SUBST);
		}

	}    #end for (my $k = 0; $k < scalar(@keys); $k++)

	print "Completed writing out to $subst_dir.\n";
}    #end writeOutSubstFiles

sub generateOutput() {
	print "Generating Single Records... takes a long time\n";
	my ( $meta_href, $dbdir, $gene_transcripts_file, $fasta_dir, $subst_dir ) =
	  @_;
	print "about to prepare $dbdir\n";

	my $chr_db = Bio::DB::Fasta->new("$dbdir", -reindex=>1) || die "Died: $!\n";
#	my $chr_db = Bio::DB::Fasta->new("$dbdir") || die "Died: $!\n";
	print "done preparing $dbdir\n";
	my $meta_hash = %{$meta_href};

	my $peptide_check_file = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"PEPTIDE_LOG"};
	my $cds_file = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"FASTA_DIR"} . "/cds.fasta";

#	my $ncbi_genetic_code = 1; #standard genetic code . Pauline

	open( PEPTIDE, ">$peptide_check_file" ) || die "Unable to open $peptide_check_file for writing.\n"; # this contains only those with ENS, for checking with ensPep
	my $num_transcripts = 0;
	my $num_of_correct_peptide_seqs = 0; # For checking, I want to know how many transcripts translate to Mxxxx* sequences
	open( IN_TX, $gene_transcripts_file ) || die "can't open $gene_transcripts_file";
	my $tx_line;

	# At this point, I have fasta sequence
	# Take gene transcript and create GeneTranscript

	while ( $tx_line = <IN_TX> ) {
		chomp($tx_line);
		my ($chr,            $t_orn,     $cds_start,
			$cds_end,        $num_exons, $adj_exons_starts,
			$adj_exons_ends, $txid,      $geneid,
			$genename,       $proteinid, $transcript_name,
			$biotype
		) = split( ":", $tx_line );
	
		# get codon table particular for chr (X or nuclear)
		# Pauline
		my $ncbi_genetic_code = 1; # default to standard
		if (chr_is_mito ($chr)) {
			$ncbi_genetic_code = $meta_hash{"MITO_GENETIC_CODE_TABLE"};
		} elsif (chr_is_plastid ($chr)) {
                	$ncbi_genetic_code = $meta_hash{"PLASTID_GENETIC_CODE_TABLE"};
	        } else {
			$ncbi_genetic_code = $meta_hash{"GENETIC_CODE_TABLE"};
		}
#		print "ncbi genetic code is $ncbi_genetic_code\n";
		my %codon_to_amino_gencode_adjusted =
			adjust_to_nonstandard_code ($ncbi_genetic_code, %codon_to_amino);

		# Let's start with exons only first.
		my $nt_sequence = "";    # This should be multiple of 3.

		#	print "transcript type $transcript_type\n";

		# incomplete_cdsstart this didn't work
		#	if ($transcript_type eq "INCOMPLETE_CDSSTART") {
		#		$nt_sequence = "NNN";
		#	}
		my @adj_starts = split( ",", $adj_exons_starts );
		my @adj_ends   = split( ",", $adj_exons_ends );

		my $num_adj_start_exons = scalar(@adj_starts);

		#	print "NUMBER OF EXONS AFTER ADJUSTMENT: $num_adj_start_exons.\n";

		for ( my $i = 0 ; $i < scalar(@adj_starts) ; $i++ ) {

			my $exon_start = $adj_starts[$i];    # This gives me the correct nucleotide
			my $exon_end = $adj_ends[$i] - 1;

			#pauline trying to fix coordinates
			$exon_start += 1;
			$exon_end   += 1;

		  	#my $chrstring = "Homo_sapiens.GRCh37.74.dna.chromosome." . $chr;
		  	#print "trying to retrieve $txid $chr:$exon_start,$exon_end\n";
			my $fasta_subseq = $chr_db->seq("$chr:$exon_start,$exon_end");

		    #my $fasta_subseq = $chr_db->seq($chr,$exon_start => $exon_end);
			$nt_sequence .= $fasta_subseq;

#			print $nt_sequence . "\n";

		}    #end for loop

		# Store associated coordinates, needed for coord1, coord2 in generating SingleRecords
		my @array_of_coordinates = ();    # Need to store nt coordinates corresponding to the base.
		for ( my $i = 0 ; $i < scalar(@adj_starts) ; $i++ ) {
			my $exon_start = $adj_starts[$i];    # This gives me the correct nucleotide
			my $exon_end = $adj_ends[$i];
			for ( my $pos = $exon_start ; $pos < $exon_end ; $pos++ ) {
				push @array_of_coordinates, $pos;
			}
		}   #end for loop

		# Now store exon number
		my @exon_numbers_arr = ();    # For region: EXON.1, EXON.2, etc.
		my $currExonNumber   = 1;
		if ( $t_orn == -1 ) {
			$currExonNumber = scalar(@adj_starts);
		}    # Exon numbers are flipped if transcribed in negative direction.
		for ( my $i = 0 ; $i < scalar(@adj_starts) ; $i++ ) {
			my $exon_start     = $adj_starts[$i];
			my $exon_end       = $adj_ends[$i];
			my $length_of_exon = $exon_end - $exon_start;

			# Store exon number
			for ( my $j = 0 ; $j < $length_of_exon ; $j++ ) {
				push @exon_numbers_arr, $currExonNumber;
			}
			if ( $t_orn == -1 ) {
				$currExonNumber--;
			} # If transcribed in negative direction, last exon_start/end is the first exon.
			else { $currExonNumber++; }
		}
		my $num_in_exon_numbers_arr = scalar(@exon_numbers_arr);

		#print "Number of exons numbers: $num_in_exon_numbers_arr\n";

		# CHECKING, the number of coordinates stored must be of the same size as the length of the nt sequence.
		# Checked 2011-06-22: The implementation is correct. I got same number of coordinates and nucleotide length
		my $num_of_coords = scalar(@array_of_coordinates);
		my $nt_len        = length($nt_sequence);
#		print MTDEBUG ">$txid DNA SEQUENCE (length: $nt_len):\n$nt_sequence\n";
		print ">$txid DNA SEQUENCE (length: $nt_len):\n$nt_sequence\n";
		if ( $nt_len < 20 ) {
			print "short sequence  $txid $nt_sequence\n";
		}
		if ( $num_of_coords != $nt_len ) {
			print "Coordinates array has $num_of_coords, but nucleotide sequence is $nt_len. Should be the same.\n";
			print $txid . " " . $nt_sequence . "\n";
		}

		#else {
			# print "Correct number of coordinates stored for associated nucleotide sequence.\n";
		#}

		my $length_of_nt_seq = length($nt_sequence);

		# print "Length of nt sequence for protein = $length_of_nt_seq\n";
		my $remainder = $length_of_nt_seq % 3;

		# Sept 15, 2016. Pauline added back check that length is 3
		# as ensembl_gene_format_to_ucsc.pl should have already fixed
		# this
		if ($remainder != 0) {
			print "Not multiple of 3, remainder: $remainder\n";
			print INVALID ">$txid\n$nt_sequence\n\n";
		} else {

		# Create Amino Acid Sequence
		my $aa_seq = "";
		my @nts = split( "", $nt_sequence );
		my $number_of_bases = scalar(@nts);
		foreach ( my $i = 0 ; $i < $number_of_bases ; $i += 3 ) {
			# SWARNA EDIT (if loop) 07-sept-2015
			# if (i + remainder) = length  then we just ignore it
			if (($i + $remainder ) == $number_of_bases){
				next;
			}
			my $nt1 = $nts[$i];
			my $nt2 = $nts[ $i + 1 ];
			my $nt3 = $nts[ $i + 2 ];
			my $cdn = $nt1 . $nt2 . $nt3;
			if ( $t_orn == -1 ) {
				  my $first  = &complement_base($nt3);
				  my $second = &complement_base($nt2);
				  my $third  = &complement_base($nt1);
				  $cdn = $first . $second . $third;
			}
			my $codon = uc($cdn);
			#Pauline
			my $aa    = $codon_to_amino_gencode_adjusted{$codon};
  			# If it is a stop codon but NOT the last one, change to X for SIFT predictions
			my $coord_of_third_base = $i + 4;
			if (defined($aa) && $aa eq "*" && ( $coord_of_third_base < $number_of_bases ) ){   
				  # Then this is NOT the last codon
				  $aa = "X";
			}
			if (!defined($aa) || $aa eq "Z" ) {
				  $aa_seq .= "X";    # To check with Pauline
			}
			else {
				  $aa_seq .= $aa;
			}
		}
		if ( $t_orn == -1 ) {
			  $aa_seq = reverse($aa_seq);
		}

		#print "$aa_seq\n"; # This is to be written out into fasta file
		#print LOG "$aa_seq\n";

		if ( $aa_seq =~ /^M.*\*/ ) {
			  $num_of_correct_peptide_seqs++;
		}

		print PEPTIDE "$proteinid\n$aa_seq\n";

		my $aaseqlength = length($aa_seq);
#		print MTDEBUG ">$txid $proteinid AA SEQUENCE (length: $aaseqlength):\n$aa_seq\n\n";

		# 2011-06-27: Don't chop off first, do that only when your're writing out fasta sequence
		# Now, if aa_seq ends with stop (*), chop it off	
		# if ($aa_seq =~ /.*\*$/) { chop($aa_seq); }

	 	# At this point, we have what we require to start generating single records
	 	# Now create single records for each position.
	 	# Inputs:
	 	# $id (gene id)
	 	# $aa_seq (for aa1)
	 	# $chr (chr)
	 	# $t_orn (for deciding flip)
	 	# @array_of_coordinates (ntpos1, ntpos2)
	 	# @nts (nt1 and generation of nt2, codon1, codon2)
	 	# @exon_numbers_arr
	 	#
	 	# Outputs:
		# Array of Single Records
	 	# Array of array of Subst Files

		### GETTING AN IDEA OF WHAT WAS PASSED TO createSingleRecords
		#	    print "id: $id\nchr: $chr\norn: $t_orn\n";
		#	    print "AA Sequence: $aa_seq\n";
		#	    if (length($aa_seq) > 5) {
		#		for(my $i = 0; $i < 5; $i++) {
		#		    print "Coord at [$i]: $array_of_coordinates[$i]\n";
		#		    print "NTS at [$i]: $nts[$i]\n";
		#		    print "Exon number at [$i]: $exon_numbers_arr[$i]\n";
		#		}
		#	    }
		# Pauline Jan 2016 adjust_to_nonstandard_code
		my ($single_records_aref) = &singleRecordsMaker( $txid, $chr, $t_orn, $aa_seq, \@array_of_coordinates, \@nts, \@exon_numbers_arr, \%codon_to_amino_gencode_adjusted );

		my $subst_fpath = $subst_dir . "/" . $txid . ".subst";
		open( SUBST, ">$subst_fpath" ) || die "Unable to open $subst_fpath for writing subst file.\n";

		my @arr_of_single_rec_arr = @{$single_records_aref};
		my $len_of_peptide  = scalar(@arr_of_single_rec_arr);

		for ( my $i = 0 ; $i < $len_of_peptide ; $i++ ) {		
			  my $single_records_aref = $arr_of_single_rec_arr[$i];
			  my @single_records_arr  = @{$single_records_aref}; # This is the 9 or less mutations for a given position
			  my $n = scalar(@single_records_arr);
			  for ( my $j = 0 ; $j < $n ; $j++ ) {
				  my $single_record = $single_records_arr[$j];
				  if ( defined($single_record) ) {
					  my $id     = $single_record->identity;
					  my $aa1    = $single_record->aa1;
					  my $aa2    = $single_record->aa2;
					  my $aapos2 = $single_record->aapos2;

					  # For substitution files, if it's U or X, we cannot have it as well as they represent stops.
					  if ( $aa1 eq "U" || $aa1 eq "X" ) { $aa1 = "*"; }
					  if ( $aa2 eq "U" || $aa2 eq "X" ) { $aa2 = "*"; }

					  # Write out if meets requirements
					  # 1. Cannot have the last stop codon, but stop codons within sequences should be tried out
					  # 2. Cannot have aa1 == aa2 (Synonymous SNPs, makes no sense to have sift scores for that)
					  if ($aa1 ne "" && $aa2 ne ""  && $aa1 ne "*" && $aa2 ne "*" && $aa1 ne "Z" && $aa2 ne "Z"){
					  	  #$aa1 ne $aa2 && We want to compute synonymous SNP scores.
						  my $substitution = $aa1 . $aapos2 . $aa2;
						  print SUBST $substitution . "\n";
					  }    # if meets requirements
				  }    #end if (defined($single_record))
			  }    #end for (my $j = 0; $j < $n; $j++)
			  # Pauline Dec 4, 2013 - edit to write out synonymous substitutions
			  my $aa_orig = substr( $aa_seq, $i, 1 );
			  if (  $aa_orig ne "X" && $aa_orig ne "" && $aa_orig ne "*" && $aa_orig ne "Z" ) {
				  print SUBST $aa_orig . ( $i + 1 ) . $aa_orig . "\n";
			  }
		}    #end for (my $i = 0; $i < scalar(@arr_of_single_rec_arr); $i++)

		close(SUBST);

		# Dec 4, 2013 Pauline changed header of fasta file to match filename
		# to be compatible with Robert's code
		# also removed X's and *'s at end of sequence to prevent seg faults
		#            print "BEFORE CHOP: $aa_seq\n";
		if ( ( $aa_seq =~ /.*\*$/ ) || ( $aa_seq =~ /.*X$/ ) ) {
			  chop($aa_seq);    # chop off the last * or X
		}    #end if ($fasta_seq =~ /.*\*$/)

		#            print "AFTER CHOP: $aa_seq\n";

		if ( $aa_seq =~ /X/ ) {
			  print "XError:$txid\n$aa_seq\n$nt_sequence\n\n";
		}

		# Pauline. Jan 13, 2014. Throw out sequences that are all X'es
		# as that causes an error
		# replace X with A because if there are a few X's, still
		# want SIFT prediction
		# Pauline Jan 27, 2014. Must throw out sequences with < 5 amino acids
		# in length, else swsharp will fail
		my $aa_seq_length = length($aa_seq);
		my $num_X = ( $aa_seq =~ tr/X// );
		if ( $aa_seq_length > 5 && $num_X < $aa_seq_length * 0.5 ) {

			  my $aas = "ACDEFGHIKLMNPQRSTVWY";
			  while ( $aa_seq =~ /[$aas]X[$aas]/ ) {
				  $_ = $aa_seq;
				  /(.*)([$aas])X([$aas])(.*)/;
				  $aa_seq = $1 . $2 . "A" . $3 . $4;
				  print "transformed $txid to $aa_seq\n";
			  }

			  my $fasta_fpath = $fasta_dir . "/" . $txid . ".fasta";
			  open( FASTA, ">$fasta_fpath" )|| die "Unable to open $fasta_fpath\n";
			  print FASTA ">$txid\n";
			  print FASTA "$aa_seq\n\n";
			  close(FASTA);
		}    # end check not too many X'es

		} #end else (is multiple of 3)
		$num_transcripts++;
	}    # end of reading  a transcript lines

	print LOG "$num_of_correct_peptide_seqs out of $num_transcripts translated to sequences that start with M and end with *\n";
	close(PEPTIDE);
	close(IN_TX);

}    #end generateOutput

sub singleRecordsMaker() {
	  my ( $uniq_key_id, $chr, $t_orn, $aa_sequence, $coords_aref, $nts_aref,$exon_numbers_aref, $adj_codon_table_ref  ) = @_;

	  #    print "Creating Single Records for $uniq_key_id\n";
	  my @aa_seq       = split( "", $aa_sequence );
	  my @coords       = @{$coords_aref};
	  my @nts          = @{$nts_aref};
	  my @exon_numbers = @{$exon_numbers_aref};

	  if (scalar(@coords) != scalar(@nts) || scalar(@nts) != scalar(@exon_numbers) || scalar(@exon_numbers) != scalar(@coords) ) {
		  print "Incompatible sizes in arrays.\n";
	  }

	  # Not all coords will be used bec. they might end up as *
	  my $num_nts = scalar(@nts);

	  #    print "Potential no. of single records for uniq_key_id $uniq_key_id: $num_nts\n";
	  # Algorithm
	  # Iterate through the nt sequence that makes up the aa sequence
	  # For each 3 bases, create codon (in transcription direction.
	  # Mutate to get 9 nt2

	  # Very important array of arrays
	  # Array of arrays, each array in a position i has potentially 9 single records
	  # Length of @array_of_single_records_array should be the length of the peptide sequence
      # (Remember to remove any * and also the last *)
	  # Create 9 files: pop 1 from each position's array
	  my @array_of_single_records_array = ();
	  my $aa_position = -1;    # This will make the aapos1 and aapos2 zero base.
	# Pauline added because transcripts not multiples of 3 are being
	# translated. Sept. 13. 2016
	  if ( $t_orn == -1 ) { $aa_position = floor ( $num_nts / 3 ); }

	  for ( my $i = 0 ; $i < $num_nts ; $i += 3 ) {
		  my $first  = $i;
		  my $second = $i + 1;
		  my $third  = $i + 2;    # these are my ntpos2

		  # EXON NUMBER
		  my $exon_num1 = $exon_numbers[$first];
		  my $exon_num2 = $exon_numbers[$second];
		  my $exon_num3 = $exon_numbers[$third];

		  # NT1
		  my $base1 = uc( $nts[$first] );
		  my $base2 = uc( $nts[$second] );
		  my $base3 = uc( $nts[$third] );

		  # COORD1
		  my $crd1 = $coords[$first];    # these are my coord2
		  my $crd2 = $coords[$second];
		  my $crd3 = $coords[$third];

		  # I now want to mutate each base, to get 9 mutations
		  my $original_codon = $base1 . $base2 . $base3;    # still positive
		       # Codons are always positive wrt to Transcription direction.
		  if ( $t_orn == -1 ) {
			  $original_codon = reverse_dna($original_codon);
		  }

		  ###########################################
		  # 2011-10-14 Bug fix: orientation issue
		  ###########################################
		  my $mutate_base1_href;
		  my $mutate_base2_href;
		  my $mutate_base3_href;

		  if ( $t_orn == 1 ) {

			  my $b1_mut_aref = &getOtherNTs($base1);    # T -> A,C,G but if -1 A -> T,C,G
			  my $b2_mut_aref = &getOtherNTs($base2);
			  my $b3_mut_aref = &getOtherNTs($base3);

			  # This gives me cTG:tTG:gTG for $base1 = A
			  #	    print "forward $uniq_key_id  ";
			  $mutate_base1_href = &mutateBase( $first, $crd1, $base1, $b1_mut_aref,"X" . $base2 . $base3, $exon_num1 );
			  $mutate_base2_href = &mutateBase( $second, $crd2, $base2, $b2_mut_aref, $base1 . "X" . $base3, $exon_num2 );
			  $mutate_base3_href = &mutateBase( $third, $crd3, $base3, $b3_mut_aref, $base1 . $base2 . "X", $exon_num3 );
			  $aa_position++;
		  }
		  else {

			  # print "We are in negative orn\n";
			  my $rev_comp_base1 = complement_base($base3);
			  my $rev_comp_base2 = complement_base($base2);
			  my $rev_comp_base3 = complement_base($base1);

			  #	    print "$base3\t$base2\t$base1\n";
			  #	    print "$rev_comp_base1\t$rev_comp_base2\t$rev_comp_base3\n";
			  my $b1_mut_aref = &getOtherNTs( $base1, $t_orn );    # This is correct, DO NOT use rev_comp_base here.
			  my $b2_mut_aref = &getOtherNTs( $base2, $t_orn );
			  my $b3_mut_aref = &getOtherNTs( $base3, $t_orn );

			  #	print "reverse $uniq_key_id ";
			  $mutate_base1_href = &mutateBase( $first, $crd1, $base3, $b1_mut_aref,$rev_comp_base1 . $rev_comp_base2 . "X", $exon_num1 );
			  $mutate_base2_href = &mutateBase( $second, $crd2, $base2, $b2_mut_aref, $rev_comp_base1 . "X" . $rev_comp_base3, $exon_num2 );
			  $mutate_base3_href = &mutateBase( $third, $crd3, $base1, $b3_mut_aref, "X" . $rev_comp_base2 . $rev_comp_base3, $exon_num3 );
			  $aa_position--;
			  #print "AAPOS: $aa_position\n";

		  }

		  # This is for cases where there is a stop codon in the sequence that's not at the end
		  # I need to let make3SingleRecords know that if $third + 1 == $num_nts, then if AA1 ends
		  # up being a stop codon, let it be a stop codon, otherwise, set it as X
		  my $isLastCodon = $FALSE;
		  if ( ( $third + 1 ) == $num_nts ) {
			  $isLastCodon = $TRUE;
		  }
		  if ( $t_orn == -1 ) {
			  if ( $first == 0 ) {
				  $isLastCodon = $TRUE;
			  }
		  }

		  my $base1_singleRecords_aref = &make3SingleRecords( $mutate_base1_href, $aa_position, $uniq_key_id, $chr, $t_orn, $original_codon, $isLastCodon, $num_nts, $adj_codon_table_ref );
		  my $base2_singleRecords_aref = &make3SingleRecords( $mutate_base2_href, $aa_position, $uniq_key_id, $chr, $t_orn, $original_codon, $isLastCodon, $num_nts, $adj_codon_table_ref );
		  my $base3_singleRecords_aref = &make3SingleRecords( $mutate_base3_href, $aa_position, $uniq_key_id, $chr, $t_orn, $original_codon, $isLastCodon, $num_nts, $adj_codon_table_ref );
		  
		  my @a1 = @{$base1_singleRecords_aref};
		  my @a2 = @{$base2_singleRecords_aref};
		  my @a3 = @{$base3_singleRecords_aref};
		  my @all_records_for_1_aa_position = ( @a1, @a2, @a3 );    # I should have 9 mutations per amino acid
		  my $size = scalar(@all_records_for_1_aa_position);

		  # So there are potentially 9 substitutions in @all_records_for_1_aa_position
		  # None of these must be in the same subst file.
		  push @array_of_single_records_array, \@all_records_for_1_aa_position;
	  }    #end for (my $ntpos = 0; $ntpos < $num_nts; $ntpos += 3)

	  return ( \@array_of_single_records_array );

}    #end singleRecordsMaker

sub getAminoAcid() {

	  #    my ($codon, $isLastCodon) = @_;
	  my ( $codon, $isLastCodon, %codon_to_amino_gencode_adjusted ) = @_;
	  my $triplet = uc($codon);

	  if ( $triplet =~ /N/ ) {
		  return "Z"; # This means the FASTA sequence has Ns, we can't do anything about this, so we do not have the codon
	  }

	  # March 22, 2014 if triplet isn't inside, assign X
	  #    my $amino_acid = "X";
	  #    if (exists ($codon_to_amino{$triplet})) {
	  #        $amino_acid = $codon_to_amino{$triplet};
	  #     }

	  my $amino_acid = $codon_to_amino_gencode_adjusted{$triplet};
	  if ( $triplet eq "TAG" && $isLastCodon == $FALSE ) {
		  $amino_acid = "U"; # This is for database, for fasta it's "X", but that's taken care of in generateOutput
	  }
	  return $amino_acid;
}    #end getAminoAcid()

sub make3SingleRecords() {
	  my ( $href, $aa_position, $identity, $chr, $orn, $orig_codon, $isLastCodon, $nts_length, $codon_to_amino_adjusted_href ) = @_;
	  my %hash = %{$href};
	  my %codon_to_amino_adjusted = %{$codon_to_amino_adjusted_href};
	  my @keys = keys(%hash); # key == $ntpos:$coord:$original_base:$mutant_base

	  #    my $num_keys = scalar(@keys); print "Number of keys is $num_keys.\n";

	  my @result = ();
	  foreach my $key (@keys) {
		  my ( $ntpos1, $coord1, $original_base, $mutant_base, $exon_num ) = split( ":", $key );
		  if ( $orn == -1 ) {
			  my $new_ntpos1 = $nts_length - $ntpos1 - 1;
			  $ntpos1 = $new_ntpos1;
		  }

		  #	print "coord2 $coord1 is supposed to be numeric ke5 is $key \n";
		  my $ntpos2        = $ntpos1 + 1;
		  my $coord2        = $coord1 + 1;
		  my $mutated_codon = $hash{$key};
		  my $aapos1        = $aa_position;
		  my $aapos2        = $aa_position + 1;

		 # Make a singleRecord, this constitutes 1 row in the variation database
		  my $singleRecord = SingleRecord->new();
		  $singleRecord->identity($identity);
		  $singleRecord->chr($chr);
		  $singleRecord->coord1($coord1);
		  $singleRecord->coord2($coord2);
		  $singleRecord->orn($orn);
		  $singleRecord->nt1( uc($original_base) );
		  $singleRecord->nt2( uc($mutant_base) );
		  $singleRecord->ntpos1($ntpos1);
		  $singleRecord->ntpos2($ntpos2);
		  $singleRecord->aapos1($aapos1);
		  $singleRecord->aapos2($aapos2);

		  if ( $orn == -1 ) {
			  #TAa ==> TAt
			  $mutated_codon =~ tr/acgt/tgca/;
		  }

		  my $orig_aa    = &getAminoAcid( $orig_codon,    $isLastCodon, %codon_to_amino_adjusted );
		  my $mutated_aa = &getAminoAcid( $mutated_codon, $isLastCodon, %codon_to_amino_adjusted );
		  $singleRecord->aa1($orig_aa); # At this point, if stop codon, we get *
		  $singleRecord->aa2($mutated_aa);
		  if ( $orig_aa eq $mutated_aa ) { $singleRecord->snp("Synonymous"); }
		  else { $singleRecord->snp("Nonsynonymous"); }

		  $singleRecord->codon1($orig_codon);
		  $singleRecord->codon2($mutated_codon);

		  $singleRecord->exon_num($exon_num);

		  # For the moment, we just store everything and filter out when we start writing out substitution files.
		  push @result, $singleRecord;

	  }    # foreach my $key (@keys)

	  return \@result;

}    #end make3SingleRecords

sub mutateBase() {
	  my ( $ntpos, $coord, $original_base, $aref, $codon, $exon_num ) = @_;
	  my %result    = ();
	  my @mutations = @{$aref};
	  foreach my $mutant_base (@mutations) {
		  my $mutant_codon = $codon;
		  $mutant_codon =~ s/X/$mutant_base/;
		  my $key = "$ntpos:$coord:$original_base:$mutant_base:$exon_num";

		  #	print "$key\n";
		  $result{$key} = $mutant_codon;
	  }    #end foreach my $mutant (@mutations)
	  return \%result;
}    #end mutateBase

sub getOtherNTs() {
	  my ($base) = @_;
	  my @result = ();
	  if ( $base =~ /A/ ) { @result = ( "c", "g", "t" ); }
	  elsif ( $base =~ /C/ ) { @result = ( "a", "g", "t" ); }
	  elsif ( $base =~ /G/ ) { @result = ( "a", "c", "t" ); }
	  elsif ( $base =~ /T/ ) { @result = ( "a", "c", "g" ); }
	  return \@result;
}    #end getOtherNTs

sub is_amino_acid {
	  my ($aa) = @_;
	  if ( $aa =~ /[ACDEFGHIKLMNPQRSTVWY]/ ) {
		  return 1;
	  }
	  else {
		  return 0;
	  }
}

__END__




