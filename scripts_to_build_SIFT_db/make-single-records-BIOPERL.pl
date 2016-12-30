#!/usr/bin/perl -w

###############################################################
# Script to create single row records for final database
#
# Author: Sim Ngak Leng
# First created: 2011-06-16
# Last modified: 2013-10-04
#
# Used BioPerl to generate nt sequence, and
# write directly out to file instead of storing
# and calling writeOutputToFilesByChromosomes
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

our $TRUE = 0;
our $FALSE = 1;

our %codon_to_amino = ('TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C', 'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
              'TTA' => 'L', 'TCA' => 'S', 'TTG' => 'L', 'TCG' => 'S', 'TGG' => 'W', 'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H',
              'CGT' => 'R', 'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R', 'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q',
              'CGA' => 'R', 'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R', 'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N',
              'AGT' => 'S', 'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S', 'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K',
              'AGA' => 'R', 'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R', 'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D',
              'GGT' => 'G', 'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G', 'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E',
              'GGA' => 'G', 'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G', 'TAA' => '*', 'TAG' => '*', 'TGA' => '*');

struct GeneTranscript => {
    id_num => '$',
    chr => '$',
    orn => '$',
    txS => '$',
    txE => '$',
    cdsS => '$',
    cdsE => '$',
    exonS => '$',
    exonE => '$',
};

struct SingleRecord => {
    chr => '$',
    orn => '$',
    coord1 => '$',
    coord2 => '$',
    identity => '$',
    nt1 => '$',
    nt2 => '$',
    ntpos1 => '$',
    ntpos2 => '$',
    snp => '$',
    aa1 => '$',
    aa2 => '$',
    aapos1 => '$',
    aapos2 => '$',
    codon1 => '$',
    codon2 => '$',
    exon_num => '$', # EXON.1, etc.
};

my $DEBUG = $FALSE; # Change to FALSE for production

if (scalar @ARGV != 1) {
    die "Usage: perl $0 <metafile> \n" .
        "Example: perl $0 metadocs/mouse.txt\n";
}
my ($metafile) = @ARGV;

my $meta_href = readMeta($metafile);
my %meta_hash = %{$meta_href};


my $outputdir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_BY_CHR_DIR"};
# Check and create directory if it does not exist. DO NOT RE-CREATE IF IT EXISTS!
if (! -d $outputdir) {
    print "Output directory does not exist, it will be created.\n";
    make_dir($outputdir);
} else {
        #clean out files
	if (glob ("$outputdir/*.singleRecords")) {
        	system ("rm $outputdir/*.singleRecords");
	}
}

my $region_in_gene = "CDS";
my  $gene_transcripts_file = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"GENE_DOWNLOAD_DEST"} . "/protein_coding_genes.txt";
###### MAIN PROGRAM

my $chr_dir =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"CHR_DOWNLOAD_DEST"};
&generateOutput($meta_href, $chr_dir, $gene_transcripts_file, $outputdir);


#####################################################################################################
sub collectFastaSequence() {
    my ($meta_href, $chr_of_interest) = @_;
    print "Getting fasta sequence for $chr_of_interest.\n";
    my %meta_hash = %{$meta_href};

    my $chr_dir = $meta_hash{"PARENT_DIR"} . "/". $meta_hash{"CHR_DOWNLOAD_DEST"};
    my $tmp_dir = $meta_hash{"PARENT_DIR"} . "/tmp";
    if (! -d $tmp_dir) { make_dir($tmp_dir); }
    my $chr_file = $chr_dir . "/" . $chr_of_interest . ".fa.gz";
    my $tmp_file = $meta_hash{"PARENT_DIR"} . "/tmp/" . $chr_of_interest . "_fasta.tmp";
    my $cmd = "zcat $chr_file > $tmp_file";
    my $result = system($cmd);

    my $sequence = "";
    if ($result == 0) {
        open(FA, "<$tmp_file") || die "Unable to open $tmp_file\n";
        while(my $line = <FA>) {
            chomp $line;
            if ($line !~ /^>/) {
                $sequence .= $line;
            }
        } #end while
        close(FA);
    }
    unlink($tmp_file);
    return $sequence;
} #end collectFastaSequence


sub generateOutput() {
    print "Generating Single Records... takes a long time\n";
#    my ($meta_href, $chr_seq, $gene_transcripts_href) = @_;
    my ($meta_href, $dbdir, $gene_transcripts_file,  $outputdir) = @_;

    open (IN_TX, $gene_transcripts_file) || die "can't open $gene_transcripts_file";
    my $tx_line;
#only happens once, and is necessary if download is incomplete an fails.
# notice this assumes we run make-single-records-BIOPERL.pl before the other scripts
    my $chr_db = Bio::DB::Fasta->new("$dbdir", -reindex=>1) || die "Died: $! indexing $dbdir\n";
#    my $chr_db = Bio::DB::Fasta->new("$dbdir") || die "Died: $!\n";

    while ($tx_line = <IN_TX>) {
	chomp ($tx_line);
        my ($chr, $t_orn, $cds_start, $cds_end, $num_exons, $adj_exons_starts, $adj_exons_ends, $txid, $geneid, $genename, $proteinid, $transcript_name, $biotype) = split(":", $tx_line);
        my $chr_file = $outputdir . "/" . $chr . ".singleRecords";
	my $chr_file_proteins = $outputdir."/".$chr.".singleRecords_proteins.fa";
	my $chr_invalid = $outputdir. "/". $chr.".invalidRecords";

    #    print "$outputdir\n";
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
#        print "ncbi genetic code is $ncbi_genetic_code\n";
        my %codon_to_amino_gencode_adjusted =
                     adjust_to_nonstandard_code ($ncbi_genetic_code, %codon_to_amino);




        open(FILE, ">>$chr_file") || die "Unable to write out to $chr_file\n";
	open(INVALID,">>$chr_invalid") || die "Unable to write out to $chr_invalid\n";
	open(FILEPROT, ">>$chr_file_proteins") || die "Unable to write out to $chr_file_proteins\n";
        if ($biotype eq "protein_coding") {
                write_out_for_protein_coding ($tx_line, \%codon_to_amino_gencode_adjusted, $chr_db);
        } else {
		print "something wrong, wrong input file? biotype $biotype for $tx_line\n";
	} 
        close (FILE);
   } # end while reading gene file
} # end generateOutput


sub write_out_for_protein_coding  {
        my ($tx_line, $codon_to_aa_hash_ref, $chr_db) = @_;
	my %adjusted_codon_to_amino = %{$codon_to_aa_hash_ref};
	chomp ($tx_line);
        my ($chr, $t_orn, $cds_start, $cds_end, $num_exons, $adj_exons_starts, $adj_exons_ends, $txid, $geneid, $genename, $proteinid, $transcript_name, $biotype) = split(":", $tx_line);

        my $nt_sequence = ""; # This should be multiple of 3.

        my @adj_starts = split(",",$adj_exons_starts);
        my @adj_ends = split(",",$adj_exons_ends);
        my $num_adj_start_exons = scalar(@adj_starts);

        for (my $i = 0; $i < scalar(@adj_starts); $i++) {
            my $exon_start = $adj_starts[$i]; # This gives me the correct nucleotide
            my $exon_end = $adj_ends[$i] - 1;

#                print "trying to access $chr $exon_start $exon_end\n";
            #pauline trying to fix coordinates
            $exon_start += 1;
            $exon_end += 1;
            my $fasta_subseq = $chr_db->seq("$chr:$exon_start,$exon_end");
#		print "seq $fasta_subseq\n";
            $nt_sequence .= $fasta_subseq;
        } #end for loop

        # Store associated coordinates, needed for coord1, coord2 in generating SingleRecords
        my @array_of_coordinates = (); # Need to store nt coordinates corresponding to the base.
        for (my $i = 0; $i < scalar(@adj_starts); $i++) {
            my $exon_start = $adj_starts[$i]; # This gives me the correct nucleotide
            my $exon_end = $adj_ends[$i];
            for (my $pos = $exon_start; $pos < $exon_end; $pos++) {
                push @array_of_coordinates, $pos;
            }
        } #end for loop

        # Now store exon number
        my @exon_numbers_arr = (); # For region: EXON.1, EXON.2, etc.
        my $currExonNumber = 1;
        if ($t_orn == -1) { $currExonNumber = scalar(@adj_starts); } # Exon numbers are flipped if transcribed in negative direction.
        for (my $i = 0; $i < scalar(@adj_starts); $i++) {
            my $exon_start = $adj_starts[$i];
            my $exon_end = $adj_ends[$i];
            my $length_of_exon = $exon_end - $exon_start;
            # Store exon number
            for (my $j = 0; $j < $length_of_exon; $j++) {
                push @exon_numbers_arr, $currExonNumber;
            }
            if ($t_orn == -1) { $currExonNumber--; } # If transcribed in negative direction, last exon_start/end is the first exon.
            else { $currExonNumber++; }
        }
        my $num_in_exon_numbers_arr = scalar(@exon_numbers_arr);

        # CHECKING, the number of coordinates stored must be of the same size as the length of the nt sequence.
        # Checked 2011-06-22: The implementation is correct. I got same number of coordinates and nucleotide length
        my $num_of_coords = scalar(@array_of_coordinates);
        my $nt_len = length($nt_sequence);
        if ($num_of_coords != $nt_len) {
            print "Coordinates array has $num_of_coords, but nucleotide sequence is $nt_len. Should be the same.\n";
        }

        my $length_of_nt_seq = length($nt_sequence);
        my $remainder = $length_of_nt_seq % 3;
        if ($remainder != 0) {
            #print "Not multiple of 3, remainder: $remainder\n";
         #   print INVALID ">$geneid"."_"."$txid"."_"."$proteinid\n$nt_sequence\n"; #<<< THIS IS ALREADY DONE IN THE PREVIOUS STEP WHEN FASTA/SUBST FILES ARE CREATED.
        } else {
            # Create Amino Acid Sequence
            my $aa_seq = "";
            my @nts = split("",$nt_sequence);
            my $number_of_bases = scalar(@nts);
            foreach (my $i = 0; $i < $number_of_bases; $i += 3) {
                my $nt1 = $nts[$i];
                my $nt2 = $nts[$i+1];
                my $nt3 = $nts[$i+2];
                my $cdn = $nt1 . $nt2 . $nt3;
                if ($t_orn == -1) {
                    my $first = &complement_base($nt3);
                    my $second = &complement_base($nt2);
                    my $third = &complement_base($nt1);
                    $cdn = $first . $second . $third;
                }
                my $codon = uc($cdn);
                my $aa = $adjusted_codon_to_amino{$codon};

                # If it is a stop codon but NOT the last one, change to X for SIFT predictions
                my $coord_of_third_base = $i + 4;
                if ($aa eq "*" && ($coord_of_third_base < $number_of_bases)) { # Then this is NOT the last codon
                    $aa = "X";
                }
#		if (!defined ($aa) || $aa eq "Z") {
#			$aa = "X";
#		} 
	        $aa_seq .= $aa;
        }
            if($remainder>0){
                $aa_seq .="X";
             }
	     print INVALID "$geneid\t$txid\t$proteinid\t$nt_sequence\n";

            if ($t_orn == -1) {
                $aa_seq = reverse($aa_seq);
            }

            # 2011-06-27: Don't chop off first, do that only when your're writing out fasta sequence
            # Now, if aa_seq ends with stop (*), chop it off
            #if ($aa_seq =~ /.*\*$/) { chop($aa_seq); }

	
            my ($single_records_aref) =
                &singleRecordsMaker($txid, $chr, $t_orn, $aa_seq, 
			 \@array_of_coordinates, 
			\@nts, \@exon_numbers_arr, $codon_to_aa_hash_ref);
            my @sr_array = @{$single_records_aref};
            foreach my $array_ref (@sr_array) {
                my @array_of_sr = @{$array_ref};
                foreach my $singleRecord (@array_of_sr) {
		   #print $txid . "," . $singleRecord->chr . "," . $singleRecord->coord1 . "," . $singleRecord->nt1 . "," . $singleRecord->nt2 . "\n" ;
                    my $line = $singleRecord->chr . "," .
                        $singleRecord->coord1 . "," . $singleRecord->coord2 . "," .
                        $singleRecord->orn . "," . "," . # RSID, to be filled up later
                        "," . # ENSG will be removed
                        $singleRecord->identity . "," . # Takes the place of ENST
                        "," . # ENSP -> will be obtained via query to GENE_TRN_2_PROTEIN table
                        "$region_in_gene," . # REGION
                        $singleRecord->snp . "," .
                        $singleRecord->nt1 . "," . $singleRecord->nt2 . "," .
                        $singleRecord->ntpos1 . "," . $singleRecord->ntpos2 . "," .
                        $singleRecord->codon1 . "," . $singleRecord->codon2 . "," .
                        $singleRecord->aa1 . "," . $singleRecord->aa2 . "," .
                        $singleRecord->aapos1 . "," . $singleRecord->aapos2 . "," .
                        "," . # CDS
                        "," . #AA1_VALID
                        "," . #ENST_VALID
                        "," . #SCORE, to be filled up later
                        "," . #MEDIAN
                        "," . #SEQ_REP
                         ",$txid,$geneid,$genename,$proteinid,$biotype";
                    print FILE "$line\n"; # Write out to file.

                } # end foreach single record
            } # end foreach array
	    print FILEPROT ">$geneid"."_"."$txid"."_"."$proteinid\n$aa_seq\n";
        } #end else (is multiple of 3)

    close(FILE);

} #end subroutine


sub singleRecordsMaker() {
    my ($uniq_key_id, $chr, $t_orn, $aa_sequence, $coords_aref, $nts_aref, $exon_numbers_aref, $adj_codon_table_ref) = @_;
    my @aa_seq = split("", $aa_sequence);
    my @coords = @{$coords_aref};
    my @nts = @{$nts_aref};
    my @exon_numbers = @{$exon_numbers_aref};

    if (scalar(@coords) != scalar(@nts)         ||
        scalar(@nts) != scalar(@exon_numbers)   ||
        scalar(@exon_numbers) != scalar(@coords)) {
        print "Incompatible sizes in arrays.\n";
    }

    # Not all coords will be used bec. they might end up as *
    my $num_nts = scalar(@nts);

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

    my $aa_position = -1; # This will make the aapos1 and aapos2 zero base.


    if ($t_orn == -1) { $aa_position = ($num_nts / 3); }

	# March 22, 2014  Add 1 for the incomplete codon in the beginning
	# April 5, 2014 this doesn't work and was deleted
#	print "transcript type $transcript_type\n";
#    if ($transcript_type eq "INCOMPLETE_CDSSTART") {
#         $aa_position += 1; 
#	print "entered in here $aa_position\n";
#    }

    for (my $i = 0; $i < $num_nts; $i += 3) {
        my $first = $i; my $second = $i+1; my $third = $i+2; # these are my ntpos2

        # EXON NUMBER
        my $exon_num1 = $exon_numbers[$first];
        my $exon_num2 = $exon_numbers[$second];
        my $exon_num3 = $exon_numbers[$third];

        # NT1
        my $base1 = uc($nts[$first]); # Positive wrt to reference
        my $base2 = uc($nts[$second]);
        my $base3 = uc($nts[$third]);

        # COORD1
        my $crd1 = $coords[$first]; # these are my coord2
        my $crd2 = $coords[$second];
        my $crd3 = $coords[$third];


        # I now want to mutate each base, to get 9 mutations
        my $original_codon = $base1 . $base2 . $base3;
        # Codons are always positive wrt to Transcription direction.
        if ($t_orn == -1) {
            $original_codon = reverse_dna($original_codon);
        }


        ###########################################
        # 2011-10-14 Bug fix: orientation issue
        ###########################################
        my $mutate_base1_href;
        my $mutate_base2_href;
        my $mutate_base3_href;

        if ($t_orn == 1) {

            my $b1_mut_aref = &getOtherNTs($base1, $t_orn); # T -> A,C,G but if -1 A -> T,C,G
            my $b2_mut_aref = &getOtherNTs($base2, $t_orn);
            my $b3_mut_aref = &getOtherNTs($base3, $t_orn);

            # This gives me cTG:tTG:gTG for $base1 = A
            $mutate_base1_href = &mutateBase($first, $crd1, $base1, $b1_mut_aref, "X".$base2.$base3, $exon_num1);
            $mutate_base2_href = &mutateBase($second,$crd2, $base2, $b2_mut_aref, $base1."X".$base3, $exon_num2);
            $mutate_base3_href = &mutateBase($third, $crd3, $base3, $b3_mut_aref, $base1.$base2."X", $exon_num3);

            $aa_position++;
        } else {
            # print "We are in negative orn\n";
            my $rev_comp_base1 = complement_base($base3);
            my $rev_comp_base2 = complement_base($base2);
            my $rev_comp_base3 = complement_base($base1);

#           print "$base3\t$base2\t$base1\n";
#           print "$rev_comp_base1\t$rev_comp_base2\t$rev_comp_base3\n";
            my $b1_mut_aref = &getOtherNTs($base1, $t_orn); # This is correct, DO NOT use rev_comp_base here.
            my $b2_mut_aref = &getOtherNTs($base2, $t_orn);
            my $b3_mut_aref = &getOtherNTs($base3, $t_orn);

            $mutate_base1_href = &mutateBase($first, $crd1, $base1, $b1_mut_aref, $rev_comp_base1.$rev_comp_base2."X", $exon_num1); # ATg, ATc, ATa
            $mutate_base2_href = &mutateBase($second,$crd2, $base2, $b2_mut_aref, $rev_comp_base1."X".$rev_comp_base3, $exon_num2);
            $mutate_base3_href = &mutateBase($third, $crd3, $base3, $b3_mut_aref, "X".$rev_comp_base2.$rev_comp_base3, $exon_num3);


            $aa_position--;
            #print "AAPOS: $aa_position\n";
        }

        # This is for cases where there is a stop codon in the sequence that's not at the end
        # I need to let make3SingleRecords know that if $third + 1 == $num_nts, then if AA1 ends
        # up being a stop codon, let it be a stop codon, otherwise, set it as X
        my $isLastCodon = $FALSE;

        if (($third + 1) == $num_nts) {
            $isLastCodon = $TRUE;
        }

        if ($t_orn == -1) {
            if ($first == 0) {
                $isLastCodon = $TRUE;
            }
        }

        my $base1_singleRecords_aref = &make3SingleRecords($mutate_base1_href, $aa_position, $uniq_key_id,  # chr:orig -> array of mutations CTc,Cta,
                                                           $chr, $t_orn, $original_codon, $isLastCodon, $num_nts, $adj_codon_table_ref);
        my $base2_singleRecords_aref = &make3SingleRecords($mutate_base2_href, $aa_position, $uniq_key_id,
                                                           $chr, $t_orn, $original_codon, $isLastCodon, $num_nts, $adj_codon_table_ref);
        my $base3_singleRecords_aref = &make3SingleRecords($mutate_base3_href, $aa_position, $uniq_key_id,
                                                           $chr, $t_orn, $original_codon, $isLastCodon, $num_nts, $adj_codon_table_ref);

        my @a1 = @{$base1_singleRecords_aref};
        my @a2 = @{$base2_singleRecords_aref};
        my @a3 = @{$base3_singleRecords_aref};
        my @all_records_for_1_aa_position = (@a1, @a2, @a3); # I should have 9 mutations per amino acid
        my $size = scalar(@all_records_for_1_aa_position);


        # So there are potentially 9 substitutions in @all_records_for_1_aa_position
        # None of these must be in the same subst file.
        push @array_of_single_records_array, \@all_records_for_1_aa_position;

    } #end for (my $ntpos = 0; $ntpos < $num_nts; $ntpos += 3)

    return (\@array_of_single_records_array);

} #end singleRecordsMaker


sub getAminoAcid() {
#    my ($codon, $isLastCodon) = @_;
    my ($codon, $isLastCodon, $type, %codon_to_amino_gencode_adjusted ) = @_;
    my $triplet = uc($codon);

# March 22, 2014 if triplet isn't inside, assign X - this doesn't work, commented out
#    my $amino_acid = "X";
#    if (exists ($codon_to_amino{$triplet})) {

	my $amino_acid = $codon_to_amino_gencode_adjusted{$triplet};

#     }
# Mar 16, 2014. For it to be a U, it has to be a reference sequence
# not a new amino acid

    if ($triplet eq "TAG" && $isLastCodon == $FALSE && $type eq "REF") {
        $amino_acid = "U"; # This is for database, for fasta it's "X", but that's taken care of in generateOutput

    }
    return $amino_acid;
} #end getAminoAcid()


sub make3SingleRecords() {
    my ($href, $aa_position, $identity, $chr, $orn, $orig_codon, $isLastCodon, $nts_length, $codon_to_amino_adjusted_href) = @_;
    my %hash = %{$href};
    my %codon_to_amino_adjusted = %{$codon_to_amino_adjusted_href};
    my @keys = keys(%hash); # key == $ntpos:$coord:$original_base:$mutant_base

    my @result = ();
    foreach my $key (@keys) {
        my ($ntpos1, $coord1, $original_base, $mutant_base, $exon_num) = split(":", $key);

        if ($orn == -1) {
            my $new_ntpos1 = $nts_length - $ntpos1 - 1;
            $ntpos1 = $new_ntpos1;
        }

        my $ntpos2 = $ntpos1 + 1;
        my $coord2 = $coord1 + 1;
        my $mutated_codon = $hash{$key};
        my $aapos1 = $aa_position;
        my $aapos2 = $aa_position + 1;

        # Make a singleRecord, this constitutes 1 row in the variation database
        my $singleRecord = SingleRecord->new();
        $singleRecord->identity($identity);
        $singleRecord->chr($chr);
        $singleRecord->coord1($coord1);
        $singleRecord->coord2($coord2);
        $singleRecord->orn($orn);
        $singleRecord->nt1(uc($original_base));
        $singleRecord->nt2(uc($mutant_base));
        $singleRecord->ntpos1($ntpos1);
        $singleRecord->ntpos2($ntpos2);
        $singleRecord->aapos1($aapos1);
        $singleRecord->aapos2($aapos2);

        if ($orn == -1) {
            $mutated_codon =~ tr/acgt/tgca/;
        }

        my $orig_aa = &getAminoAcid($orig_codon, $isLastCodon, "REF",  %codon_to_amino_adjusted);
        my $mutated_aa = &getAminoAcid($mutated_codon, $isLastCodon, "NEW",  %codon_to_amino_adjusted);
        $singleRecord->aa1($orig_aa);   # At this point, if stop codon, we get *
        $singleRecord->aa2($mutated_aa);
        if ($orig_aa eq $mutated_aa) { $singleRecord->snp("Synonymous"); }
        else { $singleRecord->snp("Nonsynonymous"); }

        $singleRecord->codon1($orig_codon);
        $singleRecord->codon2($mutated_codon);

        $singleRecord->exon_num($exon_num);

        # For the moment, we just store everything and filter out when we start writing out substitution files.
        push @result, $singleRecord;

    } # foreach my $key (@keys)

    return \@result;

} #end make3SingleRecords

sub mutateBase() {
    my ($ntpos, $coord, $original_base, $aref, $codon, $exon_num) = @_;

    my %result = ();
    my @mutations = @{$aref};
    foreach my $mutant_base (@mutations) {
        my $mutant_codon = $codon;
        $mutant_codon =~ s/X/$mutant_base/;
        my $key = "$ntpos:$coord:$original_base:$mutant_base:$exon_num";
        $result{$key} = $mutant_codon;
    } #end foreach my $mutant (@mutations)
    return \%result;
} #end mutateBase

sub getOtherNTs() {
    my ($base, $orn) = @_;
    my @result = ();
    if ($base =~ /A/) { @result = ("c","g","t","a"); }
    elsif ($base =~ /C/) { @result = ("a","g","t","c"); }
    elsif ($base =~ /G/) { @result = ("a","c","t","g"); }
    elsif ($base =~ /T/) { @result = ("a","c","g","t"); }
    return \@result;
} #end getOtherNTs



__END__

