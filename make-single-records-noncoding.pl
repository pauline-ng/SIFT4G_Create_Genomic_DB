#!/usr/bin/perl -w

# Pauline Ng
# Last modified: 2013-10-04
#
# Used BioPerl to generate nt sequence, and
# write directly out to file instead of storing
# and calling writeOutputToFilesByChromosomes
###############################################################

use strict;
use DBI;
use Class::Struct;
use Bio::DB::Fasta;
#require 'common-utils.pl';
#require 'dna_protein_subs.pl';
use File::Basename;
use Cwd qw(abs_path);
my $directory_of_script = dirname(abs_path(__FILE__));
require $directory_of_script . '/common-utils.pl';
require $directory_of_script . '/dna_protein_subs.pl';

our $TRUE = 0;
our $FALSE = 1;

if (scalar @ARGV != 1) {
    die "Usage: perl $0 <metafile> \n" .
        "Example: perl $0 metadocs/mouse.txt\n";
}
my ($metafile) = @ARGV;

my $meta_href = readMeta($metafile);
my %meta_hash = %{$meta_href};

my $outputdir = $meta_hash{"PARENT_DIR"} . "/". $meta_hash{"SINGLE_REC_BY_CHR_DIR"};
# Check and create directory if it does not exist. DO NOT RE-CREATE IF IT EXISTS!
if (! -d $outputdir) {
    print "Output directory does not exist, it will be created.\n";
    make_dir($outputdir);
} else {
   # remove existing files because this program appends to files
   if (glob ("$outputdir/*.singleRecords_noncoding")) {
	   system ("rm $outputdir/*.singleRecords_noncoding");
   }
} 

my $chr_dir =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"CHR_DOWNLOAD_DEST"};
my  $noncoding_transcripts_file =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"GENE_DOWNLOAD_DEST"} . "/noncoding.txt";
###### MAIN PROGRAM

&generateOutput($meta_href, $chr_dir, $noncoding_transcripts_file, $outputdir);



sub generateOutput() {
    print "Generating Single Records... takes a long time\n";
    my ($meta_href, $dbdir, $gene_transcripts_file,  $outputdir) = @_;

    open (IN_TX, $gene_transcripts_file) || die "can't open $gene_transcripts_file";
    my $tx_line;
#    my $chr_db = Bio::DB::Fasta->new("$dbdir", -reindex=>1) || die "Died: $!\n";
    my $chr_db = Bio::DB::Fasta->new("$dbdir") || die "Died: $!\n";

    while ($tx_line = <IN_TX>) {

        my ($chr, $biotype, $tx_start, $tx_end,  $geneid, $transcriptid, $genename, $transcript_name) = split("\t", $tx_line);

        # Open file for writing.
        my $chr_file = $outputdir . "/" . $chr . ".singleRecords_noncoding";
        #print "$outputdir\n";

        open(FILE, ">>$chr_file") || die "Unable to write out to $chr_file\n";

        write_out_for_noncoding ($tx_line, $chr_db);
        close (FILE);
   } # end while reading gene file
} # end generateOutput

sub write_out_for_noncoding {
        my ($tx_line, $chr_db) = @_;

        my ($chr, $biotype, $exon_start , $exon_end, $geneid, $txid, $genename, $txname) = split("\t", $tx_line);

       my $fasta_subseq = $chr_db->seq("$chr:$exon_start,$exon_end");
#	print $tx_line;

       my @nucs = ("A", "C", "G", "T");
	my $length = $exon_end - $exon_start + 1;

	if (defined($fasta_subseq)) { # some chr can't retrieve sequence{
       for (my $pos = 0; $pos < $length; $pos++) {
          foreach my $nt2 (@nucs) {
		my $genome_pos = $pos + $exon_start;
                         my $line = $chr . "," .
                            ($genome_pos- 1) . "," . ($genome_pos) . "," .
                             "," . "," . # RSID, to be filled up later
                            "," . # ENSG will be removed
                            $txid . "," . # Takes the place of ENST
                            "," . # ENSP -> will be obtained via query to GENE_TRN_2_PROTEIN table
                            "$biotype," . # REGION
                             "," .
                             substr ($fasta_subseq , $pos, 1). "," . $nt2  				. "," x 15 . 
                        ",$txid,$geneid,$genename,,$txname"; 
 		print FILE $line . "\n";;    
         } # end for each nucleotide
        } #end for $pos loop
	} else { # end couldn't retrieve sequence
		print "skipping because couldn't $tx_line retrieve sequence";
	}
}

