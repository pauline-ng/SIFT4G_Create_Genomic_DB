#!/usr/bin/perl -w

# create the SIFT databases, must be run in scripts directory because it's lookin for metadocs

use strict;
#require 'common-utils.pl'; 
#use Getopt::Std;
use File::Basename;
use Cwd qw(abs_path);
my $directory_of_script = dirname(abs_path(__FILE__));
require $directory_of_script . '/common-utils.pl';

use Getopt::Long qw (GetOptions);

my %options  = ();
#getopts ("mhf:", \%options);

my $help = "";
my $ensembl_download = "";
my $metafile = "";
GetOptions ('help|h' => \$help, 
	    'ensembl_download' => \$ensembl_download,
	    'config=s' => \$metafile );

if ($help || !$metafile || !$metafile) {
	print "make-SIFT-db-all.pl\n";
	print " -h  help\n";
	print " -config  [required config file]\n";
	print " -ensembl_download \n";
}

my $meta_href = readMeta($metafile);
my %meta_hash = %{$meta_href};

make_dir ($meta_hash{"PARENT_DIR"});
make_dir ($meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"ORG_VERSION"});

for my $key (keys %meta_hash) {
	my $value = $meta_hash{$key};
#	print $value . "\n";
	if ($key =~ /_DIR$/ && $key ne "PARENT_DIR") {
#		print "making $value\n";
		my $dir = $meta_hash{"PARENT_DIR"} . "/" . $value;
		if (! -d $dir) { mkdir($dir, 0775) or die $dir . $!; }
	}
}
my $siftscore_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"};
if (! -d $siftscore_dir) { mkdir ($siftscore_dir, 0775); } 		

if ($ensembl_download) {
	# not manual, this is from Ensembl and needs to download
	print "downloading gene annotation\n";
	`perl download-annotation-srcs.pl $metafile`; 
	print "done downloading gene annotation\n"; 

	print "downloading fasta files\n"; 
	`perl download-fasta.pl $metafile`; 
	print "done downloading DNA fasta sequences"; 

	print "download dbSNP files\n";
	`perl download-dbSNP-files.pl $metafile`; 
	if (exists ($meta_hash{"DBSNP_VCF_FILE"}) && $meta_hash{"DBSNP_VCF_FILE"} ne "") 
	{
		print "splitting dbSNP file\n";
		`perl split-dbSNP-by-chr.pl $metafile`;
	}
} # end if downloading files from Ensembl

print "converting gene format to use-able input\n";
`perl gff_gene_format_to_ucsc.pl $metafile`; 
#`perl ensembl_gene_format_to_ucsc.pl $metafile`;
print "done converting gene format\n";

# decompress chromosome files so that Bio::DB can process them
my $fasta_dir =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"CHR_DOWNLOAD_DEST"} ;
if (glob ("$fasta_dir/*.gz")) {
        $fasta_dir .= "/*.gz";
        `gunzip $fasta_dir`;
        # check that all DNA files finished downloading, otherwise do it again.
        # because all future programs won't work properly otherwise.
        if ($?) {
                die "DNA files do not exist or did not unzip properly\n";
        }
} elsif (glob ("$fasta_dir/*.fa") || glob ("$fasta_dir/*.fasta") ||  glob ("$fasta_dir/*.fna")) {
        # *.fa or *.fasta files already exist
} else {
        die "no DNA fasta files in $fasta_dir\n";
}

print "making single records file\n"; 
`perl make-single-records-BIOPERL.pl $metafile`; 
print "done making single records template\n"; 

print "making noncoding records file\n";
`perl make-single-records-noncoding.pl $metafile`;
print "done making noncoding records\n";

print "make the fasta sequences\n"; 
`perl generate-fasta-subst-files-BIOPERL.pl $metafile`; 
print "done making the fasta sequences\n";

print "start siftsharp, getting the alignments\n"; 

# remove all_prot.fasta if it already exists
my $all_prot_fasta  =  $meta_hash{"PARENT_DIR"} . "/all_prot.fasta";
if ( -e  $all_prot_fasta ) {
	unlink $all_prot_fasta;
}

my $combine_prot_fasta_command = "for file in " .  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"FASTA_DIR"} . "/*.fasta ; do cat \"\$file\" >> " . $all_prot_fasta . "; done";

`$combine_prot_fasta_command`;


my $sift4g_command = $meta_hash{"SIFT4G_PATH"} .  " -d " . $meta_hash{"PROTEIN_DB"} . " -q " . $all_prot_fasta . " --subst " .  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SUBST_DIR"} . " --out " .  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SIFT_SCORE_DIR"} . " --sub-results " ;

print $sift4g_command . "\n";

#`$siftsharp_command`;
`$sift4g_command`;

if ($? != 0) {
        exit (-1);
}

print "done getting all the scores\n";

print "populating databases\n";
# check quality
# Pauline to do , by chromosame in script

`perl make-sift-scores-db-batch.pl $metafile`;

print "checking the databases\n";
system ("perl check_genes.pl $metafile");

my @chromosomes =  getChr ($meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"GENE_DOWNLOAD_DEST"});

foreach my $chr (@chromosomes) {

        my $db_file = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"} . "/" . $chr .  "_scores.Srecords.with_dbSNPid.sorted";
        my $check_outfile = $meta_hash{"PARENT_DIR"}  .  "/" . $meta_hash{"ORG_VERSION"} . "/" . $chr .  "_SIFTDB_stats.txt";
        `python check_SIFTDB.py $db_file $check_outfile`;
}

my $pep_check_file = $meta_hash{"PARENT_DIR"} . "/ENS_peptide_check.txt";
`perl checkENSPep.pl $metafile ENS $pep_check_file`; 

my $final_outfolder = $meta_hash{"PARENT_DIR"} ."/". $meta_hash{"ORG_VERSION"};
system ("cp  $metafile $final_outfolder"); 

# make some space
my $zip_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"CHR_DOWNLOAD_DEST"} . "/*" ;
print "zipping up $zip_dir\n";
system ("gzip $zip_dir");
my $rm_dir = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"SINGLE_REC_WITH_SIFTSCORE_DIR"} . "/*";
system ("rm $rm_dir");

print "All done!\n";
exit (0);


