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


