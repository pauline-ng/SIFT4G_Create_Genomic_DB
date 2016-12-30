#!/usr/bin/perl
use strict;
use feature qw(say);
use Switch;
use LWP::Simple;

# NOTICE Please make sure release is set for common organisms. Line 202
# (there is no current-release and has to be set.)
# #

# common, fungi, metazoa, bacteria, plants, protists 

my ($ensembl_root, $var_dir, $dbsnp_format, $collection_dir) = &set_according_to_type ("common", "");
#my ($ensembl_root, $var_dir, $dbsnp_format, $collection_dir) = &set_according_to_type ("plants", "");
#my ($ensembl_root, $var_dir, $dbsnp_format, $collection_dir) = &set_according_to_type ("fungi", "");
#my ($ensembl_root, $var_dir, $dbsnp_format, $collection_dir) = &set_according_to_type ("protists", "");
#my ($ensembl_root, $var_dir, $dbsnp_format, $collection_dir) = &set_according_to_type ("metazoa", "");
#my ($ensembl_root, $var_dir, $dbsnp_format, $collection_dir) = &set_according_to_type ("bacteria", "");
#my ($ensembl_root, $var_dir, $dbsnp_format, $collection_dir) = &set_according_to_type ("bacteria", "/bacteria_11_collection/");

my $computer = "GIS-KATNISS";
print "collection dir $collection_dir\n";
my $out_root = "/mnt1/SIFT_databases/";
my $out_folder = "metadocs/";
my @orgs = &get_organisms ($ensembl_root, $collection_dir);

foreach my $org_name (@orgs) {
	print "making metadoc for $org_name\n";
	my ($gc_table, $gc_tablename, $mgc_table, $mgc_tablename) = get_geneticCode($org_name);
	print_out_metadoc_for_organism ($computer, $ensembl_root, $collection_dir, $var_dir, $dbsnp_format, $org_name, $out_root, $out_folder, $out_folder . "/summary.txt", $gc_table, $gc_tablename, $mgc_table, $mgc_tablename);
}

sub get_organisms
{
	my ($ensembl_root, $collection_dir) = @_;
	my @orgs;

	my $wget_dir = $ensembl_root . "/gtf/$collection_dir/";
        # remove any existing file else the next download  will be renamed 
        # index.html1
        #print $wget_dir . "\n";
        system ("rm index.html");
        `wget $wget_dir`;
        open (INDEX, "index.html" ) || die "can't open index.html in get_version";
        while (my $line = <INDEX>) {
                if ($line =~ /Directory/) {
			print "reading index\n";
			if ($collection_dir eq "") {
#print "coollection dir is empty\n";
				my $loc = index ($line, "gtf");
       		                 my $phrase = substr ($line, $loc+5, 100);
                        print $phrase . "\n";
       		                 my $rloc = index ($phrase, "/");
				my $org = substr ($phrase, 0, $rloc);
				print $org . "\n";
				push (@orgs, $org);
			} else { # bacteria
				 my $loc = index ($line, $collection_dir);
#                        print "loc $loc\n";
                                my $phrase = substr ($line, $loc+length ($collection_dir) + 1, 100);
                                print $phrase . "\n";
                                my $rloc = index ($phrase, "/");
                                my $org = substr ($phrase, 0, $rloc);
                                print $org . "\n";
                                push (@orgs, $org);
			} # end else
                }
        }
        close (INDEX);
        system ("rm index.html");
	return (@orgs);	
}

sub print_out_metadoc_for_organism 
{
	my ($computer, $ensembl_root, $collection_dir, $var_dir, $dbsnp_format, $common_name, $out_root, $outfolder, $summary_file, $gc_table, $gc_tablename, $mgc_table, $mgc_tablename) = @_;
	my $version = &get_version  ($ensembl_root, $common_name, $collection_dir);
	#print "version $version\n";
	my $date = `date +%F`;
	chomp ($date);
	`echo "$common_name\t$version\t$date" >> $summary_file`;
	my $outfile = $outfolder . "/" . $common_name . "-" . $date . ".txt"; 
	open (OUT, ">$outfile") || die "can't open $outfile\n";
	&say_external_websites ($ensembl_root, $collection_dir, $common_name, $version, $var_dir, $dbsnp_format, $gc_table, $gc_tablename, $mgc_table, $mgc_tablename);
	&say_out_local_outputs ($out_root, $common_name, $version);
	&say_out_MOSIFT ($computer);
	&say_out_constants ();
	close (OUT);
}

sub get_version {
	my ($ensembl_root, $common_name, $collection_dir) = @_;

	print "IN HERE with $ensembl_root $common_name \n";
        my $common_name_ucl = ucfirst ($common_name);
        my $wget_dir = $ensembl_root . "/gtf/$collection_dir/" . $common_name . "/";
	# remove any existing file else the next download  will be renamed 
	# index.html1
	my $ver;
	print $wget_dir . "\n";
	system ("rm index.html");
	`wget $wget_dir`;
	open (INDEX, "index.html" ) || die "can't open index.html in get_version";
	while (my $line = <INDEX>) {
#Pauline added abinitio on Sept 3, because we don't want predicted gene sets, but the offical gene sets
		print "reading lin $line\n";
		if ($line =~ /gtf\.gz/ && !($line =~ /abinitio/) && !($line =~ /chr/)) {
# Sep 3, Pauline added period to rindex, so that 
# Felis_catus.Felis_catus_6.2.81 would work
# March 7, 2016. Pauline added ! match .chr. so that small extra scaffolds
# and chromosomes would be looked up
			my $loc = rindex ($line, $common_name_ucl . ".");
			my $phrase = substr ($line, $loc, 100);
			print "found " . $phrase . "\n";
			$_ = $phrase;
			/$common_name_ucl\.(.*)\.gtf\.gz/;
			print $1;
			$ver= $1;
			last;
		}
	}
	close (INDEX);
	system ("rm index.html");
	print "Version $ver\n";
	return $ver;
}

sub say_external_websites {
        my ($ensembl_root, $collection_dir, $common_name, $version, $var_dir, $dbsnp_format, $gc_table, $gc_tablename, $mgc_table, $mgc_tablename) = @_;

        my $common_name_ucl = ucfirst ($common_name);
        my $gene = $ensembl_root . "/gtf/$collection_dir" . $common_name . "/" . 
			$common_name_ucl . "." . $version . ".gtf.gz";
        my $pep = $ensembl_root . "/fasta/$collection_dir" . $common_name . "/pep/" .
			$common_name_ucl . "." . $version . ".pep.all.fa.gz";
        my $chr = $ensembl_root . "/fasta/$collection_dir" . $common_name . "/dna/";
        my $dbsnp = $ensembl_root . "/" . $var_dir . "/" . $common_name . "/";
        say OUT "";
        say OUT "GENE_DOWNLOAD_SITE=" . $gene ;
        say OUT "PEP_FILE=" . $pep;
        say OUT "CHR_DOWNLOAD_SITE=" . $chr ;
	if ($dbsnp_format ne "NO PRINT") { 
        say OUT "DBSNP_ORGANISM_DOWNLOAD_SITE=" . $dbsnp;
        	if ($dbsnp_format eq "UPPER CASE") { 
			say OUT "DBSNP_VCF_FILE=" . $common_name_ucl . ".vcf.gz";
		} else {

			say OUT "DBSNP_VCF_FILE=" . $common_name . ".vcf.gz";
		
		}	
	}
	say OUT "GENETIC_CODE_TABLE=". $gc_table;
	say OUT "GENETIC_CODE_TABLENAME=".$gc_tablename; 
	say OUT "MITO_GENETIC_CODE_TABLE=".$mgc_table;
	say OUT "MITO_GENETIC_CODE_TABLENAME=".$mgc_tablename;  # end do not print
}

sub say_out_local_outputs
{
        my ($out_root, $org_name, $version) = @_;

        my $PARENT_DIR = $out_root . "/" . $org_name ;
        say OUT "";
        say OUT "PARENT_DIR=" . $PARENT_DIR;
        say OUT "ORG=" . $org_name;
        say OUT "ORG_VERSION=" . $version  ;

        say OUT ""
}

sub say_out_constants
{
	say OUT "# Sub-directories, don't need to change";
        say OUT "GENE_DOWNLOAD_DEST=gene-annotation-src";
        say OUT "CHR_DOWNLOAD_DEST=chr-src";
        say OUT "LOGFILE=Log.txt";
        say OUT "ZLOGFILE=Log2.txt";
        say OUT "FASTA_DIR=fasta";
        say OUT "SUBST_DIR=subst";
        say OUT "ALIGN_DIR=SIFT_alignments";
        say OUT "SIFT_SCORE_DIR=SIFT_predictions";
        say OUT "SINGLE_REC_BY_CHR_DIR=singleRecords/";
        say OUT "SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores";
        say OUT "DBSNP_DIR=dbSNP";


        say OUT "";
        say OUT "# Doesn't need to change";
        say OUT "FASTA_LOG=fasta.log";
        say OUT "INVALID_LOG=invalid.log";
        say OUT "PEPTIDE_LOG=peptide.log";
        say OUT "ENS_PATTERN=ENS";
        say OUT "SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord";
        say OUT "";
}

sub set_according_to_type 
{
	my ($type, $collection_dir) = @_;
	my $ensembl_dir;
	my $var_dir;
	my $dbsnp_file; 

	switch ($type) {
       	 case "common" {
		$ensembl_root = "ftp://ftp.ensembl.org/pub/release-83/";
		$var_dir = "/variation/vcf/";
		$dbsnp_file = "UPPER CASE";
		 last;
	}
       	 case "fungi" {
		$ensembl_root = "ftp://ftp.ensemblgenomes.org/pub/fungi/current";
		$var_dir = "";
		$dbsnp_file = "NO PRINT";
                last;
	}
       	 case "metazoa" {
		$ensembl_root = "ftp://ftp.ensemblgenomes.org/pub/metazoa/current/";
		$var_dir = "/vcf/";
		$dbsnp_file = "LOWER CASE";
	        last;
	}
       	 case "bacteria" { 
		$ensembl_root = "ftp://ftp.ensemblgenomes.org/pub/bacteria/current/";
		$var_dir = "";
#		$collection_dir = "/bacteria_11_collection/";
		$dbsnp_file = "NO PRINT";
                last;
	}
       	 case "plants" {
		$ensembl_root = "ftp://ftp.ensemblgenomes.org/pub/plants/current";
		$var_dir = "/vcf/";
		$dbsnp_file = "LOWER CASE";
               last;
	}
       	 case "protists" { 
		$ensembl_root = "ftp://ftp.ensemblgenomes.org/pub/protists/current/";
		$var_dir = "/vcf/";
		$dbsnp_file = "LOWER CASE";
                last;
	}
       	 default {
                print "Please specify type\n";
                exit (-1);
	}
	} #end switch statement
#	print "collection " . $collection_dir . "\n";
	return ($ensembl_root, $var_dir, $dbsnp_file, $collection_dir);
}

sub say_out_MOSIFT 
{
	my ($computer) = @_;

	my $sift4g_path;
	my $protein_db;

	if ($computer eq "HOME-BUFFY") {
		$sift4g_path = "/mnt1/SIFT4G_nonGPU/sift4g/bin/sift4g";
		$protein_db = "/mnt1/protein_db/uniprot90_Jul2016/uniref90.fasta";
	} elsif ($computer eq "GIS-KATNISS") {
		$sift4g_path = "/mnt1/SIFT4G_2.0.0/bin/sift4g";
		$protein_db = "/mnt1/protein_db/uniprot90_Jul2016/uniref90.fasta"; 
	}

        say OUT "";
        say OUT "#Running SIFT 4G";
#        say OUT "#make sure <PROTEIN_DB>_5_hash.bin and <PROTEIN_DB>_5_index.bin are already created";
	
        say OUT "SIFT4G_PATH=$sift4g_path";
        say OUT "PROTEIN_DB=$protein_db";
	say OUT "COMPUTER=$computer\n";
        say OUT "";
}

#Swarna: Used NCBI eutils for extracting the taxonomy code and efetch to extract the genetic codes
sub get_geneticCode {
        my ($org) = @_;
        my $DB = 'taxonomy';
        my $retmax=10;
        my $utils = 'http://www.ncbi.nlm.nih.gov/entrez/eutils';
        my $esearch = "$utils/esearch.fcgi?db=$DB&retmax=1&usehistory=y&term=";
        my $esearch_res = get($esearch . $org);
        $esearch_res =~ m/<Count>(\d+)<\/Count>.*<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>.*<Id>(\d+)<\/Id>/s;
        my $count = $1;
        if($count>1){
                die("More TaxonIds found for $org in NCBI Taxonomy\n");
        }
        my $QueryKey = $2;
        my $WebEnv = $3;
        my $txId = $4;
        my $efetch = "$utils/efetch.fcgi?retmode=xml&retstart=0&retmax=$retmax&db=$DB&query_key=$QueryKey&WebEnv=$WebEnv";
        my $efetch_res = get($efetch);
        $efetch_res =~ m/<GCId>(\d+)<\/GCId>/s;
        my $GC_table = $1;
        $efetch_res=~m/<GCName>(\S+)<\/GCName>/s;
        my $GC_tablename = $1;
        $efetch_res =~m/<MGCId>(\d+)<\/MGCId>/s;
        my $MGC_table = $1;
        $efetch_res =~m/<MGCName>(.*)<\/MGCName>/s;
        my $MGC_tablename = $1;
        return($GC_table, $GC_tablename, $MGC_table, $MGC_tablename);
}
