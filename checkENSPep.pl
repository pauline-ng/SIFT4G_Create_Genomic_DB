#!/usr/bin/perl -w
use Bio::SeqIO;

########################################################
#
# We want to ensure the Peptide sequence is correct
# Take the <model_org>_peptide.log file
# Compare the peptides against ensPep.
#
# Author: Sim Ngak Leng
# First created: 2011-06-27
# Last modified: 2011-06-28
#
########################################################


use strict;
#require './common-utils.pl';
use File::Basename;
use Cwd qw(abs_path);
my $directory_of_script = dirname(abs_path(__FILE__));
require $directory_of_script . '/common-utils.pl';



if (scalar(@ARGV) != 3) {
    die "Usage: perl $0 <metafile> <ens pattern> <output file>\n" .
	"Example: perl $0 human.txt ENST checkAgainstENS.txt\n";
}

my ($metafile, $pattern, $outFile) = @ARGV;
my $meta_href = readMeta($metafile);
my %meta_hash = %{$meta_href};

my $ensPepZipDir =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"GENE_DOWNLOAD_DEST"} . "/" ;
my $peptideFile = $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"PEPTIDE_LOG"};

my  $ensPepZipFile = &getPeptideFile ($ensPepZipDir);
my ($numEnsSeq, $ensPep_href) = &getENSPepZipFileData($ensPepZipFile); # The peptide sequences provided by Ensembl
my ($numOurEnsSeq, $ourPep_href) = &getOurPeptideData($peptideFile); # The fasta sequences we generated

my ($matched, $total_num, $unmatched_href) = &performChecks($ensPep_href, $ourPep_href);
&writeOutToFile($matched, $total_num, $unmatched_href, $outFile);



### SUB ROUTINE ########################################
sub getPeptideFile () {
	my ($dir) = @_;
	my @files = glob ("$dir/*pep.all.fa.gz");
	if (@files > 0) {
		chomp ($files[0]);
		return $files[0];
	} else {
		print "no peptide file to compare it against";
		exit (-1);
	}
}

sub writeOutToFile() {
    my ($matched, $total_num, $href, $output) = @_;
    my $unmatched = $total_num - $matched;
    my %unmatched_hash = %{$href};
    print "Writing out to file: $output\n";

    my @keys = keys(%unmatched_hash);

    open(OUT, ">$output") || die "Unable to open $output for writing.\n";
    print OUT "$matched out of $total_num peptide sequences matched\n";
    print OUT "$unmatched out of $total_num peptide sequences from our generated sequences did not match.\n\n";
    foreach my $key (@keys) {
	my $seqs = $unmatched_hash{$key};
	print OUT "$key\n$seqs\n\n";
    }
    close(OUT);
    print "Completed writing to $output\n";
} #end writeOutToFile


sub performChecks() {
    my ($ens_href, $our_href) = @_;
    my %ens = %{$ens_href};
    my %our = %{$our_href};
    my @ensKeys = (keys %our);

    my %unmatched = ();
    my $numberOfMatches = 0;

    foreach my $key (@ensKeys) {
	if (defined $ens{$key}) {
	    my $ensSeq = $ens{$key};
	    my $ourSeq = $our{$key};
#	    print "$ensSeq\n$ourSeq\n";
	    if (uc($ensSeq) eq uc($ourSeq)) {
		$numberOfMatches++;
	    } else {
		$unmatched{$key} = $ensSeq . "\n" . $ourSeq;
	    }
	}
    } #end foreach
    my $totalNumber = scalar(@ensKeys);
    print "$numberOfMatches out of $totalNumber peptide sequences matched\n";
    return ($numberOfMatches, $totalNumber, \%unmatched);

} #end performChecks




sub getENSPepZipFileData() {
    my $ensZip = $_[0];
    print "Getting ens peptide data from $ensZip\n";
    my $tmpdir = $meta_hash{"PARENT_DIR"} . "/tmp";
    if (! -d $tmpdir) { mkdir($tmpdir, 0755); }
    if (! -e $ensZip) { die "$ensZip does not exist, unable to continue.\n"; }

    my %results = ();

    my $tmpfile = $tmpdir . "/checkAgainstEnsemblPeptide.fa";
    my $cmd = "zcat $ensZip > $tmpfile";

    my $res = system($cmd);
    if ($res == 0) {
	my $seqio = Bio::SeqIO->new (-file => $tmpfile , '-format' => 'Fasta');
	while (my $seq = $seqio->next_seq) {
		my $seq_id = $seq->id;
		my $sequence = $seq->seq;
		$results{$seq_id}  = $sequence
	}
	#open(TMP, "<$tmpfile") || die "Unable to open $tmpfile for reading.\n";
	#while(my $line = <TMP>) {
	#    chomp $line;
	#    my ($ens_key, $seq) = split(/\s+/, $line);
	#    $results{$ens_key} = $seq;
	#}
	#close(TMP);
    }
    unlink($tmpfile);

    
    my $numEnsSeq = scalar(keys(%results));

    print "Completed reading in $numEnsSeq Ensembl sequences.\n";
    return ($numEnsSeq, \%results);

} #end getENSPepZipFileData

sub getOurPeptideData() {
    my $peptideFile = $_[0];
    print "Getting our peptide data from $peptideFile\n";
    my %results = ();
    open(FILE, "<$peptideFile") || die "Unable to open $peptideFile for reading.\n";
    my $ens_key = "";
    my $seq = "";
    while(my $line = <FILE>) {
	chomp $line;
	if ($line ne "") {
	    if ($line =~ /^$pattern/) {
		$ens_key = $line;
	    } else {
		if ($line =~ /\*$/ || $line =~ /\X$/) {
		    chop $line;
		    $seq = $line;
		}
	    }
	    if ($ens_key ne "" && $ens_key =~ /$pattern/ && $seq ne "") {
		$results{$ens_key} = $seq;
		$ens_key = "";
		$seq = "";
	    }
	}
    }
    close(FILE);
    
    my $numEntries = scalar(keys(%results));
    print "Completed reading in our generated peptides that we know exists in Ensembl.\n";
    print "There are $numEntries sequences generated.\n";

    return ($numEntries, \%results);

} #end getOurPeptideData




__END__
