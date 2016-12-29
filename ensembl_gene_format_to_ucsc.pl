#!/usr/bin/perl
use strict;
use Class::Struct;
use Switch;
require 'readMeta.pl';

struct Transcript => {
        transcript_id =>  '$',
        biotype => '$',
	gene_name => '$',
	gene_id => '$',
        transcript_name => '$',
	protein_id => '$',
        type => '$',
        chrom => '$',
        orientation => '$',
        tx_start => '$',
        tx_end => '$',
        cds_starts => '$',
        cds_ends => '$',
	exon_starts => '$',
	exon_ends => '$',
	utr_starts => '$',
	utr_ends => '$',
        #ensembl fields
        start_beg => '$',
        start_end => '$',
        stop_beg =>'$',
        stop_end => '$',
	frame => '$', #only used if start_codon is not available 
	# calculated by me because not every gene has start_codon and stop_codon 
	cds_start => '$',
	cds_end => '$',  
	transcript_type => '$',
};
my $TRUE = 1;
my $FALSE = 0;

my $file = $ARGV[0]; 
my $outfile = $ARGV[1]; 

sub adjustForCDSStartAndEnd {
#    print "Adjust for CDS Starts and End...\n" .
#       "Remove exons/exon parts NOT within " .
#       "coding regions starts and ends.\n";
    my ($tx) = @_;

#    my ($e_starts_aref, $e_ends_aref, $cdsS, $cdsE) = @_;
    my @exon_starts = @{$tx->exon_starts} ; 
    my @exon_ends = @{$tx->exon_ends}; 
    my $cdsS ; 
    my $cdsE ; 
    my $orientation;
    assign_cds_start_end ($tx);
    if ($tx->orientation eq "-") {
        $orientation = -1;
    }
    $cdsS = $tx->cds_start;
    $cdsE = $tx->cds_end;


    my @adjusted_exon_starts = (); # to store exons after adjustments
    my @adjusted_exon_ends = ();

    for (my $i = 0; $i < scalar(@exon_starts); $i++) {
        my $exon_start = $exon_starts[$i];
        my $exon_end = $exon_ends[$i];

       # Cases:
        # (1) exon end < coding start: discard the exon
        # (2) exon start < coding start but coding start < exon end < coding end: trim it
        # (3) coding start < exon start < coding end but coding end < exon endL trim it
        # (4) coding end < exon start: discard the exon

        my $include = $TRUE;
        if ($exon_end < $cdsS || $cdsE < $exon_start) {
            $include = $FALSE; # discard because exon is outside of coding region (alternate splicing?)
		# utr
		push (@{$tx->utr_starts}, $exon_start);
                push (@{$tx->utr_ends}, $exon_end);

        } else {
            # Bec. cases totally outside the coding region have been taken into account of
            # we now know we need to include the exon, the only thing left to do is to
            # decide whether to trim or not
            $include = $TRUE;
            if ($exon_start < $cdsS) { 
		push (@{$tx->utr_starts}, $exon_start);
                push (@{$tx->utr_ends}, $cdsS);
		$exon_start = $cdsS;
	     }
            if ($cdsE < $exon_end) { 
		push (@{$tx->utr_starts}, $cdsE);
                push (@{$tx->utr_ends}, $exon_end);
		$exon_end = $cdsE;
	    }
        }
 if ($include == $TRUE) {
            push @adjusted_exon_starts, $exon_start;
            push @adjusted_exon_ends, $exon_end;
        }
    } #end for
	

    $tx->cds_starts (\@adjusted_exon_starts); 
    $tx->cds_ends (\@adjusted_exon_ends); 

} #end adjustForCDSStartAndEnd


sub assign_cds_start_end 
{
	my ($tx) = @_;

	# use start_codon and stop_codon to get CDS starts and ends
	# but if it's not there, use exon starts and ends
       my @cds_start_array = sort {$a <=> $b}  @{$tx->cds_starts};
       my $chr_start = $cds_start_array[0];
       my @cds_end_array = sort {$a <=> $b} @{$tx->cds_ends};
       my $chr_end = $cds_end_array[$#cds_end_array];

	my $bases_to_trim_off = 0;
	if ($tx->frame == 2) {
		$bases_to_trim_off = 2;
	} elsif ($tx->frame == 1) {
		$bases_to_trim_off = 1;
	}

# Pauline
# March 22, 2014 is CDS start incomplete, just take 3 bp upstream so that 
# coordinate is correct. otherwise a check with VEP shows we say L2 while 
# VEP has L3 (because we threw out the first amino acid
# first amino acid will be incorrect, but positioning should be OK
# April 5, 2014, this doesn't work and was removed

	$tx->transcript_type ("");
	if ($tx->orientation eq "-") {
                if ($tx->start_end > 0) {
			$tx->cds_end ( $tx->start_end);
		} else {
			$tx->cds_end ($chr_end - $bases_to_trim_off); 
		}
		if ($tx->stop_beg > 0) {
	                $tx->cds_start  ($tx->stop_beg);
		} else {
			$tx->cds_start ($chr_start );
		}
        } elsif ($tx->orientation eq "+") {
                if ($tx->start_beg > 0) {
			$tx->cds_start ( $tx->start_beg);
                } else {
			$tx->cds_start ($chr_start + $bases_to_trim_off); 
		}
		if ($tx->stop_end > 0) {
			$tx->cds_end ( $tx->stop_end);
        	} else {
			$tx->cds_end ($chr_end);
		}
	}

		
}

sub print_out_utr 
{
	my ($tx) = @_;
	my @lines_to_print; 
	my @utr_starts = @{$tx->utr_starts};
	my @utr_ends = @{$tx->utr_ends};
	if (scalar @utr_starts == 0) {
#		print "no utr's at all " . $tx->transcript_id . "\n"; 
		return @lines_to_print;
	}
	for (my $i = 0; $i < @utr_starts; $i++) {
		my $utr_type = "UTR";
		if ($utr_ends[$i] <= $tx->cds_start) {
			if ($tx->orientation eq "+") {
				$utr_type = "UTR_5";
			} elsif ($tx->orientation eq "-") {
				$utr_type = "UTR_3";
			}
		} elsif ($utr_starts[$i] >= $tx->cds_end) {
			if ($tx->orientation eq "+") {
                                $utr_type = "UTR_3";
                        } elsif ($tx->orientation eq "-") {
                                $utr_type = "UTR_5";
                        }
		} else {
			print "dont understand utrs " . $tx->transcript_id . "\n";
			print $utr_ends[$i] . " " .  $tx->cds_start . " " .
			 $utr_starts[$i] . " " .  $tx->cds_end . "\n";	
#			exit (-1);
		}
		my $line =  $tx->chrom . "\t" . $utr_type . "\t" .
			($utr_starts[$i] + 1) . "\t" . $utr_ends[$i] .
			"\t" .  
			 $tx->gene_id . "\t" .
                        $tx->transcript_id . "\t" .
                        $tx->gene_name . "\t" .
                        $tx->transcript_name . "\t" . $tx->orientation; 
		push (@lines_to_print, $line);
	}
	return (@lines_to_print);
}

sub print_out_transcript
{
	my ($tx) = @_;
	my $orientation = 1;
	my $cds_start = $tx->cds_start;
	my $cds_end = $tx->cds_end;
	if ($tx->orientation eq "-") {
		$orientation = -1;
	} 	
        my $line =  $tx->chrom . ":" . $orientation . ":" . $cds_start 
			 . ":". $cds_end ;

	my $num_exons = scalar (@{$tx->cds_starts});

        my $start_string = join ("," , sort {$a <=> $b} @{$tx->cds_starts}); 
        my $end_string = join (",", sort {$a <=> $b} @{$tx->cds_ends}); 
        #ensembl fields
	$line .= ":" . $num_exons . ":" . $start_string . ":" . $end_string;
	$line .= ":" . $tx->transcript_id . ":" . $tx->gene_id . ":" .
			$tx->gene_name . ":" .  $tx->protein_id. ":" .
			$tx->transcript_name . ":" . $tx->biotype; 
	print $line . "/n";
	return $line ;
}

sub trim {
        my $s = shift;
        $s =~ s/^\s+|\s+$//g;
        return $s;
}

sub get_region_info {
        my ($gene_info_string) = @_;
	#print $gene_info_string;
        my @gene_desc = split (";", $gene_info_string);
        my %info;
        for (my $i = 0; $i <@gene_desc; $i++) {
                my $_ = $gene_desc[$i];
                /(.*)"(.*)"/;
                my $key = trim ($1);
                my $val = trim ($2);
		# June 27, 2014. Added below because some gene files had
		# colons in gene id's, and that was screwing up later
		# parsing as colons are used as a delimiter elsewhere
		$val =~ s/\:/\./g;
                $info{$key} = $val;
#		print "key $key val $val \n";    
    }
	return %info;
}

sub initialize_transcript {

        my ($tx, %info_hash) = @_;
        $tx->gene_id ( );
        $tx->protein_id ("") ;
	$tx->transcript_id ( $info_hash{"transcript_id"} ); 
        $tx->gene_id($info_hash{"gene_id"}); 
        $tx->gene_name ($info_hash{"gene_name"}); 
	$tx->transcript_name($info_hash{"transcript_name"});
	my @array1; my @array2; my @array3; my @array4; my @array5; my @array6;
        $tx->cds_starts (\@array1); 
        $tx->cds_ends (\@array2); 
	$tx->exon_starts (\@array3); 
	$tx->exon_ends (\@array4);
	$tx->utr_starts (\@array5);
	$tx->utr_ends (\@array6);
        $tx->tx_start (0); 
        $tx->tx_end (0); 
	$tx->start_beg (0); 
        $tx->start_end (0); 
        $tx->stop_beg (0);
        $tx->stop_end (0);
	$tx->protein_id ("");
}

sub check_chr_orientation
{
        my ($tx, $chr, $orientation) = @_;
        if ($tx->chrom == "") {
                $tx->chrom ($chr);
        } elsif ($tx->chrom != $chr) {
                print "chrs have changed \n";
		print_out_transcript ($tx);
                exit (-1);
        }
        if ($tx->orientation == "") {
                $tx->orientation ($orientation);
        } elsif ($orientation != $tx->orientation) {
                print "orientations have changed \n";
		print_out_transcript ($tx);
                exit (-1);
        }
}

sub add_info_to_transcript
{
        my ($tx, $line) = @_;
        my @fields = split (/\t/, $line);
	# adjust to 0-count coordinate system
        my $beg = $fields[3] - 1;
        my $end = $fields[4];
        check_chr_orientation ($tx, $fields[0], $fields[6]);
        switch ($fields[2]) {
                case "exon" {
                        push (@{$tx->exon_starts}, $beg);
                        push (@{$tx->exon_ends}, $end);
                        #print "entered exon $line \n";         
                 }
                case "CDS" {
                        push (@{$tx->cds_starts}, $beg);
                        push (@{$tx->cds_ends}, $end);
			if ($line =~ /exon_number \"1\"/) {
				$tx->frame ($fields[7]);
                        #print "entered CDS $line \n";  
			}
                }
                case "stop_codon" {
                        $tx->stop_beg ( $beg);
                        $tx->stop_end ( $end);
                        #print "entered stop codon $line \n";  
                }
                case "start_codon" {
                        $tx->start_beg ( $beg);
                        $tx->start_end ($end);
                        #print "stop codon $line \n";  
                }
                else {
                        #print "did not process $line";
		}
        } # end switch statement
}

sub get_biotype_from_geneinfo{
	my ($geneinfo) = @_;
	$geneinfo =~ /gene_biotype\s\"(\S+)\"\;/;
	my $type = $1;
	#print $type . "\n";
	return($type);
}

## GET FILE  
my ($metafile ) = @ARGV;
my $meta_href = readMeta($metafile);
my %meta_hash = %{$meta_href};

my $infile_dir =  $meta_hash{"PARENT_DIR"} . "/" . $meta_hash{"GENE_DOWNLOAD_DEST"}; 
print $infile_dir . "\n";
my @files = `ls $infile_dir`;
print $infile_dir . "\n";
my $infile = "";
my $outfile = $infile_dir  . "/protein_coding_genes.txt";
my $outfile_noncod = $infile_dir . "/noncoding.txt";
for (my $i=0; $i < @files; $i++) {
	chomp ($files[$i]);
	if ($files[$i] =~ /\.gtf$/ || $files[$i] =~ /\.gtf\.gz$/) {
		$infile = $infile_dir . "/" . $files[$i];
	}
}

if ($infile =~ /\.gz$/) {
        open (IN, "gunzip -c $infile |") || die "can't open pipe to $infile";
    } else {
        open(IN, "<$infile") || die "Unable to open $infile for reading\n";
   }
print "infile " . $infile . "\n";
#### END GET FILE

#my $infile = "Homo_sapiens.GRCh37.74.gtf"; 
#open (IN, $infile) || die "can't open $infile";

open (OUT, ">$outfile") || die "can't open $outfile";
open (NONCOD, ">$outfile_noncod") || die "can't open $outfile_noncod";

my %transcript_hash;
my $line;
GENE_LINE: while ($line = <IN>) {
	if($line =~ /^#/) {
		next GENE_LINE;
	}
        my ($chr,  $dummy, $region_type, $start, $end, $d1, $orientation, $frame, $geneinfo)
                                 = split (/\t/, $line);
        my $biotype = get_biotype_from_geneinfo ($geneinfo);             
        my %geneinfo_hash = get_region_info ($geneinfo);
        my $tx_id = $geneinfo_hash{"transcript_id"};
#	print $tx_id . " for $line\n";
	if ($tx_id ne "") {
	if ($biotype eq "protein_coding") {
		if (!exists ($transcript_hash{$tx_id})) {
       	        	 my $new_gene = Transcript->new();
       	        	 initialize_transcript ($new_gene, %geneinfo_hash);
	               	 $new_gene->chrom ($chr);
			 $new_gene->biotype ($biotype);
			 $new_gene->orientation ($orientation);	
			 $transcript_hash{$tx_id} = $new_gene;
       		 }
	        my $tx = $transcript_hash{$tx_id};
		# add ENSP
		if ($region_type eq "CDS" && $tx->protein_id eq "") {
			%geneinfo_hash = get_region_info ($geneinfo);
			$tx->protein_id ($geneinfo_hash{"protein_id"});
		}
                add_info_to_transcript ($tx, $line);
        } elsif ($biotype =~ /RNA/) {
		print NONCOD "$chr\t$biotype\t" .  $start .
			"\t" .  $end . "\t" . 
			$geneinfo_hash{"gene_id"} . "\t" . 
			$geneinfo_hash{"transcript_id"} . "\t" .
			$geneinfo_hash{"gene_name"} . "\t" .	
			$geneinfo_hash{"transcript_name"} . "\t" .
			$orientation . "\n";
	}
	} # end check that tx_id is not empty
#	print_out_transcript ($tx);
}
close (IN);

foreach my $tx_id (keys %transcript_hash) {
        my $tx = $transcript_hash{$tx_id};
        adjustForCDSStartAndEnd ($tx);
	print OUT print_out_transcript ($tx) . "\n";

	my @utr_lines = print_out_utr ($tx) ;
	if (scalar(@utr_lines) >= 1) {
		for (my $i = 0; $i < @utr_lines; $i++) {
			print NONCOD $utr_lines[$i] . "\n";
		}
	}
}
close (OUT);
close (NONCOD); 
exit (0); 
