#package DNA_PROT;
use strict;
#use Exporter;
#our @EXPORT_OK = qw ( chr_is_mito);
use POSIX;

my %codons = ('TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C', 'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
              'TTA' => 'L', 'TCA' => 'S', 'TTG' => 'L', 'TCG' => 'S', 'TGG' => 'W', 'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H',
              'CGT' => 'R', 'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R', 'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q',
              'CGA' => 'R', 'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R', 'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N',
              'AGT' => 'S', 'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S', 'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K',
              'AGA' => 'R', 'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R', 'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D',
              'GGT' => 'G', 'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G', 'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E',
              'GGA' => 'G', 'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G', 'TAA' => '*', 'TAG' => '*', 'TGA' => '*');

sub chr_is_plastid {
	my ($chr) = @_;
	#Arabidopsis has Pt	
	if ($chr =~ /^P[tT]/ || $chr =~ /chrP[tT]/ || $chr =~ /chloroplast/) {
		return 1;
	} 
	return 0;
}

sub chr_is_mito {
	my ($chr)  = @_;

	# Drosophila has "dmel_mitochondrion"
	# human vertebrate has "MT"
	# one organism had Mito
	# C. elegans has Mt 
	if ($chr =~ /^M[tT]/ || $chr eq "Mito" || $chr =~ /mito/i || $chr =~ /chrM[tT]/) {
		return 1;
	}
	return 0;
}


sub adjust_to_nonstandard_code {
	my ($ncbi_code, %cod_table) = @_;

	# taken from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
	if ($ncbi_code == 1) {
		return %cod_table;
#	} elsif (chr_is_mito ($chr) && $ncbi_code == 2) {
	} elsif ($ncbi_code == 2) { 
		#vertebrate mitochondrial code
		$cod_table{'AGA'} = '*';
		$cod_table{'AGG'} = '*';
		$cod_table{'ATA'} = 'M';
		$cod_table{'TGA'} = 'W';
	} elsif ($ncbi_code == 3) { 
#	} elsif (chr_is_mito ($chr)  && $ncbi_code == 3) {
                #yeast mitochondrial code
                $cod_table{'ATA'} = 'M';
                $cod_table{'CTT'} = 'T';
                $cod_table{'CTC'} = 'T';
                $cod_table{'CTA'} = 'T';
		$cod_table{'CTG'} = 'T';
                $cod_table{'TGA'} = 'W';
		$cod_table{'CGA'} = '';
                $cod_table{'CGC'} = '';
	} elsif ($ncbi_code == 4) {

#	} elsif (chr_is_mito ($chr) && $ncbi_code == 4) {
                # Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
                $cod_table{'TGA'} = 'W';
		# NOTE website says the following which I did not code for#		:
		# (I am assuming same amino acid rather than a start = M
#		Alternative Initiation Codons:
#		Trypanosoma: UUA, UUG, CUG
#		Leishmania: AUU, AUA
#		Tertrahymena: AUU, AUA, AUG
#		Paramecium: AUU, AUA, AUG, AUC, GUG, GUA(?) 
#		Alternate start codons are still translated as Met when 
#		they are at the start of a protein (even if the codon encodes a different amino acid otherwise).
#		To code thi properly, need to know amino acid pos = 1, and ncbi_code
# ERROR need to separate mold, protozoan, coelenterate mito from mycoplase/spiroplasma!!!!
# I checked mold, I don't think Mt is assembled right now
	} elsif ($ncbi_code == 5) {
#	} elsif (chr_is_mito ($chr) && $ncbi_code == 5) {
                # Invertebrate Mitochondrial Code (transl_table=5)
#                $cod_table{'AGG'} = ''
# the codon AGG is absent in Drosophila, but it doesn't say about other invertebrates, so leaving it in, just in case;
                $cod_table{'AGA'} = 'S';
                $cod_table{'AGG'} = 'S';
                $cod_table{'ATA'} = 'M';
                $cod_table{'TGA'} = 'W';
	} elsif ($ncbi_code == 6) {
		# Ciliate, Dasycladacean and Hexamita Nuclear Code 
		$cod_table{'TAA'} = 'Q';
                $cod_table{'TAG'} = 'Q';
	} elsif ($ncbi_code == 9) {
#	}  elsif (chr_is_mito ($chr) && $ncbi_code == 9) {
		# Echinoderm and Flatworm Mitochondrial Code
		$cod_table{'AAA'} = 'N';
                $cod_table{'AGA'} = 'S';
                $cod_table{'AGG'} = 'S';
                $cod_table{'TGA'} = 'W';
	} elsif ($ncbi_code == 10) {
		# Euplotid Nuclear Code
		$cod_table{'TGA'} = 'C';	
	} elsif (  $ncbi_code == 11) {
		# change start position
#		The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
		# different initiation codons	
# mixed -- I think it's the entire Bacterial  & Archael genome, but just plant plastids. 
# what does the NCBI file say?

	} elsif ($ncbi_code == 12) {
		$cod_table{'CTG'} = 'S';
	} elsif ($ncbi_code == 13) {
#	} elsif (chr_is_mito ($chr) && $ncbi_code == 13) {
#		The Ascidian Mitochondrial Code (transl_table=13)
		$cod_table{'AGA'} = 'G';
                $cod_table{'AGG'} = 'G';
                $cod_table{'ATA'} = 'M';
                $cod_table{'TGA'} = 'W';
	} elsif($ncbi_code == 14) {
#	} elsif (chr_is_mito ($chr) && $ncbi_code == 14) {
#		The Alternative Flatworm Mitochondrial Code (transl_table=14)
		$cod_table{'AAA'} = 'N';
                $cod_table{'AGA'} = 'S';
                $cod_table{'TGG'} = 'S';
                $cod_table{'TAA'} = 'Y';
		$cod_table{'TGA'} = 'W';
	} elsif ($ncbi_code == 16) {
#	}  elsif (chr_is_mito ($chr) && $ncbi_code == 16) {
		# Chlorophycean Mitochondrial Code (transl_table=16)
		$cod_table{'TAG'} = 'L';
	} elsif ($ncbi_code == 21) {
#	} elsif (chr_is_mito ($chr) && $ncbi_code == 21) {
		# Trematode Mitochondrial Code (transl_table=21)
		$cod_table{'TGA'} = 'W';
                $cod_table{'ATA'} = 'M';
                $cod_table{'AGA'} = 'S';
                $cod_table{'AGG'} = 'S';
                $cod_table{'AAA'} = 'N';
	} elsif ($ncbi_code == 22) {
#	} elsif (chr_is_mito ($chr) && $ncbi_code == 22) {
#	22. Scenedesmus obliquus Mitochondrial Code (transl_table=22)
		$cod_table{'TCA'} = '*';
                $cod_table{'TAG'} = 'L';
	} elsif ($ncbi_code == 23) {
#	} elsif (chr_is_mito ($chr) && $ncbi_code == 23) {
# 	23. Thraustochytrium Mitochondrial Code (transl_table=23)
		#start code is different
		$cod_table{'TTA'} = '*';
	} elsif ($ncbi_code == 24) {
#	} elsif (chr_is_mito ($chr) && $ncbi_code == 24) {
	# 24. Pterobranchia Mitochondrial Code (transl_table=24)
		$cod_table{'AGA'} = 'S';
                $cod_table{'AGG'} = 'K';
                $cod_table{'TGA'} = 'W';
	} elsif ($ncbi_code == 25) {
#		 Candidate Division SR1 and Gracilibacteria Code (transl_table=25)
		$cod_table{'TGA'} = 'G';
	}  


	return %cod_table;
}	

sub dna_2_codon_2_aa {
    my @dna = @_;
    my $cds = join('', @dna); 
    my @codonArray;
    my $counter = 0;    
    my $cds_length = scalar(@dna)/3;

    if ((scalar(@dna) % 3) != '0'){
	print "INVALID: CDS is not divisible by 3\n";
    }
#    print "CDS: $cds_length..$cds\n";

    while ($cds =~ m/(\w{3})/g){ #splits coding sequence by 3 letter groups
	#print "1: $1\n";
	$counter ++;
	push(@codonArray, $1) if ($counter == $cds_length);
	push(@codonArray, $1) if (($1 !~ m/tga|taa|tag/ig) && ($counter < $cds_length)) ;
	#push(@codonArray, $1);
    }
#    print "COD: $counter..$codonArray[0]..$codonArray[1]..$codonArray[3]\n";
    my $translatedCDS = '';
    my $count;
    foreach my $codon (@codonArray){
	#print "inside: $codon..$codons{$codon}\n";
	$codon = uc($codon);
	$count++;
	if ($codon =~ m/tga|taa|tag/ig) {
	    print "STOP CODON: $count\n";
	}
	else{
	    $translatedCDS = $translatedCDS.$codons{$codon};
	}
    }
    
#    print "PROT: $translatedCDS\n";
    my @aa = split("",$translatedCDS);
    return (\@codonArray,\@aa);    
} #end dna_2_codon_2_aa


sub compare_arrays {
    my ($first, $second) = @_;
    no warnings;  # silence spurious -w undef complaints
    return 0 unless @$first == @$second;
    for (my $i = 0; $i < @$first; $i++) {
	return 0 if $first->[$i] ne $second->[$i];
    }
    return 1;
}  

sub retrieve_seq_from_chr_file {
    my $gz_chr_file = $_[0];
    my $start_pos = $_[1]+1;
    my $end_pos = $_[2];
    #my ($gz_chr_file, $start_pos, $end_pos) = @_;
    my $dna_seq;
    my $totalbases = ($end_pos-$start_pos)+1;

    #print "CHRO FILE: $gz_chr_file $start_pos..$end_pos\n";

    my $seqtmp  = `zcat $gz_chr_file | head -2  | tail -1`;
    $_ = $seqtmp;
    my $linewidth = tr/[a-zA-Z]//;

    #print "CHRO FILE: $start_pos..$end_pos\n";
    if ($start_pos < 0) {$start_pos = 0;}

    if ($linewidth == 0) { return ""; }

    my $linestart = floor ($start_pos / $linewidth) + 1;  # add 1 to take
    # into account  the header > line;
    if ($start_pos % $linewidth == 0) {
	$linestart = $linestart - 1;
    }

    my $numlines =  ceil ($totalbases / $linewidth) + 1;
    $linestart += $numlines;

    my @seq = `zcat $gz_chr_file | head -$linestart | tail -$numlines`;
    my $basestart = $start_pos % $linewidth;
    if ($basestart == 0) { $basestart = $linewidth;}

    my $basecount = 0; my $line = 0;
    while ($basestart <= $linewidth && $basecount < $totalbases) {
	$dna_seq .=  substr ($seq[$line], $basestart - 1, 1);
	$basestart++;
	$basecount++;
    }
    $line++;
    while ( $basecount < $totalbases){
	if ($basecount + $linewidth <= $totalbases) {
	    chomp ($seq[$line]);
	    $dna_seq .= $seq[$line];
	    $basecount += $linewidth;
	    $line++;
        } else {
	    my $index = 0;
	    while ($basecount < $totalbases) {
		$dna_seq .=  substr($seq[$line], $index, 1);
		$index++;
		$basecount++;
	    }
        }
    }
    return $dna_seq;
}



sub trim_last_codon_if_stop {
    my ($seq1_ref, $seq2_ref) = @_;
    my $seq1 = $$seq1_ref;
    my $seq2 = $$seq2_ref;
    
    $_ = $seq1; 
    my $length = tr/A-Za-z//;
    my $seq1lastcodon = substr ($seq1, $length - 3, 3);
    my $seq2lastcodon = substr ($seq2, $length - 3, 3);
    $seq1lastcodon = uc ($seq1lastcodon); 
    $seq2lastcodon = uc ($seq2lastcodon);
#    print "last codon is $seq1lastcodon $seq2lastcodon\n";
    if ($seq1lastcodon eq "TGA" || $seq2lastcodon eq "TGA" ||
	$seq1lastcodon eq "TAG" || $seq2lastcodon eq "TAG" ||
	$seq1lastcodon eq "TAA" || $seq2lastcodon eq "TAA") {
	$$seq1_ref = substr ($seq1, 0, $length - 3);
	$$seq2_ref = substr ($seq2, 0, $length - 3);
    }
}


sub mutants {
    my($codon) = $_[0];
    my @mutants;
    my @mut_aa;
    my @mut_base;

    #my(@nucleotides) = ('A', 'C', 'G', 'T');
    my(@nucleotides) = ('a', 'c', 'g', 't');
    

    for (my $i =0; $i < 3; $i++){
	my $mut_codon = $codon;
	foreach my $nuc (@nucleotides){
	    substr($mut_codon,$i,1,$nuc);
	    if ($mut_codon !~ /$codon/i){
		my $aa = $codons{uc($mut_codon)};
		my $base = uc($nuc);
		#print "Trying this:$codon..$mut_codon..$base..$aa\n";
		push (@mutants, $mut_codon);
		push (@mut_aa, $aa);
		push (@mut_base, $base);
		
	    }	

	}
    }
    # Insert the random nucleotide into the random position in the DNA
    # The substr arguments mean the following:
    # In the string $dna at position $position change 1 character to
    # the string in $newbase
    return (\@mut_base,\@mutants,\@mut_aa);
}


sub split_fasta_60 {

    use vars qw ($ofile $line $f_cnt );
    my $fil = $_[0];
    
    $ofile = $_[1];

#    print "split_fasta_60: Reading in from $fil\n";

    open(INFILE,"<$fil") || die(" could not open $fil");

    #$ofile = ""; 
    $f_cnt = 0;
    my @lin;
    while ($line = <INFILE> ) {
	chomp($line);
	if (substr($line,0,1) eq '>') {
	    if ($ofile) {close (OUTFILE);}
#	    $ofile = "temp.fasta";
#	    print "Writing out to $ofile\n";
	    $f_cnt += 1;     # $f_cnt ++;
	    open(OUTFILE,">$ofile");

	    print OUTFILE "$line";
	}
	else{
	    my @seq = split ("",$line);
	    my $i;
	    my $j=length($line);
	    for ($i=0;$i<$j;$i++){
		if ($i % 60 == 0){
		    print OUTFILE "\n$seq[$i]";
		}

		else{
		    print OUTFILE $seq[$i];
		}
	    }
	}
    }
    print OUTFILE "\n";
    close OUTFILE;
    close INFILE;
}

sub reverse_dna {

    my $DNA = $_[0];
    my $revcom = reverse $DNA;

    $revcom =~ tr/ACGTacgt/TGCAtgca/;
    #print "Here is the reverse complement DNA:\n\n$revcom\n";
    
    return $revcom;

}

#1;

