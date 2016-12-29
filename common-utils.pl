#!/usr/bin/perl -w

########################################################
#
# Common functions
# Author: Sim Ngak Leng
# First created: 2011-05-16
# Last modified: 2011-05-16
#
########################################################

use strict;

sub getChr () {
	my ($gene_dir) = @_;
	my $gene_file = `ls $gene_dir/*.gtf.gz`;
	chomp ($gene_file);

	my @chr_array = `zcat $gene_file | cut -f1 | uniq | grep -v "^#"| sort | uniq`;
	foreach my $chr (@chr_array) {
		chomp ($chr); 
	}
#	my ($chr_dir) = @_;
#	my @chr_files = `ls $chr_dir/*.fa.gz`;
#	my @chr_array;
#	if (@chr_files == 0) {
#		@chr_files = `ls $chr_dir/*.fa`;
#	}	

#	for (my $i= 0 ; $i < @chr_files; $i++) {
#		$_ = $chr_files[$i];
#		if ($chr_files[$i] =~ /chromosome/) {
#			/chromosome\.(.*)\.fa/;
#			push (@chr_array, $1);
#		} elsif ($chr_files[$i] =~ /genome/) {   
#                        /(.*)\.dna\.genome/;
#                        push (@chr_array, $1);
#		} 
#	}
#	print join ('\n', @chr_array);
	return @chr_array;
}

sub readMeta() {
    my ($metafile) = @_;
    my %results = ();
    open(META, "<$metafile") || die "Unable to open $metafile\n";
    while(my $line = <META>) {
	chomp $line;
	if ($line ne "" && $line !~ /^#/) {
	    my ($key, $value) = split("=", $line);
	    $results{$key} = $value;
	}
    } #end while
    close(META);
    return \%results;
} #end readMeta


sub trim() {
    my ($str) = @_;
    if ($str ne "") {
	$str =~ s/^\s+//;
	$str =~ s/\s+$//;
    }
    return $str;
}

sub deldir {

    my ($dirtodel) = @_;
    print "Deleting $dirtodel\n";
    my $sep = '/';
    opendir(DIR, $dirtodel);
    my @files = readdir(DIR);
    closedir(DIR);
    
    @files = grep { !/^\.{1,2}/ } @files;
    @files = map { $_ = "$dirtodel$sep$_"} @files;
    @files = map { (-d $_)?deldir($_):unlink($_) } @files;
    
    rmdir($dirtodel);
}

sub complement_base() {
    my ($orig_base) = @_;
    if (uc($orig_base) eq "A") { return "T"; }
    elsif (uc($orig_base) eq "T") { return "A"; }
    elsif (uc($orig_base) eq "G") { return "C"; }
    elsif (uc($orig_base) eq "C") { return "G"; }
}

sub make_dir {
	my ($dir) = @_; 
	if (! -d $dir) { 
		print "entered mkdir $dir\n";
		mkdir($dir, 0775) or die $!; 
        }

}

1;
