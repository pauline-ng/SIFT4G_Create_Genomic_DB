#!/usr/bin/perl

# need to put this in /etc/cron.hourly

my $dbs_to_make_list = "/mnt1/scripts_ensembl/orgs_to_run";
my $done_file = $dbs_to_make_list . ".completed";
my $logfile = $dbs_to_make_list . ".log";
my $scripts_dir = "/mnt1/scripts_ensembl/"; 

open (LOG, ">>$logfile") ;

my $date = `date`;
print LOG $date;
my $running = `ps -eo \"%a\" | grep make-SIFT-db-all | grep ^perl`;
chomp ($running);
print LOG "$running\n";
if ($running ne "") {
	print LOG "something is running $running\n";
} else {
	# start a new job
	# no locks on the files!
	my $org_file = get_first_line ($dbs_to_make_list, $done_file);
	my $commandline = "cd $scripts_dir; perl /mnt1/scripts_ensembl/make-SIFT-db-all.pl metadocs/$org_file";
	print LOG "running $commandline\n";
	system ("$commandline");
}
close (LOG);
 
exit (0);

sub
get_first_line
{
	my ($inlist, $outlist) = @_;
	open (IN, $inlist) || die "can't open $inlist";
	my $line = <IN>;
	chomp ($line);
	print $line;
	my @fields = split (/\s+/, $line);
	my $file_to_run = $fields[0];
	my $tmpfile = $inlist . ".tmp";
	# remove first line from $inlist because it has been executed
	system ("tail -n +2 $inlist > $tmpfile"); 
	system ("mv $tmpfile $inlist");
	# add line to the outlist
	system ("echo \"$line\" >> $outlist");
	close (IN);
	return ($file_to_run);
}
