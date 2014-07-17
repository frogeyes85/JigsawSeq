#!/usr/bin/perl

# Convert candidate contigs into fasta format to mapping raw data
# last modified: Jul-04-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./contigs2fa.pl [input: contigs] [output: fasta]\n";
my ($in_fname, $out_fname) = @ARGV;
die $usage unless ($#ARGV == 1);

my $t_begin = new Benchmark;
open(IN, "<$in_fname") or die "[Error] Can't open $in_fname.\n";
open(OUT, ">$out_fname");
my $num_lines=1;
while(<IN>){ 
	next if (substr($_, 0, 1) eq "#");
	my ($temp, $status, $seq, $init_seed, $term_seed,) = split /\s+/, $_;
	next unless ($status eq "T");		# consider the contigs with terminal seeds only
	print OUT ">$in_fname\_$num_lines\n";
	print OUT "$init_seed$seq$term_seed\n";
	$num_lines++;
}
close(IN);
close(OUT);
my $t_end = new Benchmark; 
print "[Report:contigs2fa] ", $num_lines-1, " contigs were recorded in a fasta format file($out_fname).\n";
print "[Process:contigs2fa] Converting was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
exit;
