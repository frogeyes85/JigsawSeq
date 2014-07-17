#!/usr/bin/perl

# Converted K-mer into fasta format file to detect initial/terminal seeds later. 
# last modified: Jul-04-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./Kmer2fa.pl [input: graph] [output: fasta]\n";
die $usage unless ($#ARGV == 1);
my ($in_fname, $out_fname) = @ARGV;

print "[Program:Kmer2fa] input: $in_fname -> output: $out_fname\n";

my $t_begin = new Benchmark;
open(IN, "<$in_fname") or die "[Error] Can't open $in_fname.\n";
open(OUT, ">$out_fname");
my $num_lines=1;
while(<IN>){ 
	if (substr($_, 0, 1) eq "#"){
		next;
	}
#	if ($num_lines % 1000000 == 0) print $num_lines, "K-mers were processed.\n";
	my ($str, $DP,) = split /\s+/, $_;
	print OUT ">$num_lines\_$DP\n";
	print OUT "$str\n";
	$num_lines++;
}
close(IN);
close(OUT);

my $t_end = new Benchmark; 
print "[Report:Kmer2fa] $num_lines K-mers were recorded to a fasta format file.\n";
print "[Process:Kmer2fa] Converting was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;
