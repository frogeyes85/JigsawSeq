#!/usr/bin/perl

# Detect and filter out the reads from others (hg, pcr amplicons)
# Requirement: bwa 0.7.9a-r786
# last modified: Jul-25-2014
# Developed by Jung-Ki Yoon

# bwa index -a is ref.fa
# bwa mem ref.fa TMP.fa

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./filter_reads_from_other.pl [input: forward] [input: reverse] [input: reference] [output: forward] [output:reverse]\n";
die $usage unless ($#ARGV == 4);
my ($in_for, $in_rev, $ref_fname, $out_for, $out_rev) = @ARGV;
my $t_begin = new Benchmark;
my $sam_fname = "TMP_" . $out_for . ".sam";
print "[Program:filter_reads_from_other] input: $in_for, $in_rev --(align to)--> $ref_fname output: $out_for, $out_rev\n";

JigsawSeq::call_sys("./bwa mem $ref_fname $in_for $in_rev > $sam_fname");
my $t_end = new Benchmark; 
print "[Process:detect_seeds] Alignment was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";

# MG00HS19:389:HA8RHADXX:1:2216:20195:101006      77      *       0       0       *       *       0       0       CGGATGAACAGGCAGACATCTGTGAATCGCTTCACGACCACGCTGATGAGCTTTACCGCAGCTGCCTCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGG      @@@F;?BDDHHHGFGIIEEHGIIIEGIGFHGGIIIHGIGGIIIIIBHHIEHIIIEIIGHADFCEECCCD>?@D>B@CDD??@CDDCDD@?2:>:ABCDDDD:CCCDCDDDDDDDD@@BBDDDBBBDCCDDCACDDDDCCCDDDDD>DDDDD AS:i:0  XS:i:0

open(IN, "<$sam_fname");
open(OUT1, ">$out_for");
open(OUT2, ">$out_rev");
my $num_lines=0;
while(<IN>){ 
	next if (substr($_, 0, 1) eq "@");
	$num_lines++;
	if ($num_lines % 1000000 == 0) {print "[Process:detect_seeds] $num_lines lines were processed.\n";}
	my ($qname, $flag, $rname, $pos, $mapQ, $CIGAR, $rnext, $pnext, $tlen, $str, $qual, $NM, $MD, $AS, $XS,) = split /\s+/, $_;
	if ($pos > 0){
	    next;  # aligned to hg
	}else{
	    if ($flag < 100){ # for
		print OUT1 "@". "$qname 1:N:0:AGTTCC\n$str\n+\n$qual\n";
	    }else{ #rev
		print OUT2 "@". "$qname 2:N:0:AGTTCC\n$str\n+\n$qual\n";
	    }
	}
}
close(IN);
close(OUT1);
close(OUT2);

my $t_end = new Benchmark; 
print "[Process:detect_seeds] Detection of seeds were completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;

