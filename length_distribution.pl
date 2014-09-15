#!/usr/bin/perl
use strict;

=pro
@HWI-D00236:112:H94UCADXX:2:1101:1569:2182
AGGTTACTCATATATACTTTAGATTGATTTAAAACTTCATTTTTAATTTAAAAGGATCTAGGTGAAGATCCTTTTTGATAATCTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAA
+
@CCDFFFFHGHHHIJJJJJJIIHIJJHJJJJIIJJJJJJJJJIJJJJJJJJJJJJIJJJJJJIIIIIJJJIJIJJJJJJJJJHHHHHHHFFFFDDECCDDDCDDDCBDDADDDDDBDDDEDDDDDDDBDDDDDDDDDDDDDDDDDDDDDD
@HWI-D00236:112:H94UCADXX:2:1101:1967:2120
TGAATGAACTGCAAGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTGCTCGACGTTGTCACTGAAGCGGGAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC
+
    BDDDDFHFHHEEAHIB9FCDBEFDGBAFHIIDC=;E?A:>?>5>A?@B@BBB8555<AC>C45555>>C:>:4>?9B>08A<:>CCCC>@CBBBB5<8>>8<2<@<<<@8<89BAC834>@>C@ACC@CCCC::?0288@<889?:
=cut
my $usage = "./length_distribution.pl [for_fastq] [rev_fastq] [output]\n";
if ($#ARGV != 2){die $usage;}
my ($in_for, $in_rev, $out_fname) = @ARGV;
my %L;
my $num=0;
my @bin = qw(75 90 105 120 135 150);
my @S;
open(IN, "<$in_for") or die "Can't open $in_for.\n";
while(<IN>){
    if (substr($_, 0, 1) ne "\@"){
	die "Error: $_"; 
    }
    my ($seq, ) =split /\s+/, <IN>;
    $L{length($seq)}++;
    <IN>;
    <IN>;
    $num++;
}
close(IN);
open(IN, "<$in_rev") or die "Can't open $in_rev.\n";
while(<IN>){
    if (substr($_, 0, 1) ne "\@"){
	die "Error: $_"; 
    }
    my ($seq, )=split /\s+/, <IN>;
    $L{length($seq)}++;
    <IN>;
    <IN>;
    $num++;
}
close(IN);

open(OUT, ">$out_fname");
print OUT "Total_reads= $num\n";
my $total_sum=0;
foreach my $k (keys %L){
    if ($L{$k} eq ""){ $L{$k} = 0; }
    $total_sum += $L{$k};
    for (my $i=0; $i<=$#bin; $i++){
	if ($k >= $bin[$i]){ $S[$i]+=$L{$k};}
    }
    print OUT join("\t", $k, $L{$k}), "\n";
}
for (my $i=0; $i<=$#bin; $i++){
    print join("\t", $bin[$i], $S[$i], ($S[$i]/$total_sum)), "\n";
}
close(OUT);
exit;
