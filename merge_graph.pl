#!/usr/bin/perl

# This script is to merge de Bruijn graph files
# developed by Jung-Ki Yoon
# last modified: Jul-17-2014

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./merge_graph.pl [input: graph] [input:graph] [output: graph]\n";
die $usage unless ($#ARGV == 2);
my ($in_fname1, $in_fname2, $out_fname,) = @ARGV;

print "[Program:merge_graph] $in_fname1 + $in_fname2 -> $out_fname\n";

my $t_begin = new Benchmark;
my %data;
my $num_kmer=my $num_merged=0;
open(IN1, "<$in_fname1") or die "[Error] Can't open $in_fname1.\n";
while(<IN1>){
    my @foo = split /\s+/, $_;
    $data{$foo[0]} = join("\t", @foo[1..4]);
    $num_kmer++;
}
close(IN1);
print "[Report:merge_graph] $num_kmer kmer were loaded\n";

open(IN2, "<$in_fname2") or die "[Error] Can't open $in_fname2.\n";
while(<IN2>){
    my ($cur_kmer, @info2) = split /\s+/, $_;
    if (exists $data{$cur_kmer}) {
	# merge
	my %BackNode;
	my @info1 = split /\t/, $data{$cur_kmer};
	my @oriBackL = split /\,/, $info1[2];
	my @oriBackD = split /\,/, $info1[3];
	for(my $i=0; $i<=$#oriBackL; $i++){ $BackNode{$oriBackL[$i]} = $oriBackD[$i]; }

	my @curBackL = split /\,/, $info2[2];
	my @curBackD = split /\,/, $info2[3];
	for(my $i=0; $i<=$#curBackL; $i++){ $BackNode{$curBackL[$i]} += $curBackD[$i]; }

	$info1[2] = $info1[3] = "";
	$info1[1] = 0;
	foreach my $k (keys %BackNode){
	    $info1[1]++;
	    $info1[2] .= ($k . ",");
	    $info1[3] .= ($BackNode{$k} . ",");
	}
	$info1[0] += $info2[0];

	$data{$cur_kmer} = join("\t", @info1[0..3]);
	$num_merged++;
    }else{
	$data{$cur_kmer} = join("\t", @info2[0..3]);
	$num_kmer++;
    }
}
close(IN2);
my $t_end = new Benchmark; 
print "[Report:merge_graph] in2: $num_merged k-mer data were merged and $num_kmer k-mer were recorded.\n";
print "[Process:merge_graph] Merging was completed; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";

open(OUT, ">$out_fname");
foreach my $k (keys %data){
    print OUT join("\t", $k, $data{$k}), "\n";
}
close(OUT);

my $t_end = new Benchmark; 
print "[Process:merge_graph] Merging was completed; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
exit;
