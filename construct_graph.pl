#!/usr/bin/perl

# Construct de Bruijn graph from Jigsaw-Seq fastQ data. 
# last modified: Jul-15-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./construct_graph.pl [input: fastq format] [k: k-mer length] [s: step size] [output: graph]\n";
die $usage unless ($#ARGV == 3);
my ($in_fname, $in_kmer, $step_size, $out_fname,) = @ARGV;
die "[Error] k-mer length must be 9 <= k-mer <=93.\n" unless ((9<=$in_kmer)&&($in_kmer<=93));
die "[Error] step size must be 1, 2, or 3.\n" unless ((1<=$step_size)&&($3<=$step_size));
die "[Error] k-mer length must be divisible by step size.\n" if ($in_kmer % $step_size);

my $t_begin = new Benchmark;my $t_end;
my %node; my $a_node = \%node;
my $num_node=my $num_read=my $num_filter_node=0;
print "[Program: construct_graph] input file: $in_fname, k-mer length: $in_kmer, step size: $step_size, output file: $out_fname\n";

index_kmer($in_fname);

$t_end = new Benchmark;
print "[Process:construct_graph] (k,s) = ($in_kmer, $step_size) $num_read reads were processed; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
print "[Report:construct_graph] With k=$in_kmer and s=$step_size, $num_node nodes were detected.\n";

open(OUT, ">$out_fname");
for my $c_node (keys %node){
	unless ($num_filter_node % 5000000){
		$t_end = new Benchmark;
		print "[Process:construct_graph] $num_filter_node nodes were recorded; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
	}
#	next if ($node{$c_node}{'dep_node'}<$minNodeDepth);
	$num_filter_node++;	
	my $str_back_case="";
	my $str_back_case_depth="";
	my $num_edge=0;
	for my $c_edge (keys %{$a_node->{$c_node}}){
		next if (($c_edge eq "dep_node") || ($c_edge eq "num_edge"));
		$num_edge++; 
		$str_back_case .= ($c_edge . ",");
		$str_back_case_depth .= ($node{$c_node}{$c_edge} . ",");
	}
	print OUT join("\t", $c_node, $node{$c_node}{'dep_node'}, $num_edge, $str_back_case, $str_back_case_depth), "\n";
}
close(OUT);
print "[Report:construct_graph] With k=$in_kmer and s=$step_size, $num_filter_node nodes were recorded.\n";

$t_end = new Benchmark;
print "[Process:construct_graph] Contruction of de Bruijn graph was completed; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;

sub index_kmer(){
	my $fname = shift;
	open(IN, "<$fname") or die "[Error] Can't open $fname.\n";
	print "[Process:construct_graph] $fname was opened.\n";
	while(<IN>){          # 1st line: header
		$num_read++;
		unless ($num_read % 100000){
			$t_end = new Benchmark;
			print "[Process:construct_graph] (k,s) = ($in_kmer, $step_size) $num_read reads were processed; Mem. Used = ", JigsawSeq::memcheck(), " Gb; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
		}
		die "[Error] header was not detected. Probably input file was truncated or not fastq file.\n$_" unless (substr($_, 0, 1) eq "@");
		my $seq = <IN>;   # 2nd line; sequence
		chop($seq);
		<IN>;             # 3rd line; +
		<IN>;             # 4th line; score

		my $rev_seq = JigsawSeq::rev_comp($seq);
		next if ($rev_seq eq "Error"); # if the reverse complementary of this read could not be retrieved, skip the read.

		# Index K-mer and build de Bruijn graph
		for (my $i=0; $i<=(length($seq)-$in_kmer); $i++){
#		for (my $i=0; $i<=(length($seq)-$in_kmer); $i+=$step_size){
			my $cur_kmer = substr($seq, $i, $in_kmer);
			my $back_node = substr($cur_kmer, $step_size, $in_kmer-$step_size);
			my $back_case = substr($cur_kmer, -$step_size);
			my $front_node = substr($cur_kmer, 0, $in_kmer-$step_size);
			$num_node++ if ($node{$front_node}{'dep_node'} eq "");
			$node{$front_node}{'dep_node'}++;		
			$node{$front_node}{$back_case}++;

			# for reverse complementary sequences
			$cur_kmer = substr($rev_seq, $i, $in_kmer);	
			$back_node = substr($cur_kmer, $step_size, $in_kmer-$step_size);
			$back_case = substr($cur_kmer, -$step_size);
			$front_node = substr($cur_kmer, 0, $in_kmer-$step_size);
			$num_node++ if ($node{$front_node}{'dep_node'} eq "");
			$node{$front_node}{'dep_node'}++;		
			$node{$front_node}{$back_case}++;
		}
	}
	close(IN);
}
