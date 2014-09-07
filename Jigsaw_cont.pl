#!/usr/bin/perl

# Main module for analyzing JigsawSeq data 
# last modified: Jul-09-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./Jigsaw.pl [configuration file]\n";
my $cut_size = 12000000;

die $usage unless ($#ARGV == 0);
my $in_conf = shift;

my $t_begin = new Benchmark;

open(IN, "$in_conf") or die "[Error] Can't open $in_conf.\n";
my %conf;
while(<IN>){
	my ($a, $b,) = split /\s+/, $_;
	$conf{$a}=$b;
}
close(IN);

#input_F intein_CFNerr96_1.fastq
#input_R intein_CFNerr96_2.fastq
#vector_seq pBR322_vector.fasta
#exp_contig_size 420
#k-mer_len 30
#step_size 3
#min_depth 2
#cutoff_ratio 100
#output_prefix k30s3

die "[Error] Some required fields are missing in $in_conf.\n" unless ( (exists $conf{'input_F'})&&(exists $conf{'input_R'})&&(exists $conf{'vector_seq'})&&(exists $conf{'exp_contig_size'})&&(exists $conf{'k-mer_len'})&&(exists $conf{'step_size'})&&(exists $conf{'min_depth'})&&(exists $conf{'cutoff_ratio'})&&(exists $conf{'output_prefix'}) );
die "[Error] k-mer length must be 9 <= k-mer <=150.\n" unless ((9<=$conf{'k-mer_len'})&&($conf{'k-mer_len'}<=150));
die "[Error] step size must be 1, 2, or 3.\n" unless ((1<=$conf{'step_size'})&&($3<=$conf{'step_size'}));
die "[Error] k-mer length must be divisible by step size.\n" if ($conf{'k-mer_len'} % $conf{'step_size'});

#print "[Run:main] Start cutting input files.\n";
my $num_files=my $num_lines=0;
my $cur_fname;

=pro
open(IN, "<$conf{'input_F'}") or "[Error] Can't open $conf{'input_F'}";
while(<IN>){
	if (($num_lines % $cut_size) == 0){
		close(OUT) if ($num_files > 0);
		$cur_fname = "TMP_input_" . $conf{'output_prefix'} . ".$num_files" . ".fastq";
		open(OUT, ">$cur_fname");
		$num_files++;
	}
	print OUT $_;
	$num_lines++;
}
close(IN);
open(IN, "$conf{'input_R'}") or "[Error] Can't open $conf{'input_R'}";
while(<IN>){
	if (($num_lines % $cut_size) == 0){
		close(OUT);
		$cur_fname = "TMP_input_" . $conf{'output_prefix'} . ".$num_files" . ".fastq";
		open(OUT, ">$cur_fname");
		$num_files++;
	}
	print OUT $_;
	$num_lines++;
}
close(IN);
close(OUT);
print "[Process:main] $num_files temporary files were ready.\n";

$num_files = 3;
$cur_fname = "TMP_input_" . $conf{'output_prefix'};
for (my $i=0; $i<$num_files; $i++){
	JigsawSeq::call_sys("./construct_graph.pl $cur_fname\.$i\.fastq $conf{'k-mer_len'} $conf{'step_size'} $cur_fname\.$i\.graph");
}

if ($num_files > 1){
	JigsawSeq::call_sys("./merge_graph.pl $cur_fname\.0.graph $cur_fname\.1.graph $cur_fname\.0.merged");
	for (my $i=1; $i<$num_files-1; $i++){
		JigsawSeq::call_sys("./merge_graph.pl $cur_fname\." . ($i-1) . ".merged $cur_fname\." . ($i+1) . ".graph $cur_fname\.$i\.merged");
	}
	JigsawSeq::call_sys("mv $cur_fname\." . ($num_files-2) . ".merged $conf{'output_prefix'}\.graph");
} else{
	JigsawSeq::call_sys("mv $cur_fname.0.graph $conf{'output_prefix'}\.graph");
}
my $t_end=new Benchmark;
print "[Process:main] Merging was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
=cut

JigsawSeq::call_sys("./NTS_graph.pl $conf{'output_prefix'}\.graph 1000 > $conf{'output_prefix'}\.graph.NTS");
JigsawSeq::call_sys("./cleanup_graph.pl $conf{'output_prefix'}\.graph $conf{'cutoff_ratio'} $conf{'min_depth'} $conf{'output_prefix'}\.graph.clean");

JigsawSeq::call_sys("sort -r -n -k 2 $conf{'output_prefix'}\.graph.clean > $conf{'output_prefix'}\.graph.clean.sort");
my $t_end=new Benchmark;
print "[Process:main] Sorting was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";

JigsawSeq::call_sys("./Kmer2fa.pl $conf{'output_prefix'}\.graph.clean.sort $conf{'output_prefix'}\.kmer.fa");
JigsawSeq::call_sys("./detect_seeds.pl $conf{'output_prefix'}\.kmer.fa $conf{'vector_seq'} $conf{'k-mer_len'} $conf{'step_size'} $conf{'cutoff_ratio'} $conf{'output_prefix'}\.seeds");
JigsawSeq::call_sys("./explore_graph.pl $conf{'output_prefix'}\.graph.clean.sort $conf{'output_prefix'}\.seeds $conf{'k-mer_len'} $conf{'step_size'} $conf{'exp_contig_size'} $conf{'output_prefix'}\.contigs");
JigsawSeq::call_sys("./contigs2fa.pl $conf{'output_prefix'}\.contigs $conf{'output_prefix'}\.contigs.fa");
JigsawSeq::call_sys("./mapping_contigs.pl $conf{'output_prefix'}\.contigs.fa $conf{'input_F'} $conf{'input_R'} $conf{'k-mer_len'} $conf{'output_prefix'}\.contigs");
JigsawSeq::call_sys("./DP_analysis.pl $conf{'output_prefix'}\.contigs.DP $conf{'output_prefix'}\.contigs.len $conf{'k-mer_len'} $conf{'step_size'} $conf{'output_prefix'}\.contigs.stat");
JigsawSeq::call_sys("./select_contigs.pl $conf{'output_prefix'}\.contigs.stat $conf{'output_prefix'}\.contigs.fa $conf{'k-mer_len'} $conf{'step_size'} $conf{'output_prefix'}\.contigs");
JigsawSeq::call_sys("wc -l $conf{'input_F'} $conf{'input_R'} $conf{'output_prefix'}\.graph.clean");

my $t_end=new Benchmark;
print "[Process:main] Whole processes of JigsawSeq analysis were completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";

exit;
