#!/usr/bin/perl

# Mapping pair-end read (Jigsaw-Seq data) into candidate contigs. 
# Requirement: bwa 0.7.9a-r786     
#              samtools 0.1.19-44428cd
# last modified: Jul-04-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./mapping_contigs.pl [input: contigs fasta] [input: fastq format (forward)] [input: fastq format (reverse)] [k: k-mer length] [output: prefix of sam/bam]\n";
die $usage unless ($#ARGV == 4);
my ($in_contigs, $in_fnameF, $in_fnameR, $in_kmer, $out_fname,) = @ARGV;
die "[Error] k-mer length must be 9 <= k-mer <=150.\n" unless ((9<=$in_kmer)&&($in_kmer<=150));
my $sam_fname = "TMP_" . $out_fname . ".sam";
my $len_fname = $out_fname . ".len";

my $t_begin = new Benchmark;
JigsawSeq::call_sys("./bwa index $in_contigs");
JigsawSeq::call_sys("./bwa mem -t 3 -B100 -O100 -E50 $in_contigs $in_fnameF $in_fnameR > $sam_fname");

my $t_end = new Benchmark; 
print "[Process:mapping_contigs] Alignment was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";

open(IN, "<$sam_fname") or die "[Error] Can't open $sam_fname.\n";
open(OUT, ">$out_fname\.sam");
open(LEN, ">$len_fname");
my %len_contigs;
my $num_lines=my $num_EM=0;
while(<IN>){ 
	if (substr($_, 0, 1) eq "@"){
		print OUT $_;
		if (substr($_, 0, 3) eq "\@SQ"){        # retrieve length of (contig + seeds)
		    my ($temp, $name, $length,) = split /\s+/, $_;
		    substr($name, 0, 3) = "";
		    substr($length, 0 ,3 ) = "";
		    $len_contigs{$name}=$length;
		    print LEN "$name\t$length\n";
		}
		next;
	}
	$num_lines++;
	unless ($num_lines % 5000000){
		$t_end = new Benchmark;
		print "[Process:mapping_contigs] $num_lines lines were processed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
	}

	my ($qname, $flag, $rname, $pos, $mapQ, $CIGAR, $rnext, $pnext, $tlen, $str, $qual, $NM, $MD, $AS, $XS,) = split /\s+/, $_;
	my @d = split /[MIDNSHPX=]/, $CIGAR; # number in CIGAR
	my @l = split /[0-9]+/, $CIGAR;      # letter in CIGAR
	shift @l;
	substr($MD, 0, 5) = ""; 			 #MD:Z:150 -> 150

	# skip the fusion reads
	next unless (($rnext eq "=")||($rnext eq "*"));

	# exact match in contigs.
	if (($#l==0)&&($l[0] eq "M")&&($MD eq $d[0])) {
		print OUT $_;
		$num_EM++;
	}

	# exact match with soft clip at the begin or end of contigs.
#	if (($#l==1)&&($l[0] eq "S")&&($l[1] eq "M")){
#		if (($pos == 1)&&($MD eq $d[1])&&($d[1] >= $in_kmer)) {
#		if (($pos == 1)&&($MD eq $d[1])) {
#			print OUT $_;
#			next;
#		}
#	}
#	if (($#l==1)&&($l[0] eq "M")&&($l[1] eq "S")){
#		if (($pos == ($len_contigs{$rname}-$d[0]))&&($MD eq $d[0])&&($d[0] >= $in_kmer)){
#		if (($pos == ($len_contigs{$rname}-$d[0]))&&($MD eq $d[0])) {
#			print OUT $_;
#			next;
#		}
#	}

	# non-unique reads (XS = AS)
}
close(IN);
close(OUT);
close(OUT);
$t_end = new Benchmark; 
print "[Report:mapping_contigs] $num_EM reads were aligned to contigs exactly.\n";

JigsawSeq::call_sys("samtools view -bS $out_fname\.sam -o TMP_$out_fname\.bam");
JigsawSeq::call_sys("samtools sort TMP_$out_fname\.bam $out_fname");
JigsawSeq::call_sys("samtools index $out_fname\.bam");
JigsawSeq::call_sys("samtools depth $out_fname\.bam > $out_fname\.DP");

print "[Process:mapping_contigs] mapping raw data to contigs was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;
