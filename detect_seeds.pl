#!/usr/bin/perl

# Detect initial/terminal seeds from K-mers.
# Requirement: bwa 0.7.9a-r786
# last modified: Jul-04-2014
# Developed by Jung-Ki Yoon

# bwa index -a is ref_seed.fa
# bwa mem -t 3 -O2 -E1 ref_seed.fa k90s3.fa

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./detect_seeds.pl [input: Kmer fasta] [input: vector fasta] [k: k-mer length] [s: step size] [r: cutoff ratio] [output: seeds]\n";
die $usage unless ($#ARGV == 5);
my ($in_fname, $vector_fname, $in_kmer, $step_size, $cutoff_ratio, $out_fname) = @ARGV;
die "[Error] k-mer length must be 9 <= k-mer <=150.\n" unless ((9<=$in_kmer)&&($in_kmer<=150));
die "[Error] step size must be 1, 2, or 3.\n" unless ((1<=$step_size)&&($3<=$step_size));
die "[Error] k-mer length must be divisible by step size.\n" if ($in_kmer % $step_size);
my $t_begin = new Benchmark;

my $len_ref_seed = int($in_kmer*1.5);				# allow a half size of k-mer as insertion on seeds
my $ref_seed_fname = "TMP_k$in_kmer\_" . $vector_fname;
my $seeds_sam_fname = "TMP_" . $in_fname . ".sam";

print "[Program:detect_seeds] input: $in_fname vector: $vector_fname k-mer_len: $in_kmer cutoff_ratio: $cutoff_ratio ref_seed_fa: $ref_seed_fname ref_seed_sam: $seeds_sam_fname len_ref_seed: $len_ref_seed output: $out_fname\n";

open(IN, "<$vector_fname") or die "[Error] Can't open $vector_fname.\n";
open(OUT, ">$ref_seed_fname");
<IN>;my ($seq,)=split /\s+/, <IN>;
print OUT ">initial\n", uc(substr($seq,-$len_ref_seed)), "\n";
print OUT ">terminal\n", uc(substr($seq,0,$len_ref_seed)), "\n";
close(OUT);
close(IN);

JigsawSeq::call_sys("./bwa index $ref_seed_fname");
if ($in_kmer > 40){    # ./bwa mem only work for kmer longer than xx.
	JigsawSeq::call_sys("./bwa mem -t 3 -O2 -E1 $ref_seed_fname $in_fname > $seeds_sam_fname");
} else{
	JigsawSeq::call_sys("./bwa aln -t 3 $ref_seed_fname $in_fname > TMP_$in_fname\.sai"); # may need further opimization for detecting seeds with indels.
	JigsawSeq::call_sys("./bwa samse $ref_seed_fname TMP_$in_fname\.sai $in_fname > $seeds_sam_fname");
}
my $t_end = new Benchmark; 
print "[Process:detect_seeds] Alignment was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n";
#my $seeds_sam_fname = "TMP_KanR_N_96.k60s3.kmer.fa.sam";

my $num_lines=my $num_align_init=my $num_init=my $num_align_term=my $num_term=0;
my $max_DP_init=my $max_DP_term=-1;
open(IN, "<$seeds_sam_fname") or die "[Error] Can't open $seeds_sam_fname.\n";
open(OUT, ">$out_fname");
while(<IN>){ 
	next if (substr($_, 0, 1) eq "@");
	$num_lines++;
	if ($num_lines % 1000000 == 0) {print "[Process:detect_seeds] $num_lines lines were processed.\n";}
	my ($qname, $flag, $rname, $pos, $mapQ, $CIGAR, $rnext, $pnext, $tlen, $str, $qual, $NM, $MD, $AS, $XS,) = split /\s+/, $_;
	next if ($flag != 0);
	my @d = split /[MIDNSHPX=]/, $CIGAR;
	my @l = split /[0-9]+/, $CIGAR;
	shift @l;
	my $sum_d=0;
	if ($rname eq "initial"){
		$num_align_init++;
		($qname, my $DP) = split /_/, $qname;
		if ($max_DP_init == -1){$max_DP_init=$DP;}
		next if ($max_DP_init / $DP > $cutoff_ratio);
		for(my $i=0; $i<=$#d; $i++){
			$sum_d+=$d[$i];
		}
		if (($pos + $sum_d) == ($len_ref_seed+1) ) { # find initial node!
		    if ($l[$#l] eq "S"){
			next;
		    }
		    print OUT ">initial $DP\n$str\n";
		    $num_init++;
		}
	}elsif ($rname eq "terminal"){
		$num_align_term++;
		($qname, my $DP) = split /_/, $qname;
		if ($max_DP_term == -1){$max_DP_term=$DP;}
		next if ($max_DP_term / $DP > $cutoff_ratio);
		for(my $i=0; $i<=$#d; $i++){
			if (($l[$i] eq "M")||($l[$i] eq "I")){
				$sum_d+=$d[$i];
			}
		}
		if (($pos == 1)&&($sum_d==($in_kmer-$step_size))){
			print OUT ">terminal $DP\n$str\n";
			$num_term++;
		}
	}
}
close(IN);
close(OUT);

print "[Report:detect_seeds] Cutoff ratio: $cutoff_ratio\n[Report] Max_inti_DP: $max_DP_init\tMax_term_DP: $max_DP_term\n";
print "[Report:detect_seeds] $num_lines K-mers were processed.\n";
print "[Report:detect_seeds] $num_align_init K-mers were alinged to initial node.\t$num_init were considered as initial nodes.\n";
print "[Report:detect_seeds] $num_align_term K-mers were aligned to terminal node.\t$num_term were considered as terminal nodes.\n";
my $t_end = new Benchmark; 
print "[Process:detect_seeds] Detection of seeds were completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;

