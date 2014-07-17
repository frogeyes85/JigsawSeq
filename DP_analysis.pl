#!/usr/bin/perl
use strict;

# Mapping pair-end read (Jigsaw-Seq data) into candidate contigs. 
# last modified: Jul-16-2014
# Developed by Jung-Ki Yoon

# contigs_1	1	1407
# contigs_1	2	1441
# contigs_1	3	1466
# contigs_1	4	1499
# contigs_1	5	1543
my $len_read = 125;

my $usage = "./DP_analysis.pl [input:samtools depth] [input:contig length] [k: k-mer length] [s: step size] [output: stat]\n";

my ($in_fname, $len_fname, $in_kmer, $step_size, $out_fname) = @ARGV;
die $usage if ($#ARGV != 4);

my %len;
open(LEN, "<$len_fname") or die "Can't open $len_fname.\n";
while(<LEN>){
    my ($cur_con, $cur_len,) = split /\s+/, $_;
    $len{$cur_con} = $cur_len - ($in_kmer - $step_size) * 2;
}
close(LEN);

open(IN, "<$in_fname") or die "Can't open $in_fname.\n";
open(OUT, ">$out_fname");
my $prev_con="";
my %DP;
my $a_node=\%DP;
my $n=0;
my @DPs;
my $cal = ($len_read - ($in_kmer - $step_size))*2;
while(<IN>){
	my ($cur_con, $pos, $depth) = split /\s+/, $_;
#	print join("\t", ($in_kmer-$step_size), ($len{$cur_con}-($in_kmer-$step_size)), $cur_con, $pos, $depth);
	if (($pos <= ($len_read))||($pos > ($len{$cur_con}-(($len_read)-($in_kmer-$step_size)*2) ))) {
	    next;
#	}else{
#	    print "\n";
	}
   # remove initial/terminal seeds
	if ($cur_con ne $prev_con){
	    if ($prev_con ne ""){
		die "Can't find length data of $prev_con in $len_fname\n" unless (exists $len{$prev_con});
		for (my $i=0; $i<$len{$prev_con}-$cal-($#DPs+1); $i++){
		    push @DPs, 0;
#		    print $n++;
		}
#		print $n, "\n";
		print OUT join("\t", $prev_con, $len{$prev_con}, calc_stat(@DPs)), "\n";
#		print join(" ", @DPs), "\n";
#		exit;
		@DPs = qw();
	    }
	    $prev_con = $cur_con;
	}
	push @DPs, $depth;
}
close(IN);
for (my $i=0; $i<$len{$prev_con}-$cal-($#DPs+1); $i++){
    push @DPs, 0;
#    $n++;
}
print OUT join("\t", $prev_con, $len{$prev_con}, calc_stat(@DPs)), "\n";
#print join (" ", @DPs), "\n";
#print $n, "\n";
close(OUT);
exit;

sub calc_stat(){
    my @D = @_;
    my $c_sum=0;
    for (my $i=0; $i<=$#D; $i++){
	$c_sum += $D[$i];
    }
    my $c_ave = $c_sum / ($#D+1);
    my $c_std=0;
    for (my $i=0; $i<=$#D; $i++){
	$c_std += ($D[$i]-$c_ave)*($D[$i]-$c_ave);
    }
    $c_std = sqrt($c_std/($#D+1));
    
#	my $c_std = sqrt($c_sq_ave - $c_ave*$c_ave);
    return ($c_ave, $c_std, $c_std/$c_ave, $c_std*$c_std/$c_ave);
}
