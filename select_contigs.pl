#!/usr/bin/perl
use strict;

# select real contigs from stat data
# last modified: Jul-31-2014
# requirement: perl version 5.10.1 or higher
# Developed by Jung-Ki Yoon

#k60s1_5k2.contigs_1     414     29.0496453900709        4.44561500605408        0.153035086878326       0.680335078679077
#k60s1_5k2.contigs_2     418     9.11538461538461        2.34833098572253        0.257622808560278       0.604983623970965

my $usage = "./select_contigs.pl [input: stat] [input: contigs.fa] [input: k-mer length] [input:step-size] [output_prefix: contigs]\n";
my $criteria = 0.2163;

die $usage if ($#ARGV != 4);
my ($in_stat, $in_fa, $in_kmer, $step_size, $out_fname,) = @ARGV;
my $out_stat = $out_fname . ".cv.stat";
my $out_fa = $out_fname . ".cv.fa";

open(IN, "<$in_stat") or die "Can't open $in_stat.\n";
open(OUT, ">$out_stat");
my %pass;
my $num_process=my $num_pass=0;
my $sum_ave=0;
while(<IN>){
    my ($name, $len, $ave, $std, $cv, $vmr,) = split /\s+/, $_;
    if ($cv <= $criteria){
	print OUT $_;
        $pass{$name} = "$ave\t$cv";
	$sum_ave += $ave;
        $num_pass++;
    }
    $num_process++;
}
close(IN);
close(OUT);
print $num_process, " contigs were processed, and ", $num_pass, " contigs were passed (cutoff: CV <= $criteria)\nAverage depth of passed contigs = ", ($sum_ave / $num_pass), "\n";

open(IN, "<$in_fa") or die "Can't open $in_fa.\n";
open(OUT, ">$out_fa");
$num_process=0;
my $last="";; 
while(my $header = <IN>){
    die "[Error] File format might be wrong\n" unless (substr($header, 0, 1) eq ">");
    substr($header, 0, 1)  = "";
    chop($header);
    $last = $header;
    (my $seq,) = split /\s+/, <IN>;

    substr($seq, 0, ($in_kmer - $step_size))="";
    substr($seq, -($in_kmer-$step_size)) = "";
    if (exists $pass{$header}){
        print OUT ">$header\t$pass{$header}\n$seq\n";
        $num_process++;
    }else{
#       die "Can't find $header in $in_fa.\n";
    }
}
print "Last: $last\n";
close(IN);
close(OUT);
exit;
