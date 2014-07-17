#!/usr/bin/perl
use strict;

# select real contigs from stat data
# last modified: Jul-17-2014
# requirement: perl version 5.10.1 or higher
# Developed by Jung-Ki Yoon

#k60s1_5k2.contigs_1     414     29.0496453900709        4.44561500605408        0.153035086878326       0.680335078679077
#k60s1_5k2.contigs_2     418     9.11538461538461        2.34833098572253        0.257622808560278       0.604983623970965

my $usage = "./select_contigs.pl [input: stat] [input: contigs.fa] [output: contigs.fa]\n";
my $criteria = 0.306;

die $usage if ($#ARGV != 2);
my ($in_stat, $in_fa, $out_fname,) = @ARGV;

open(IN, "<$in_stat") or die "Can't open $in_stat.\n";
my %pass;
my $num_process=my $num_pass=0;
while(<IN>){ 
    my ($name, $len, $ave, $std, $cv, $vmr,) = split /\s+/, $_;
    if ($cv <= $criteria){
	$pass{$name} = "";
	$num_pass++;
    }
    $num_process++;
}
close(IN);
print $num_process, " contigs were processed, and ", $num_pass, " contigs were passed (cutoff: CV <= $criteria)\n";

open(IN, "<$in_fa") or die "Can't open $in_fa.\n";
open(OUT, ">$out_fname");
$num_process=0;
while(my $header = <IN>){
    die "[Error] File format might be wrong\n" unless (substr($header, 0, 1) eq ">");
    substr($header, 0, 1)  = "";
    chop($header);
    my $seq = <IN>;
    if (exists $pass{$header}){	
	print OUT ">$header\n$seq";
	$num_process++;
    }else{
#	die "Can't find $header in $in_fa.\n";
    }
} 
close(IN);
close(OUT);
if ($num_process < $num_pass){
    print "[Error] Can't find passed contigs in $in_fa.\n";
}else{
#    system("bwa index -a is $out_fname");
}
exit;
