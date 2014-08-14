#!/usr/bin/perl
use strict;

# last modified: Jul-04-2014
# Developed by Jung-Ki Yoon

my $usage = "./check_Sanger.pl [contigs1.fa] [contigs2.fa] [output]\n";
die $usage if ($#ARGV != 2);
my ($in_cont, $in_cont2, $out_fname,) = @ARGV;
my $num_con1=my $num_con2=my $num_union=0;
my $r1; my $r2;
my %Con;

open(IN, "<$in_cont") or die "Can't open $in_cont\n";
while(my $h=<IN>){ 
    substr($h, 0, 1) = "";
    chop($h);
    (my $seq,) = split /\s+/, <IN>;
    my @chunks = split /\s+/, $_;
    $Con{$seq}{'1'} = $h;
    $num_con1++;
}
close(IN);
open(IN2, "<$in_cont2") or die "Can't open $in_cont2\n";
while(my $h=<IN2>){ 
    substr($h, 0, 1) = "";
    chop($h);
    (my $seq, ) = split /\s+/, <IN2>;
    my @chunks = split /\s+/, $_;
    $Con{$seq}{'2'} = $h;
    $num_con2++;
}
close(IN2);
my $num=0;
open(OUT, ">$out_fname");
foreach my $k (keys %Con){
    if (($Con{$k}{'1'} ne "")&&($Con{$k}{'2'} ne "")) {$num_union++;}
    if ($Con{$k}{'1'} eq "") { 
	$r1 = "x";
    }else {
	$r1 = "o"; #r1 = $Con{$k}{'1'};
    }
    if ($Con{$k}{'2'} eq "") { 
	$r2 = "x";
    }else {
	$r2 = "o"; #r2 = $Con{$k}{'2'};
    }
    $num++;
    print OUT join("\t", $num, $r1, $r2, length($k), $k), "\n";
}
close(OUT);
print "Con1: $num_con1, Con2: $num_con2\nCon1_only: ", $num_con1-$num_union, ", Union: $num_union, Con2_only: ", $num_con2-$num_union, "\n";
exit;
