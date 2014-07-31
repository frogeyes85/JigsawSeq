#!/usr/bin/perl
use strict;

# ATAAAGGAATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT    105738  4       C,T,G,A,        27204,26871,26613,25050,

my $usage = "./NTS_graph.pl [graph] [DP crit]\n";
if ($#ARGV != 1){die $usage;}
my ($in_fname, $crit) = @ARGV;
my @bin = qw(0 1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384 32768 65536 131078);
my @data;

open(IN, "<$in_fname");
while(<IN>){
    my @chunks = split /\s+/, $_;
    if ($chunks[2] <= 1){next;}
    if ($chunks[1] <= $crit) {next;}
    my @l = split(/\,/, $chunks[4]);
    @l = sort {$b <=> $a} @l;
    my $cur_ratio = $l[0]/$l[1];
    for(my $i=1; $i<=$#bin; $i++){
	if (($bin[$i-1]<$cur_ratio)&&($bin[$i]>=$cur_ratio)){
	    $data[$i]++;
	}
    }
}
close(IN);

for(my $i=0; $i<$#bin; $i++){
    print "$bin[$i]\t$data[$i]\n";
}

exit;
