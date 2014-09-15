#!/usr/bin/perl

# select the unique sequences from sanger results 
# last modified: Sep-06-2014
# Developed by Jung-Ki Yoon

use strict;
my $usage = "./filter_uniq.pl [input: Sanger.fasta] [output]\n";
die $usage unless ($#ARGV == 1);
my ($in_fname, $out_fname,) = @ARGV;

open(IN, "<$in_fname") or die "[Error] Can't open $in_fname.\n";
my %D;
while(<IN>){
    my ($ID, ) = split /\s+/, $_;
    my ($seq, ) = split /\s+/, <IN>;
    if ($seq eq "EGFP") {next;}
    if (!(exists ($D{$seq}))) {
	$D{$seq} = $ID;
    }
}
close(IN);

open(OUT, ">$out_fname");
foreach my $k (keys %D){
    if ($k eq ""){next;}
    print OUT "$D{$k}\n$k\n";
}
close(OUT);
exit;
