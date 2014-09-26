#!/usr/bin/perl

use strict;
use JigsawSeq;

my $usage = "[Usage] ./chop_graph.pl [input: graph]\n";
die $usage unless ($#ARGV == 0);
my ($in_fname,) = @ARGV;
my %out;
my @L = qw(A T C G);
print "[Program:chop_graph] $in_fname -> ";
open(A, ">$in_fname.A");
open(T, ">$in_fname.T");
open(C, ">$in_fname.C");
open(G, ">$in_fname.G");

=pro
for (my $i=0; $i<=$#L; $i++){
    $out{$L[$i]} = $in_fname . ".$L[$i]";
    print "$out{$L[$i]}\t";
    open($L[$i], ">$out{$L[$i]}");
}
print "\n";
=cut

open(IN, "<$in_fname") or die "Can't open $in_fname\n";
while(<IN>){
    my $c = uc(substr($_, 0, 1));
    if ($c eq "A"){
	print A $_;
    }elsif ($c eq "T"){
	print T $_;
    }elsif ($c eq "C"){
	print C $_;
    }elsif ($c eq "G"){
	print G $_;
    } else{
	die "Error! $c\n$_";
    }
}
close(IN);
close(A);
close(T);
close(C);
close(G);
exit;
