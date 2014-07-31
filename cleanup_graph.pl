#!/usr/bin/perl

# Clean up the de Bruijn graph to speed up the searching step.
# last modified: Jul-24-2014
# Developed by Jung-Ki Yoon

use strict;
use JigsawSeq;
use Benchmark ':hireswallclock';

my $usage = "[Usage] ./cleanup_graph.pl [input: graph] [r: cutoff ratio] [m: min node/edge depth] [output: graph]\n";
die $usage unless ($#ARGV == 3);
my ($in_fname, $cutoff_ratio, $minNodeDepth, $out_fname,) = @ARGV;
my $minEdgeDepth = $minNodeDepth;
my $polyA = "AAAAAAAA";   # due to illumina seq with short insertion size sequences
my $strange1 = "CATAAACAGTAATACAAGGGGTGTTGCGGAGTGTATACTGGCTTA";    
my $strange2 = "CATAAACAGTAATACAAGGGGTGTTACATAAACAGTAATACAAGGGGTGTT";    # due to Gibbs cloning breakpoint
substr($strange1, 0, 6) = substr($strange1, -6) = substr($strange2, 0, 8) = substr($strange2, -8) = "";

print "[Program:cleanup_graph] input: $in_fname cutoff_ratio: $cutoff_ratio min.Node/EdgeDepth: $minNodeDepth output:$out_fname\n";

my $t_begin = new Benchmark;
open(IN, "<$in_fname") or die "[Error] Can't open $in_fname.\n";
open(OUT, ">$out_fname");
my @edges;my @DPs;
my $num_cleaned=my $num_DP_skipped=my $num_polyA_skipped=my $num_Gibbs_skipped=0;
print OUT "#r: $cutoff_ratio\n";
while(<IN>){
    if (substr($_, 0, 1) eq "#") {
        print OUT $_;
        next;
    }
    chop($_);
    (my $seq, my $NodeDepth, my $num_edge, my $edge_str, my $DP_str,) = split(/\t/, $_);
    if ($NodeDepth<$minNodeDepth){
        $num_DP_skipped++;
        next;
    }
    if (($seq =~ $polyA)){
        $num_polyA_skipped++;
        next;
    }
    if (($seq =~ $strange1)||($seq =~ $strange2)){
        $num_Gibbs_skipped++;
        next;
    }

    @edges = split /\,/, $edge_str;
    @DPs = split /\,/, $DP_str;

    my $max=-1;
    for (my $i=0; $i<=$#DPs; $i++){
        $max = $DPs[$i] if ($max < $DPs[$i]);
    }

    my @new_DP_edges;
    my $new_edge_str = "";
    my $new_DP_str = "";

    for (my $i=0; $i<=$#DPs; $i++){
        if ( (($max/$DPs[$i])<=$cutoff_ratio) && ($DPs[$i]>=$minEdgeDepth) ){
            push @new_DP_edges, ($DPs[$i] . " " . $edges[$i]);
        }else{
            $num_cleaned++;
        }
    }
#    print "Cleaned\n";
#    print join("|", @new_DP_edges), " $#new_DP_edges\n";
    @new_DP_edges = sort { &csort($b, $a) } @new_DP_edges;
    $NodeDepth=0;
    for (my $i=0; $i<=$#new_DP_edges; $i++){
        my ($c_DP, $c_edge) = split /\s/, $new_DP_edges[$i];
        $new_edge_str .= ($c_edge . ",");
        $new_DP_str .= ($c_DP . ",");
        $NodeDepth += $c_DP;
    }
#    print "Sorted\n";
#    print join("|", @new_DP_edges), " $#new_DP_edges\n";
    if ($#new_DP_edges >= 0){
#        print join("\t", $seq, $NodeDepth, ($#new_DP_edges+1), $new_edge_str, $new_DP_str), "\n";
        print OUT join("\t", $seq, $NodeDepth, ($#new_DP_edges+1), $new_edge_str, $new_DP_str), "\n";
    }
#    <STDIN>;
}
close(OUT);
close(IN);
my $num_skipped = $num_DP_skipped + $num_polyA_skipped + $num_Gibbs_skipped;
my $t_end = new Benchmark;
print "[Report:cleanup_graph] Depth of node is less than $minNodeDepth: $num_DP_skipped.\n";
print "[Report:cleanup_graph] Poly-A node ($polyA): $num_polyA_skipped.\n";
print "[Report:cleanup_graph] Problematic vector node ($strange1, $strange2): $num_Gibbs_skipped.\n";
print "[Report:cleanup_graph] Total $num_skipped nodes were cleaned up.\n";
print "[Report:cleanup_graph] $num_cleaned edges were cleaned up (maxEdgeDepth/EdgeDepth > $cutoff_ratio or EdgeDepth is less than $minEdgeDepth.)\n";
print "[Process:cleanup_graph] Cleanup was completed; Processed Time = ", timestr(timediff($t_end, $t_begin)), "\n\n";
exit;

sub csort{
    my $ta = shift;
    my $tb = shift;
    (my $ca, ) = split /\s/, $ta;
    (my $cb, ) = split /\s/, $tb;
    $ca <=> $cb;
}
