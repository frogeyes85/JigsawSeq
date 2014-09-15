#!/usr/bin/perl

# run bwa and summarize the output of it
# last modified: Sep-06-2014
# Developed by Jung-Ki Yoon

use strict;
my $usage = "./compare_w_bwa.pl [input: JigsawSeq_contigs.cv.fa] [input: JigsawSeq_contigs.stat] [input: Sanger.uniq.fasta] [output: prefix]\n";
die $usage unless ($#ARGV == 3);
my ($in_Jig, $in_stat, $in_Sang, $out_fname,) = @ARGV;

`./bwa index $in_Jig`;
`./bwa mem $in_Jig $in_Sang > temp.result`;
#@SQ     SN:cardinal_N_96.k60s3.contigs_1        LN:727
#cardinal_N_66   0       cardinal_N_96.k60s3.contigs_542 1       60      735M    *       0       0       ATGGTTTCGAAAGGGGAAGAACTTATAAAGGAGAACATGCACATGAAACTTTACATGGAAGGCACCGTAAATAACCACCACTTCAAGTGTACCACTGAGGGTGAGGGGAAACCATACGAGGGGACCCAGACCCAGCGGATTAAAGTCGTTGAGGGCGGCCCCCTACCCTTTGCCTTTGATATCCTTGCGACGTGTTTTATGTATGGGTCCAAGACCTTCATAAATCACACGCAGGGTATCCCCGACTTCTTCAAGCAGTCTTTCCCGGAAGGGTTCACGTGGGAACGGGTCACGACTTACGAGGATGGCGGCGTCCTCACCGTGACACAGGATACCTCCCTCCAGGACGGCTGTCTCATATACAATGTAAAACTGCGAGGCGTCAACTTCCCTTCTAACGGGCCTGTGATGCAGAAGAAGACGCTAGGCTGGGAGGCGACTACTGAGACCCTGTATCCAGCGGACGGAGGACTGGAGGGCCGGTGCGATATGGCCCTCAAGCTTGTTGGCGGGGGGCATCTTCACTGCAATCTTAAAACCACGTATCGCTCAAAAAAACCCGCGAAAAACCTGAAAATGCCCGGCGTTTACTTTGTTGATCGTCGTCTTGAACGCATAAAGGAAGCTGATAATGAAACATATGTTGAACAACATGAGGTCGCCGTCGCTCGGTACTGCGACCTACCCTCGAAACTTGGGCACAAACTCAACGGAATGGACGAACTATACAAGTAA  *       NM:i:0  AS:i:735        XS:i:0
my %J;
my $num_S=my $num_J=my $num_Yes=0;
open(IN, "<$in_Jig") or die "[Error] Can't open $in_Jig\n";
while(<IN>){
    if (substr($_, 0, 1) eq ">"){
	my @chunks = split /\s+/, $_;
	substr($chunks[0], 0, 1) = "";
	$J{$chunks[0]}->{'DP'} = $chunks[1];
	$J{$chunks[0]}->{'CV'} = $chunks[2];
	$num_J++;
    }
}
close(IN);

open(IN, "<temp.result") or die "[Error] Can't open temp.result\n";
open(OUT, ">$out_fname\.Sanger");
my %S;
my $NM = my $CIGAR = my $S_ID = my $J_ID = "";
while(<IN>){    
    my @chunks = split /\s+/, $_;
    if ($chunks[0] eq "\@SQ"){
	substr($chunks[1], 0, 3) = "";
	substr($chunks[2], 0, 3) = "";
	$J{$chunks[1]}->{'len'} = $chunks[2];
	next;
    }
    $S_ID = $chunks[0]; 
    if ($S{$S_ID} == 1) {next;}
    $num_S++;
    $J_ID = $chunks[2];
    if ($J_ID eq "*"){
	$NM = $CIGAR = $chunks[5] = $chunks[11] = "n/a";
    }else{
	$CIGAR = $chunks[5]; 
	my @temp = split /M/, $CIGAR;
	$CIGAR = $temp[0];
	$NM = $chunks[11];
	my @temp = split /:/, $NM;
	if ($temp[0] ne "NM"){
	    die "Error: $_";
	}
	$NM = $temp[2];
    }

    if (($CIGAR == $J{$J_ID}->{'len'})&&($NM == 0)&&($CIGAR ne "n/a")){
	$J{$J_ID}->{'S_ID'} = $S_ID;
	print OUT join("\t", "Yes", $S_ID, $J_ID, $J{$J_ID}->{'DP'}, $J{$J_ID}->{'CV'}), "\n"; 
	$num_Yes++;
    }else{
	if ($NM==0){
	    print "[Report] please check this: $_";
	}
	print OUT join("\t", "No", $S_ID, $J_ID, $chunks[5], $chunks[11]), "\n";
    }
    $S{$S_ID} = 1;
}
close(IN);
print "[Report] #Sanger= $num_S #Jigsaw= $num_J #identical= $num_Yes\n";
close(OUT);
open(IN, "<$in_stat") or die "Can't open $in_stat\n";
open(OUT, ">$out_fname\.stat");
print OUT "#Sanger= $num_S #Jigsaw= $num_J #identical= $num_Yes\n";
while(<IN>){
    my @chunks = split /\s+/, $_;
    if ( $J{$chunks[0]}->{'S_ID'} eq "") {$J{$chunks[0]}->{'S_ID'} = "none";}
    print OUT join("\t", @chunks[0..4], $J{$chunks[0]}->{'S_ID'}), "\n";
}
close(IN);
close(OUT);
exit;
