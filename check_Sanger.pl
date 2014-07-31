#!/usr/bin/perl
use strict;
# bwa mem -B30 -O30 -E10 k60s1_96_d6.contigs.fa ../intein_CFNerr96_Sanger.fasta | less
# Mapping pair-end read (Jigsaw-Seq data) into candidate contigs. 
# g12     0       k60s1_96_d12.contigs_55 60      60      420M    *       0       0       TGCCTTTCGTATGAGACGGAGATCCTTACTGTGGAATACGGTCTACTACCTATTGGTAAAATCGTGGAGAAACGCATCGAGTGCACTGTGTATTCAGTGGATAACAACGGGAATATCTATACCCAACCCGTGGCGCAGTGGCATGATCGGGGTGAGCAGGAGGTCTTCGAATACTGCCTTGAGGACGGGTCACTGATCCGCGCAACCAAAGATCACAAGTTCATGACGGTTGACGGACAAATGCTGCCAATTGACGAAATCTTCGAGCGAGAACTCGACCTAATGCGGGTCGATAATCTCCCGAACATTAAAATCGCCACGCGCAAATACCTAGGTAAACAAAATGTCTACGACATTGGTGTTGAACGCGATCACAACTTTGCTCTAAAGAATGGCTTCATCGCCTCTAATTGCTTTAAC   *       NM:i:1  AS:i:415        XS:i:54

# last modified: Jul-04-2014
# Developed by Jung-Ki Yoon

# contigs_1	1	1407
# contigs_1	2	1441
# contigs_1	3	1466
# contigs_1	4	1499
# contigs_1	5	1543

my $usage = "./check_Sanger.pl [contigs.fa] [Sanger.fa]\n";
if ($#ARGV != 1){die $usage;}
my ($in_cont, $in_Sanger,) = @ARGV;

system("./bwa mem $in_cont $in_Sanger > $in_cont\.Sanger.sam");

open(IN, "<$in_cont\.Sanger.sam");
my %San;
my %Con;
my %Eq;
my $num_contig=my $num_Sanger=my $num_eq=0;
while(<IN>){ 
    my @chunks = split /\s+/, $_;
    if (substr($_, 0, 1) eq "\@"){
	$num_contig++;
	substr($chunks[1], 0, 3) = "";
	substr($chunks[2], 0, 3) = "";
	$Con{$chunks[1]}->{'len'}=$chunks[2];
	next;
    }
    my ($dig,) = split /M/, $chunks[5];
    substr($chunks[12], 0, 5) = "";
    
    if ($dig eq $chunks[12]){
#	print "|$dig|$chunks[12]|$_";
	$Eq{$chunks[0]} = $chunks[2];
	$num_eq++;
    }
    if (!exists($San{$chunks[0]})){
	$San{$chunks[0]} = "Exist";
	$num_Sanger++;
    }
}
close(IN);
print "NumCon: $num_contig, NumSan: $num_Sanger, NumEq: $num_eq\n";

exit;
