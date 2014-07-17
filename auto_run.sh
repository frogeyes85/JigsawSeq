#/bin/bash

#./construct_graph.pl 5k/TMP_intein_5k_1.fastq.0 90 3 2 TMP_k90s3_5k.graph.0
#./construct_graph.pl 5k/TMP_intein_5k_1.fastq.1 90 3 2 TMP_k90s3_5k.graph.1
#./construct_graph.pl 5k/TMP_intein_5k_2.fastq.0 90 3 2 TMP_k90s3_5k.graph.2
#./construct_graph.pl 5k/TMP_intein_5k_2.fastq.1 90 3 2 TMP_k90s3_5k.graph.3

#./construct_graph.pl 5k/TMP_intein_5k_1.fastq.0 60 3 2 TMP_k60s3_5k.graph.0
#./construct_graph.pl 5k/TMP_intein_5k_1.fastq.1 60 3 2 TMP_k60s3_5k.graph.1
#./construct_graph.pl 5k/TMP_intein_5k_2.fastq.0 60 3 2 TMP_k60s3_5k.graph.2
#./construct_graph.pl 5k/TMP_intein_5k_2.fastq.1 60 3 2 TMP_k60s3_5k.graph.3

#./cleanup_graph.pl k90s1_5k.graph 100 2 k90s1_5k.graph.clean
#./Kmer2fa.pl k60s1_5k.graph.clean k60s1_5k.kmer.fa
#./detect_seeds.pl k60s1_5k.kmer.fa pBR322_vector.fasta 60 1 100 k60s1_5k.seeds
#./explore_graph.pl k60s1_5k.graph.clean k60s1_5k.seeds 60 1 420 k60s1_5k.contigs
#./contigs2fa.pl k60s1_5k.contigs k60s1_5k.contigs.fa
#./mapping_contigs.pl k60s1_5k.contigs.fa 5k/intein_5k_1.fastq 5k/intein_5k_2.fastq 60 k60s1_5k.contigs


#./cleanup_graph.pl k90s1_5k.graph 100 2 k90s1_5k.graph.clean
./Kmer2fa.pl k90s1_5k.graph.clean k90s1_5k.kmer.fa
./detect_seeds.pl k90s1_5k.kmer.fa pBR322_vector.fasta 90 1 100 k90s1_5k.seeds
./explore_graph.pl k90s1_5k.graph.clean k90s1_5k.seeds 90 1 420 k90s1_5k.contigs
./contigs2fa.pl k90s1_5k.contigs k90s1_5k.contigs.fa
./mapping_contigs.pl k90s1_5k.contigs.fa 5k/intein_5k_1.fastq 5k/intein_5k_2.fastq 90 k90s1_5k.contigs

