#!/usr/bin/perl

package JigsawSeq;

#use Exporter;
#use vars qw(@ISA, @EXPORT_OK);
#@ISA = qw(Exporter);
#@EXPORT_OK = qw(rev_comp, codon2aa);

sub rev_comp($){
	my @in_str = split //, uc(shift);
	my %complemBase = ('A'=>'T', 'T'=>'A', 'C'=>'G', 'G'=>'C');
	my $re_str;
	for(my $i=0; $i<=$#in_str; $i++){
		return "Error" unless (exists $complemBase{$in_str[$i]});
		$re_str = $complemBase{$in_str[$i]} . $re_str;
	}
	return $re_str;
}

sub codon2aa($){
	my $codon = uc(shift);
	my %genetic_codes = ('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'*','TAG'=>'*','TGC'=>'C','TGT'=>'C','TGA'=>'*','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');;
	return "Error"	unless (exists $genetic_codes{$codon});
	return $genetic_codes{$codon};
}

sub call_sys($){
	my $sys_str = shift;
	print "[Run] $sys_str", "\n";
	system($sys_str);
}

sub memcheck{
	my $memory_used = `ps h -o size $$`;			#my $memory_used = `ps h -o sz $$`;
	chop($memory_used);
	$memory_used /= (1024*1024);
	$memory_used = substr($memory_used, 0, 4);
	return $memory_used;
}


#   memusage subroutine
#
#   usage: memusage [processid]
#
#   this subroutine takes only one parameter, the process id for 
#   which memory usage information is to be returned.  If 
#   undefined, the current process id is assumed.
#
#   Returns array of two values, raw process memory size and 
#   percentage memory utilisation, in this order.  Returns 
#   undefined if these values cannot be determined.
#sub memusage { # http://www.perlmonks.org/?node_id=115098
#    use Proc::ProcessTable;
#    my @results;
#    my $pid = (defined($_[0])) ? $_[0] : $$;
#    my $proc = Proc::ProcessTable->new;
#    my %fields = map { $_ => 1 } $proc->fields;
#    return undef unless exists $fields{'pid'};
#    foreach (@{$proc->table}) {
#        if ($_->pid eq $pid) {
#            push (@results, $_->size) if exists $fields{'size'};
#            push (@results, $_->pctmem) if exists $fields{'pctmem'};
#        };
#    };
#    return @results;
#}

1;