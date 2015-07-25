#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <*.fq.psl.gz> [outfile *.list]\n" unless @ARGV == 1;

my (%reads,%not);
if ($ARGV[0] =~ /\S+\.gz$/) {
	open (IN,"gunzip -c $ARGV[0] |") or die $!;
} else {
	open (IN,$ARGV[0]) or die $!;
}
for (my $i=0;$i<5;$i++) {
	<IN>;
}
while (<IN>) {
	my @info = split /\s+/;
	my @query = split /;/,$info[9];
	if(@query!=2){die "Error: ID contatins semicolon ';' !\n";}
	my ($scaffold,$st,$ed,$cigar,$strand,$mismatch) = split /:/,$query[1];
	next unless $info[0]>=($ed-$st+1)*0.9;
	my $head = $st-1;
	my $tail = $ed;
	if ($info[1] <= $mismatch) {
#		if (($info[8] ne $strand) || ($info[13] ne $scaffold)) {
		if ($info[13] ne $scaffold) {
			if (! exists $not{$query[0]}) {
				print "\@$info[9]\n";
				$not{$query[0]} = 1;
				next;
			}
		} else {
			if (($info[15] >= $tail) || ($info[16] <= $head)) {
				if (! exists $not{$query[0]}) {
					print "\@$info[9]\n";
					$not{$query[0]} = 1;
					next;
				}
			}
		}
	}	
}
close IN;

