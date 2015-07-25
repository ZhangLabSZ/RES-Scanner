#!/usr/bin/perl -w
use strict;
die "Usage: <gtf file>\n" unless @ARGV == 1;

my %cdsPos;
if ($ARGV[0] =~ /\.gz$/) {
	open IN, "gunzip -c $ARGV[0] | ";
} else {
	open IN, $ARGV[0];
}
while (<IN>) {
	chomp;
	next if /^#/;
	my @info = split /\t/;
	if ($info[1] !~ "pseudogene" && $info[2] eq "CDS") {
		die "Warning: Unrecognized gff format, 'transcript_id' is missing." unless $info[-1] =~ /transcript_id\s+"(\S+)\s?";/;
		my $transcript_id = $1;
		my ($geneID)= $info[-1] =~ /gene_id\s+"(\S+)\s?";/;
		($info[3], $info[4]) = sort {$a <=> $b}($info[3], $info[4]);
		push @{$cdsPos{$transcript_id}}, [$info[3], $info[4], $info[0], $info[6], $info[7], $geneID];
	}
}
close IN;

foreach my $transcript (keys %cdsPos) {
	@{$cdsPos{$transcript}} = sort {$a->[0] <=> $b->[0]} @{$cdsPos{$transcript}};
	my ($transcript_bg, $transcript_ed, $chr, $strand, $geneID) = ($cdsPos{$transcript}->[0]->[0], $cdsPos{$transcript}->[-1]->[1], $cdsPos{$transcript}->[0]->[2], $cdsPos{$transcript}->[0]->[3],$cdsPos{$transcript}->[0]->[5]);
	print "$chr\tensembl\tmRNA\t$transcript_bg\t$transcript_ed\t.\t$strand\t.\tID=$transcript;GeneID=$geneID\n";
	foreach my $p (@{$cdsPos{$transcript}}) {
		my ($cds_bg, $cds_ed, $chr, $strand, $phase) = @$p;
		print "$chr\tensembl\tCDS\t$cds_bg\t$cds_ed\t.\t$strand\t$phase\tParent=$transcript;GeneID=$geneID\n";
	}
}
