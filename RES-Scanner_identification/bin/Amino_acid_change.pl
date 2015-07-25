#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <RE_final_result.annotation> <condon.database> \n" unless @ARGV == 2;

my %CODE = (
	'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
	'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
	'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
	'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
	'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
	'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
	'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
	'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
	'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
	'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
	'ATG' => 'M',                                                                         # Methionine
	'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
	'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
	'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
	'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
	'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
	'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
	'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
	'TGG' => 'W',                                                                         # Tryptophan
	'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
	'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'                                              # Stop
);

my %hash;
open (IN,$ARGV[0]) or die $!;
<IN>;
while (<IN>) {
	chomp;
	my @info = split /\s+/;
	$hash{$info[0]}{$info[1]}{$info[2]} = (split /->/,$info[4])[1];
}
close IN;

my %result;
open (IN,$ARGV[1]) or die $!;
while (<IN>) {
	chomp;
	my @info = split /\s+/;
	next unless (exists $hash{$info[0]}{$info[1]}{$info[2]});
	my @base = split //,$info[5];
	$base[$info[4]] = $hash{$info[0]}{$info[1]}{$info[2]};
	my $codon = join "",@base;
	my $C_c = join "->",$info[5],$codon;
	$CODE{$codon} = "*" unless (exists $CODE{$codon});
	my $amino = $CODE{$codon};
	my $C_a = join "$info[7]",$info[6],$amino;
#	my $key;
#	$key = 0 if ($info[6] eq $amino);
#	$key = 1 if ($info[6] ne $amino);
#	$key = -1 if ($amino eq "U" && $info[6] ne "U");
#	$key = 2 if ($amino ne "U" && $info[6] eq "U");
#	$result{$info[0]}{$info[1]}{$info[2]} = [$C_c,$C_a,$key];
	$result{$info[0]}{$info[1]}{$info[2]} = [$C_c,$C_a];
}
close IN;

open (IN,$ARGV[0]) or die $!;
my $title = <IN>;
chomp $title;
$title .= "\tCodonChange\tAminoAcidChange\n";
print "$title";
while (<IN>) {
	chomp;
	next if (/^#/);
	my @info = split /\s+/;
	if (exists $result{$info[0]}{$info[1]}{$info[2]}) {
		my $line = join "\t",@info,@{$result{$info[0]}{$info[1]}{$info[2]}};
		print "$line\n";
	} else {
		my $line = join "\t",@info,"-","-";
		print "$line\n";
	}
}
close IN;

