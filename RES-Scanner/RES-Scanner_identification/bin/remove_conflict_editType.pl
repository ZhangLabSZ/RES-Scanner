#!/usr/bin/perl -w
use strict;
die "usage: perl $0 <RNAedit.merged.table>\n" unless @ARGV==1;
my @array;
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,$ARGV[0] or die $!;
}
while(<IN>){
	chomp;
	if($_=~/^Chromosome|^#/){
		print "$_\n";
		next;
	}
	my @A=split /\t/;
	push @array,[$A[0],$A[1],$_];
}
close IN;

@array=sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @array;

my ($chr_last,$pos_last,$info_last);
foreach my $value (@array){
	my ($chr,$pos,$info)=@$value;
	if(!defined $chr_last){
		($chr_last,$pos_last,$info_last)=($chr,$pos,$info);
	}else{
		if($chr eq $chr_last && $pos eq $pos_last){
			($chr_last,$pos_last,$info_last)=();
		}else{
			print "$info_last\n";
			($chr_last,$pos_last,$info_last)=($chr,$pos,$info);
		}
	}
}

if(defined $info_last){
	print "$info_last\n";
}
