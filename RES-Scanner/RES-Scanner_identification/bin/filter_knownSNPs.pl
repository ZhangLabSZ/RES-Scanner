#!/usr/bin/perl -w
use strict;
die "perl $0 <edit.site> <knownSNPs.gff>\n" unless @ARGV==2;

my %hash;
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,"$ARGV[0]" or die $!;
}
while(<IN>){
	chomp;
	next if $_=~/^1.Chromosome_ID/;
	my @A=split /\t/;
	$hash{$A[0]}{$A[1]}=1;
}
close IN;

if($ARGV[1]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[1] |" or die $!;
}else{
	open IN,"$ARGV[1]" or die $!;
}
while(<IN>){
	chomp;
	next if $_=~/^#/;
	my @A=split /\t/;
	if(exists $hash{$A[0]}{$A[4]}){
		delete $hash{$A[0]}{$A[4]};
	}
}
close IN;

if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,"$ARGV[0]" or die $!;
}
while(<IN>){
	chomp;
	if($_=~/^1.Chromosome_ID/){
		print "$_\n";
		next;
	}
	my @A=split /\t/;
	if(exists $hash{$A[0]}{$A[1]}){
		print "$_\n";
	}
}
close IN;



