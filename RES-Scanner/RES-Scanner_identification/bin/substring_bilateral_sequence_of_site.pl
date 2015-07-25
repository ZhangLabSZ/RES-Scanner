#!/usr/bin/perl -w
use strict;
use Getopt::Long;
die "<RES_final_result.txt> <genome.fa>\n" unless @ARGV==2;
my $extend ||= 100;
GetOptions(
		"extend:s"=>\$extend,
		);
my %site;
open IN,"$ARGV[0]" or die $!;
while(<IN>){
	chomp;
	next if $_=~/^#|^Chromosome/;
	my @A=split /\t/;
	$site{$A[0]}{$A[1]}=1;
}
close IN;

if($ARGV[1]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[1] |" or die $!;
}else{
	open IN,"$ARGV[1]" or die $!;
}
$/=">";
<IN>;
while(<IN>){
	chomp;
	my $ref=$_;
	$ref=~/(.+?)\n/;
	my $id=(split /\s+/,$1)[0];
	if(exists $site{$id}){
		$ref=~s/.+?\n//;
		$ref=~s/\s+//g;
		my $len=length $_;
		foreach my $pos (sort keys %{$site{$id}}){
			my $beg=$pos-$extend<1 ? 1 : $pos-$extend;
			my $end=$pos+$extend>$len ? $len : $pos+$extend;
			my $str=substr($ref,$beg-1,$end-$beg+1);
			print ">$id;$pos;$beg;$end\n$str\n";
		}
	}
}
$/="\n";
close IN;






