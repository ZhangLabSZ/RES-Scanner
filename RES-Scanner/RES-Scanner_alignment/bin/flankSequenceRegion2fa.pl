#!/usr/bin/perl -w
use strict;
die "perl $0 <genome.fa> <cdsPOS2junction.coordinate>\n" unless @ARGV==2;

my %genome;
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,"$ARGV[0]" or die $!;
}
$/=">";
<IN>;
while(<IN>){
	chomp;
	$_=~/(.+?)\n/;
	my $id=(split /\s+/,$1)[0];
	$_=~s/.+?\n//;
	$_=~s/\s+//g;
	$genome{$id}=uc $_;
}
$/="\n";
close IN;


open IN,"$ARGV[1]" or die $!;

while(<IN>){
	chomp;
	next if $_=~/^#/;
	my @A=split /\t/;
	my @B=split /,/,$A[2];
	my $string;
	foreach my $region (@B){
		my @R=split /\-/,$region;
		my $len=$R[1]-$R[0]+1;
		my $str=substr($genome{$A[1]},$R[0]-1,$len);
		$string.=$str;
	}
	print ">$A[0]\t$A[1]\t$A[2]\n$string\n";
}
close IN;
