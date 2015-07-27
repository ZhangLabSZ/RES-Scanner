#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use File::Basename qw(basename dirname);
use Getopt::Long;
die "usage: perl $0 <abnormal.list> <inputBam> <outputBam> <samtools>\n" unless @ARGV==4;
my $inBamfile=$ARGV[1];
my $outBamfile=$ARGV[2];

my $filename=basename $outBamfile;
my $outdir=dirname $outBamfile;

my $samtools=$ARGV[3];
die "$samtools not exists!\n" unless -e $samtools;
$samtools=abs_path $samtools;


my %abnormal;
my $count=0;
my $tag=1;

open IN,"$ARGV[0]" or die $!;
while(<IN>){
	chomp;
	my $id;
	if($_=~/^@(.+?\/[12]);.+?/){
		($id)=$_=~/^@(.+?\/[12]);.+?/;
	}else{
		$id=(split /;/)[0];
	}
	$abnormal{$id}=1;
	$count++;

	if($count>=10000000){
		&filterBam(\%abnormal,$inBamfile,"$outdir/temp.$tag.$filename");
		$count=0;
		undef %abnormal;
		$inBamfile="$outdir/temp.$tag.$filename";
		if($tag>1){
			my $n=$tag-1;
			system "rm $outdir/temp.$n.$filename";
		}
		$tag++;
	}
}
close IN;


if($count>0){
	&filterBam(\%abnormal,$inBamfile,$outBamfile);
	if($tag>1){
		my $n=$tag-1;
		system "rm $outdir/temp.$n.$filename";
	}
}elsif($count==0 && $tag>1){
	my $n=$tag-1;
	system "mv $outdir/temp.$n.$filename $outBamfile";
}else{
#	print STDERR "Log: There is no ambiguous mapping between BWA and BLAT!\n";
	system "ln -sf $inBamfile $outBamfile";
}


sub filterBam{
	my ($abnormal_hash,$inputBam,$outputBam)=@_;
	open IN,"$samtools view -h $inputBam |" or die $!;
	open OUT,"| $samtools view -bhS - > $outputBam" or die $!;

	while(<IN>){
		chomp;
		my @A=split /\t/;
		if(@A<10){print OUT "$_\n";next;}
		my %hash;
		&index($A[1],\%hash);
		my $id;
		if(exists $hash{64}){
			$id="$A[0]\/1";
		}elsif(exists $hash{128}){
			$id="$A[0]\/2";
		}else{
			$id=$A[0];
		}
		if(!exists $$abnormal_hash{$id}){
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
}



sub index{
	my $num=shift;
	my $hash=shift;
	my $temp=1;
	if($num>=1){
		while($temp<=$num){
			$temp=$temp*2;
		}
		$temp=$temp/2;
		$$hash{$temp}=1;
		my $diff=$num-$temp;
		if($diff>0){
			&index($diff,\%$hash);
		}
	}
}

