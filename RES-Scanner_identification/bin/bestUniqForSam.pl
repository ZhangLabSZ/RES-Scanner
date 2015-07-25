#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use File::Basename qw(basename);
use Getopt::Long;
die "<bam> <outdir> <samtools> --mis [int(readLength*0.06)] --ss [1|0] --mq [20] --soft [0|1]\n" unless @ARGV>=3;
my ($mis,$ss,$mq,$soft);
$ss ||= 1;
$mq ||= 20;
$soft ||= 0;
my $samtools=$ARGV[2];
die "$samtools does not exists!\n" unless -e $samtools;

GetOptions(
		"mis:s"=>\$mis,
		"ss:s"=>\$ss,
		"mq:s"=>\$mq,
		"soft:s"=>\$soft,
		);

if($ARGV[0]=~/\.sam.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}elsif($ARGV[0]=~/\.sam$/){
	open IN,"$ARGV[0]" or die $!;
}elsif($ARGV[0]=~/\.bam$/){
	open IN,"$samtools view -h $ARGV[0] |" or die $!;
}else{
	die "the suffix of input file $ARGV[0] should be .bam or .sam.gz/.sam!";
}
my $outdir = abs_path($ARGV[1]);
my $filename= basename($ARGV[0]);

if($ss){
	open OUT1,"| $samtools view -hbS - >$outdir/$filename.negative.bam" or die $!;
	open OUT2,"| $samtools view -hbS - >$outdir/$filename.positive.bam" or die $!;
}else{
	open OUT3,"| $samtools view -hbS - >$outdir/$filename.best.bam" or die $!;
}

	while(<IN>){
		chomp;
		my @array4tile=split /\t/;
		if(@array4tile<9){
			if($ss){
				print OUT1 "$_\n";
				print OUT2 "$_\n";
			}else{
				print OUT3 "$_\n";
			}
			next;
		}
		
		if( $_=~/XT:A:U\s+/ && $_=~/X0:i:1\s+/ && $_=~/X1:i:0\s+/ ){
			my @A=split /\t/;
			if($A[4]<$mq){next;}
			if($soft eq "1"){
				if($A[5]=~/D|I|S/){next;}
			}else{
				if($A[5]=~/D|I/){next;}                       ##### optional  
			}
			my ($missmatch)=$_=~/NM:i:(\d+)\s+/;
			my $seqLen=length($A[9]);
			if(!defined $mis){
				$mis=int $seqLen*0.06;
			}
			my $Tlen= $A[8] <0 ? -1*$A[8] : $A[8];
			if($missmatch>$mis || $Tlen>5000){next;}
			my $flag=$A[1];
			my %hash;
			&index($flag,\%hash);
			if(exists $hash{32} && exists $hash{16}){next;}
			if($ss){
				if(  exists $hash{64}  && exists $hash{32} && $A[8]>0  ){
					print OUT1 "$_\n";
				}elsif( exists $hash{128}  && exists $hash{16} && $A[8]<0 ){
					print OUT1 "$_\n";
				}elsif( exists $hash{64}  && exists $hash{16} && $A[8]<0 ){
					print OUT2 "$_\n";
				}elsif( exists $hash{128}  && exists $hash{32} && $A[8]>0){
					print OUT2 "$_\n";
				}elsif(exists $hash{64} && !exists $hash{16} && !exists $hash{32}){
					print OUT1 "$_\n";
				}elsif( exists $hash{128} && !exists $hash{16} && !exists $hash{32}){
					print OUT2 "$_\n";
				}elsif ($flag==16){
					print OUT2 "$_\n";
				}elsif ($flag==0){
					print OUT1 "$_\n";
				}
			}else{
				print OUT3 "$_\n";
			}
		}
	}
close IN;
if($ss){
	close OUT1;
	close OUT2;
}else{
	close OUT3;
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

