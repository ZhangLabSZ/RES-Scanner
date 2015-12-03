#!/usr/bin/perl
use strict;
use Getopt::Long;
my (@fq,$outdir);
GetOptions(
		"fq:s"=>\@fq,
		"outdir:s"=>\$outdir,
		);
if (!@fq || !$outdir){die "perl $0 --fq IN.fq1  --fq IN.fq2  --outdir ./outDir\n"}
if(@fq>2){die "--fq error!\n";}
if ($fq[0] =~ /\.gz$/) {
	open IN1,"gunzip -c $fq[0] | " or die $!;
} else {
	open IN1,$fq[0] or die $!;
}
if(@fq==2){
	if ($fq[1] =~ /\.gz$/) {
		open IN2,"gunzip -c $fq[1] | " or die $!;
	} else {
		open IN2,$fq[1] or die $!;
	}
}
my $fq_num=1;
my $length=3221225472 ;# 1G =1024 **3   length=3G
my $fq_len=0;
system("mkdir -p $outdir/fq$fq_num/");
open OUT1,"| gzip >$outdir/fq$fq_num/$fq_num.1.fq.gz" or die $!;
if(@fq==2){
	open OUT2,"| gzip >$outdir/fq$fq_num/$fq_num.2.fq.gz" or die $!;
}
while(my $line1 = <IN1>){
	my $line2 = <IN1>;
	my $line3 = <IN1>;
	my $line4 = <IN1>;
	$fq_len+=length($line1.$line2.$line3.$line4);
	print OUT1 $line1, $line2, $line3, $line4;
	if(@fq==2){
		$line1 = <IN2>;
		$line2 = <IN2>;
		$line3 = <IN2>;
		$line4 = <IN2>;
		$fq_len+=length($line1.$line2.$line3.$line4);
		print OUT2 $line1, $line2, $line3, $line4;
	}
	if ($fq_len>=$length){
		close OUT1;
		if(@fq==2){
			close OUT2;
		}
		$fq_num++;
		$fq_len=0;
		system("mkdir -p $outdir/fq$fq_num/");
		open OUT1,"| gzip >$outdir/fq$fq_num/$fq_num.1.fq.gz" or die $!;
		if(@fq==2){
			open OUT2,"| gzip >$outdir/fq$fq_num/$fq_num.2.fq.gz" or die $!;
		}
	}
}
close IN1;
if(@fq==2){
	close IN2;
}
