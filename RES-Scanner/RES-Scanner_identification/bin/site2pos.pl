#!/usr/bin/perl -w
use strict;
die "Usage:perl $0 <RNAedtingType>\n" unless @ARGV == 1;

my %count;
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open (IN,$ARGV[0]) or die $!;
}
while (<IN>) {
	next if $_=~/^#|^1.Chromosome_ID|^Chromosome/;
	chomp;
	my @info = split /\t/;
	my $id = join "#",@info;
	$count{$id} ++;
	die "Something($info[0]\t$id) was wrong!\n" if ( $count{$id} >= 2);
	my $strand = $info[2];
	my $line = join "\t",$id,$info[0],$strand,$info[1],$info[1];
	print "$line\n";
}
close IN;


