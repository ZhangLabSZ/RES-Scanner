#!/usr/bin/perl -w
use strict;
use Getopt::Long;
die "perl $0 <overlap.dir> <RE_final_result.txt>\n" unless @ARGV == 2;
my %huge;
my @path = `find $ARGV[0]\/*.overlap`;

foreach my $p (@path) {
	my $ele = (split /\./,(split /-/,(split /\//,$p)[-1])[-1])[0];
	open (IN,$p) or die $!;
	while (<IN>){
		my @info = (split /\s+/);
		next unless ($info[6]);
		my ($name);
		my @B=split /#/,$info[0];
		my ($Chromosome_ID,$Coordinate)=($B[0],$B[1]);
		my $line = join "\t",@B;

		for (my $i=7;$i<@info;$i++) {
			die unless $info[$i] =~ /(\S+),\S,\d+,\d+/;
			$name = $1;	
			push @{$huge{$Chromosome_ID}{$Coordinate}{ele}} , $ele;
			push @{$huge{$Chromosome_ID}{$Coordinate}{name}} , $name;
		}
	}
	close IN;
}
my @output;

if($ARGV[1]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[1] |" or die $!;
}else{
	open (IN,$ARGV[1]) or die $!;
}
while (<IN>) {
	chomp;
	if($_=~/^#|^Chromosome/){
		print "$_\tTargetedGenomicFeature\tTargetedFeatureID\n";
		next;
	}
	my @info = split /\s+/;
	if (exists $huge{$info[0]}{$info[1]}){
		my $ele=join(";",@{$huge{$info[0]}{$info[1]}{ele}});
		my $name=join(";",@{$huge{$info[0]}{$info[1]}{name}});
		print "$_\t$ele\t$name\n",
	}else{
		print "$_\t-\t-\n";
	}
}
close IN;

