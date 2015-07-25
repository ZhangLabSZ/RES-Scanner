#!/usr/bin/perl -w
use strict;
die "<genome.fa> <sam2base>\n" unless @ARGV==2;

my $len;
if ($ARGV[0] =~ /\.gz$/) {
	open IN, "gunzip -c $ARGV[0] |" or die $!;;
} else {
	open IN, "$ARGV[0]" or die $!;;
}
$/ = ">";
<IN>;
while (<IN>) {
	chomp;
	s/.+\n//;
	s/\s+//g;
	$len += length($_);
}
$/="\n";
close IN;
my %hash;
my $covered_base=0;
my $total_depth=0;
if($ARGV[1]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[1] |" or die $!;
}else{
	open IN,"$ARGV[1]" or die $!;
}
while(<IN>){
	chomp;
	my @A=split /\t/;
	if($A[5]>0){
		$covered_base++;
	}
	$hash{$A[5]}++;
	$total_depth+=$A[5];
}
close IN;

my @array= sort {$hash{$b}<=>$hash{$a}} keys %hash;
my $peak_depth=$array[0];
my $average_depth=sprintf "%.2f",$total_depth/$covered_base;
my $coverage=sprintf "%.2f",$covered_base/$len*100;
print "#total_bases\t#covered_bases\tcoveage(%)\tcoverage(X)\tpeak_coverage(X)\n";
print "$len\t$covered_base\t$coverage\t$average_depth\t$peak_depth\n";






