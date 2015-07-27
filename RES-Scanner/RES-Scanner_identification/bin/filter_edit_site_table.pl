#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my ($input,$DNAdepth,$RNAdepth,$editLevel,$goodRead,$readType,$p,$help,$method,$Bayesian_Posterior_Probability,$P_value_DNA_Heterozygosis,$FDR_DNA_Heterozygosis,$Non_Ref_BaseCount,$Non_Ref_BaseRatio);
$DNAdepth ||= 7;
$RNAdepth ||= 3;
$editLevel ||= 0.05;
$goodRead ||= 1;
$readType ||= 3;
$p ||= 0.05;
$method ||= "Bayesian";
$Bayesian_Posterior_Probability ||= 0.95;
$P_value_DNA_Heterozygosis ||= 0.05;
$FDR_DNA_Heterozygosis ||= 0.05;
$Non_Ref_BaseCount ||= 0;
$Non_Ref_BaseRatio ||= 0;

GetOptions (
		"input:s"=>\$input,
		"DNAdepth:s"=>\$DNAdepth,
		"RNAdepth:s"=>\$RNAdepth,
		"editLevel:s"=>\$editLevel,
		"goodRead:s"=>\$goodRead,
		"readType:s"=>\$readType,
		"method:s"=>\$method,
		"Bayesian_P:s"=>\$Bayesian_Posterior_Probability,
		"Binomial_P:s"=>\$P_value_DNA_Heterozygosis,
		"Binomial_FDR:s"=>\$FDR_DNA_Heterozygosis,
		"Frequency_N:s"=>\$Non_Ref_BaseCount,
		"Frequency_R:s"=>\$Non_Ref_BaseRatio,
		"editPvalue:s"=>\$p,
		"help"=>\$help,
		);
if (!$input || $help) {
	print <<"Usage End.";
Description:
	This script is used to filter the RNA editing sites.
		--input             The input table of raw RNA edit sites.
		--DNAdepth          The homozygous DNA depth [optional(default:7)]
		--RNAdepth          The minimum cutoff for RNA reads coverage. [3]
		--editLevel         The minimum cutoff for editing level. [0.05]
		--method            Method for detecting SNPs.'Bayesian' or 'Binomial' or 'Frequency'. [Bayesian]
		--Bayesian_P        The minimun Bayesian Posterior Probability cutoff for corresponding genotype.(force --method Bayesian) [0.95]
		--Binomial_P        The maximun P value cutoff of Binomial test for DNA Heterozygosis. (force --method Binomial) [0.05]
		--Binomial_FDR      The maximun FDR cutoff of Binomial test for DNA Heterozygosis. (force --method Binomial) [0.05]
		--Frequency_N       The maximun non-refference base count cutoff. (force --method Frequency) [0]
		--Frequency_R       The maximun non-refference base ratio cutoff. (force --method Frequency) [0]
		--readType          The least number of RNA reads that were mapped to overlapping but not identical positions in the reference genome 
		                    that supporting an RNA-editing site in a given sample. [3]
		--goodRead          The least number of RNA reads which in the middle of their length supporting an editing site 
		                    (that is, from positions 26~75 of the 100-bp read), to avoid potential false positives resulting 
		                    from mis-mapping of reads at splice junctions. [1]
		--editPvalue        Cutoff of FDR RNA editing. [0.05]
		--help              Show the help information.


Example:
	perl $0 --input editSite.table --DNAdepth 7 --RNAdepth 3 --editLevel 0.05 --readType 3 --goodRead 1 --editPvalue 0.05

Usage End.
		exit;
}

die "$input does not exist!\n" unless -e $input;

if ($input =~ /\S+\.gz$/) {
	open (IN,"gunzip -c $input |") or die $!;
} else {
	open (IN,$input) or die $!;
}

while (<IN>) {
	chomp;
	if($_!~/^1.Chromosome_ID|^#/){
		my ($Chromosome_ID,$Coordinate,$strand,$gbase,$coverage_DNA,$DNA_BaseCount,$label1,$label2,$coverage_RNA,$RNA_BaseCount,$editType,$editing_degree,$badReadNum,$goodReadNum,$ReadTypeNum,$Pvalue,$Fdr) = split /\s+/,$_;
		next unless $coverage_DNA >= $DNAdepth && $coverage_RNA >= $RNAdepth && $editing_degree >= $editLevel && $goodReadNum >=$goodRead && $ReadTypeNum>=$readType && $Fdr<$p;
		if($method eq "Bayesian"){
			my $str=$label1;
			$str=~s/$gbase//g;
			next if $str;
			next if $label2<$Bayesian_Posterior_Probability;
		}elsif($method eq "Binomial"){
			next unless $label1<$P_value_DNA_Heterozygosis;
			next unless $label2<$FDR_DNA_Heterozygosis;
		}elsif($method eq "Frequency"){
			next if $label1>$Non_Ref_BaseCount;
			next if $label2>$Non_Ref_BaseRatio;
		}else{
			die "Error: unknown --Method $method\n";
		}
	}
	print "$_\n";
}
close IN;


