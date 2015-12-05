#RES-Scanner -- RNA Editing Sites Scanner, a flexible and efficient software package that detects and annotates genome-wide RNA-editing sites using matched RNA-Seq and DNA-Seq data from the same individuals or samples.
#Copyright (C) 2015 China National GeneBank
#Zongji Wang <wangzongji@genomics.cn>
#Jinmin Lian <lianjinmin@genomics.cn>
#Qiye Li <liqiye@genomics.cn>
#Pei Zhang <zhangpei@genomics.cn>
#RES-Scanner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#RES-Scanner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my ($input,$DNAdepth,$RNAdepth,$editLevel,$p,$help,$method,$Bayesian_Posterior_Probability,$P_value_DNA_Heterozygosis,$FDR_DNA_Heterozygosis,$Non_Ref_BaseCount,$Non_Ref_BaseRatio,$extremeLevel);
$DNAdepth ||= 7;
$RNAdepth ||= 3;
$editLevel ||= 0.05;
$extremeLevel ||= 0;
$p ||= 0.05;
$method ||= "Bayesian";
$Bayesian_Posterior_Probability ||= 0.95;
$P_value_DNA_Heterozygosis ||= 0.05;
$FDR_DNA_Heterozygosis ||= 0.05;
$Non_Ref_BaseCount ||= 0;
$Non_Ref_BaseRatio ||= 0;
my $editDepth ||= 3;

GetOptions (
		"input:s"=>\$input,
		"DNAdepth:s"=>\$DNAdepth,
		"RNAdepth:s"=>\$RNAdepth,
		"editLevel:s"=>\$editLevel,
		"extremeLevel:s"=>\$extremeLevel,
		"method:s"=>\$method,
		"Bayesian_P:s"=>\$Bayesian_Posterior_Probability,
		"Binomial_P:s"=>\$P_value_DNA_Heterozygosis,
		"Binomial_FDR:s"=>\$FDR_DNA_Heterozygosis,
		"Frequency_N:s"=>\$Non_Ref_BaseCount,
		"Frequency_R:s"=>\$Non_Ref_BaseRatio,
		"editPvalue:s"=>\$p,
		"editDepth:s"=>\$editDepth,
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
		--extremeLevel      Exclude polymorphic sites with extreme degree of variation (100%) or not. 1 for yes, 0 for not. [0]
		--method            Method for detecting SNPs.'Bayesian' or 'Binomial' or 'Frequency'. [Bayesian]
		--Bayesian_P        The minimun Bayesian Posterior Probability cutoff for corresponding genotype.(force --method Bayesian) [0.95]
		--Binomial_P        The maximun P value cutoff of Binomial test for DNA Heterozygosis. (force --method Binomial) [0.05]
		--Binomial_FDR      The maximun FDR cutoff of Binomial test for DNA Heterozygosis. (force --method Binomial) [0.05]
		--Frequency_N       The maximun non-refference base count cutoff. (force --method Frequency) [0]
		--Frequency_R       The maximun non-refference base ratio cutoff. (force --method Frequency) [0]
		--editPvalue        Cutoff of FDR RNA editing. [0.05]
		--editDepth     INT     The minimum number of RNA reads supporting editing for a candidate editing site. [3]
		--help              Show the help information.


Example:
	perl $0 --input editSite.table --DNAdepth 7 --RNAdepth 3 --editLevel 0.05 --editPvalue 0.05

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
		my ($Chromosome_ID,$Coordinate,$strand,$gbase,$coverage_DNA,$DNA_BaseCount,$label1,$label2,$coverage_RNA,$RNA_BaseCount,$editType,$editing_degree,$Pvalue,$Fdr) = split /\s+/,$_;
		next unless $coverage_DNA >= $DNAdepth && $coverage_RNA >= $RNAdepth && $editing_degree >= $editLevel && $Fdr<$p;
		my @R=split /,/,$RNA_BaseCount;
		my %editHash=("A",$R[0],"C",$R[1],"G",$R[2],"T",$R[3]);
		my ($editBase)=$editType=~/\w\-\>(\w)/;
		$editBase=uc $editBase;
		if($strand eq "-"){
			$editBase=~tr/ATCGN/TAGCN/;
		}
		next unless $editHash{$editBase}>=$editDepth;
	
		if($extremeLevel){
			if($editHash{$gbase}==0 || $editing_degree==1){
				next;
			}
		}
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


