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
use File::Basename qw(basename);
use Cwd 'abs_path';
use Getopt::Long;
my ($config,$phred,$qual_cutoff,$genome,$intron,$homopolymer,$paralogous_D);
my $HomoPrior ||= 0.99;
my $rate ||= 2; #the rate of transition over transversion
my $method ||= "Bayesian";
my $ploidy ||= 2;
my $DNAdepth ||= 10;
my $RNAdepth ||= 3;
my $Bayesian_Posterior_Probability ||= 0.95;
my $FDR_DNA_Heterozygosis ||= 0.05;
my $Non_Ref_BaseCount ||= 0;
$phred ||= "33,33";
$qual_cutoff ||= 30;
my $intronic ||=6;
$paralogous_D ||=1;
$homopolymer ||=1;

GetOptions(
		"config:s"=>\$config,
		"phred:s"=>\$phred,
		"qual_cutoff:s"=>\$qual_cutoff,
		"HomoPrior:s"=>\$HomoPrior,
		"rate:s"=>\$rate,
		"method:s"=>\$method,
		"genome:s"=>\$genome,
		"ploidy:s"=>\$ploidy,
		"intron:s"=>\$intron,
		"DNAdepth:s"=>\$DNAdepth,
		"RNAdepth:s"=>\$RNAdepth,
		"Binomial_FDR:s"=>\$FDR_DNA_Heterozygosis,
		"Bayesian_P:s"=>\$Bayesian_Posterior_Probability,
		"Frequency_N:s"=>\$Non_Ref_BaseCount,
		"homopolymer:s"=>\$homopolymer,
		"intronic:s"=>\$intronic,
		"paralogous_D:s"=>\$paralogous_D,
		);
if (!$config || !$genome) {
	print <<"Usage End.";
Description:
	This script is used to combine RNA editing sites for all samples.
		--config FILE       The tab-delimited configuration table which contains two columns each line:   
		                    column 1: sample name that used for REseeker pipeline.
		                    column 2: output directory of REscanner pipeline of the corresponding sample.
		--phred             The Phred base quality for query QUALity of DNA.bam and RNA.bam files respectively, default DNA_ASCII-33,RNA_ASCII-33. [33,33]
		--qual_cutoff       Quality cutoff for BWA alignment. [30]
		--method            Method for detecting SNPs.'Bayesian' or 'Binomial' or 'Frequency'. [Bayesian]
		--HomoPrior         The prior probability of homozygous genomic positions. (force --method Bayesian) [0.99]
		--rate              The rate of transitions over transversions of the genome. (force --method Bayesian) [2]
		--Bayesian_P        The minimun Bayesian Posterior Probability cutoff for corresponding genotype, range from 0 to 1. (force --method Bayesian) [0.95]
		--Binomial_FDR      The maximun FDR cutoff of Binomial test for DNA Heterozygosis, range from 0 to 1. (force --method Binomial) [0.05]
		--Frequency_N       The maximun non-reference base count cutoff. (force --method Frequency) [0]
		--genome FILE       The refference genome.
		--ploidy            The ploidy level, 1 for monoploid , 2 for diploid, 3 for triploid, 4 for tetraploid, and so on. [2].
		--intron FILE       The intron region file with POS format. This option is called for filter editing sites locating near junctions. [null]
		--intronic     NUM  Remove intronic candidate editing sites occurring within (the number of) bases of a splice site. [6]
		--DNAdepth          Coverage of DNA base, the site which is less than the coverage will not be marked as edit site in corresponding sample. [default 10]
		--RNAdepth          Coverage of RNA base, the site which is less than the coverage will not be marked as edit site in corresponding sample. [default 3]
		--homopolymer  NUM  Whether remove candidate editing sites in homopolymer runs of >= 5 base pairs. 1 for yes, 0 for not. [default 1]
		--paralogous_D NUM  Whether discard candidate editing sites with DNA reads depth of more than twice the genome-wide peak or mean depth. 
		                    1 for yes, 0 for not. [default 1]

		Example:		
		perl $0  --config <config.file> --genome genome.fa  
Usage End.
		exit;
}
if($method ne "Bayesian" && $method ne "Binomial" && $method ne "Frequency"){
	die "Error: unknown options value '$method' for --Method [Bayesian|Binomial|Frequency]\n";
}elsif($method eq "Bayesian" && !$genome){
	die "Error: options --genome is undefined in Bayesian mode\n";
}
my $HetePrior = 1-$HomoPrior;
my %hash;
my @order;
my ($phred_DNA,$phred_RNA)=split /,/,$phred;
open IN,"$config" or die $!;
while(<IN>){
	chomp;
	my ($sample,$directory)=split /\s+/,$_;
	push @order,$sample;
	die "$directory not exists!\n" unless -e $directory;
	$directory=abs_path($directory);
	my @sam2base_RNA=glob "$directory/*RNA.sam2base.gz";
	my @site=glob "$directory/*RNA.sam2base.homo.filter.gz";
	my @sam2base_DNA=glob "$directory/*DNA.sam2base.gz";
	my @sam2base_DNA_stat=glob "$directory/*DNA.sam2base.stat";
	$hash{$sample}{3}=$sam2base_DNA_stat[0];
	if(@sam2base_DNA==1){
		$hash{$sample}{0}=$sam2base_DNA[0];
	}elsif(@sam2base_DNA>1){
		die "Error: there are more than one DNA.same2base file in $directory!\n";
	}else{
		die "Error: DNA.same2base is not exist in $directory!\n";
	}
	if(@sam2base_RNA==2 && @site==2){
		my ($positive_sam2base,$negative_sam2base,$positive_site,$negative_site);
		foreach my $file (@sam2base_RNA){
			if($file=~/positive/){
				$positive_sam2base=$file;
			}elsif($file=~/negative/){
				$negative_sam2base=$file;
			}else{
				die "filename error: $file";
			}
		}
		foreach my $file (@site){
			if($file=~/positive/){
				$positive_site=$file;
			}elsif($file=~/negative/){
				$negative_site=$file;
			}else{
				die "filename error: $file";
			}
		}
		$hash{$sample}{1.1}=$positive_sam2base;
		$hash{$sample}{1.2}=$positive_site;
		$hash{$sample}{2.1}=$negative_sam2base;
		$hash{$sample}{2.2}=$negative_site;
	}elsif(@sam2base_RNA==1 && @site==1){
		my ($sam2base,$site);
		foreach my $file (@sam2base_RNA){
			if($file=~/best/){
				$sam2base=$file;
			}else{
				die "filename error: $file";
			}
		}
		foreach my $file (@site){
			if($file=~/best/){
				$site=$file;
			}else{
				die "filename error: $file";
			}
		}
		$hash{$sample}{1.1}=$sam2base;
		$hash{$sample}{1.2}=$site;
	}else{
		die "Error: number of files of $directory/*RNA.sam2base.gz and $directory/*RNA.sam2base.homo.filter.gz are not equal or the files are not exists!";
	}
}
close IN;
my (%table,%peak_depth);
my (%several_type);
foreach my $sample (keys %hash){
	if($hash{$sample}{1.2}=~/\.gz$/){
		open IN,"gunzip -c $hash{$sample}{1.2} |" or die $!;
	}else{
		open IN,"$hash{$sample}{1.2}" or die $!;
	}
	while(<IN>){
		chomp;
		next if $_=~/^1.Chromosome_ID/;
		my @A=split /\s+/;
		my ($key1,$key2,$key3,$key4,$key5)=("$A[0],$A[1],$A[2]","$sample,DNA","$sample,RNA","$sample,SNPvalue","$sample,Editvalue");
		$table{$key1}{$key2}=$A[5];
		$table{$key1}{$key3}=$A[9];
		if($method eq "Bayesian"){
			$table{$key1}{$key4}="$A[6],$A[7]";
		}else{
			$table{$key1}{$key4}=$A[6];
		}
		$table{$key1}{$key5}=$A[12];
		if(!exists $table{$key1}{gbase}){$table{$key1}{gbase}=uc $A[3];}
		if(!exists $table{$key1}{type}){
			$table{$key1}{type}= $A[10];
		}else{
			if($table{$key1}{type} ne $A[10]){
				$several_type{$key1}=1;
			}
		}
	}
	close IN;

	open IN,"$hash{$sample}{3}" or die $!;   # storage of peak coverage of DNA data.
		while(<IN>){
			chomp;
			next if $_=~/^#/;
			my @A=split /\s+/;
			$peak_depth{$sample}=$A[4]>$A[3] ? $A[4] : $A[3];
		}
	close IN;

	if(exists $hash{$sample}{2.2}){
		if($hash{$sample}{2.2}=~/\.gz$/){
			open IN,"gunzip -c $hash{$sample}{2.2} |" or die $!;
		}else{
			open IN,"$hash{$sample}{2.2}" or die $!;
		}
		while(<IN>){
			chomp;
			next if $_=~/^1.Chromosome_ID/;
			my @A=split /\s+/;
			my ($key1,$key2,$key3,$key4,$key5)=("$A[0],$A[1],$A[2]","$sample,DNA","$sample,RNA","$sample,SNPvalue","$sample,Editvalue");
			$table{$key1}{$key2}=$A[5];
			$table{$key1}{$key3}=$A[9];
			if($method eq "Bayesian"){
				$table{$key1}{$key4}="$A[6],$A[7]";
			}else{
				$table{$key1}{$key4}=$A[6];
			}
			$table{$key1}{$key5}=$A[12];
			if(!exists $table{$key1}{gbase}){$table{$key1}{gbase}= uc $A[3];}
			if(!exists $table{$key1}{type}){
				$table{$key1}{type}= $A[10];
			}else{
				if($table{$key1}{type} ne $A[10]){
					$several_type{$key1}=1;
				}
			}
		}
		close IN;
	}
}

############### filter sites locating near junctions BEGING ################
foreach my $key (keys %several_type){
	delete $table{$key};
}

if(defined $intron && -e $intron){
	if($intron=~/\.gz$/){
		open IN,"gunzip -c $intron |" or die $!;
	}else{
		open IN,"$intron" or die $!;
	}
	while(<IN>){
		chomp;
		my @A=split /\t/;
		my ($beg,$end)=sort {$a<=>$b} ($A[3],$A[4]);
		for(my $i=$beg;$i<=$beg+$intronic-1;$i++){
			my ($key1,$key2,$key3)=("$A[1],$i,+","$A[1],$i,-","$A[1],$i,.",);
			if(exists $table{$key1}){
				delete $table{$key1};
			}elsif(exists $table{$key2}){
				delete $table{$key2};
			}elsif(exists $table{$key3}){
				delete $table{$key3};
			}
		}
		for(my $i=$end-$intronic+1;$i<=$end;$i++){
			my ($key1,$key2,$key3)=("$A[1],$i,+","$A[1],$i,-","$A[1],$i,.",);
			if(exists $table{$key1}){
				delete $table{$key1};
			}elsif(exists $table{$key2}){
				delete $table{$key2};
			}elsif(exists $table{$key3}){
				delete $table{$key3};
			}		
		}
	}
	close IN;
}
############### filter sites locating near junctions END ################

####  filter homopolymer begin ###
if($homopolymer){
my %seq;
my %length;
if($genome=~/\.gz$/){
	open IN,"gunzip -c $genome |" or die $!;
}else{
	open IN, $genome or die $!;
}
$/ = ">";
<IN>;
while (<IN>) {
    /(.+)\n/;
	my $id = (split /\s+/, $1)[0];
	s/.+\n//;
	s/\s+|>//g;
	my $len = length($_);
	$seq{$id} = $_;
	$length{$id} = $len;
}
$/ = "\n";
close IN;

foreach my $key1 (keys %table){                       # length of homopolymer is longer than 5 nt.
	my ($chr,$pos,$strand)=(split /,/,$key1);
	my $beg=$pos-4>1 ? $pos-4 : 1;
	my $end=$pos+4<$length{$chr} ? $pos+4 : $length{$chr};
	my $nt=uc (substr($seq{$chr},$beg-1,$end-$beg+1));
    my $a_mer = "A"x5;
    my $t_mer = "T"x5;
    my $c_mer = "C"x5;
    my $g_mer = "G"x5;
    if ($nt =~ /$a_mer/i || $nt =~ /$t_mer/i || $nt =~ /$c_mer/i || $nt =~ /$g_mer/i) {
		delete $table{$key1};
	}
}
}
#### filter homopolymer end ###
my %BaseContent;
if($method eq "Bayesian"){
	&BasePercent($genome,\%BaseContent);
}

my $errorRate=10**(-1*$qual_cutoff/10);       # Q=-10log(10)P
foreach my $sample (keys %hash){
	open IN,"gunzip -c $hash{$sample}{0} |" or die $!;
	while(<IN>){
		chomp;
		my @A=split /\s+/;
		my $ref_base=uc $A[2];
		foreach my $strand ("+","-","."){
			my $key1="$A[0],$A[1],$strand";
			if(exists $table{$key1}){
				my $key2="$sample,DNA";
				if(!exists $table{$key1}{$key2}){
					my $key_SNPvalue="$sample,SNPvalue";
					if($A[5]==0){
						$table{$key1}{$key_SNPvalue}="NA";
						$table{$key1}{$key2}="0,0,0,0";
					}else{
						my ($dna_info,$seq,$qua) = &dealCoverage(@A,$phred_DNA);
						$table{$key1}{$key2}=$dna_info;
						if(!$seq){
							$table{$key1}{$key_SNPvalue}="NA";
							next;
						}
						my @result;
						if($method eq "Bayesian"){
							@result = &SNPvalue_Bayesian($seq,$qua,$ploidy,\%BaseContent);
							$table{$key1}{$key_SNPvalue}="$result[0],$result[1]";
						}elsif($method eq "Binomial"){
							@result = &SNPvalue_Binomial($dna_info,$ref_base,$ploidy);
							$table{$key1}{$key_SNPvalue}=$result[0];
						}elsif($method eq "Frequency"){
							@result = &SNPvalue_Frequency($seq,$ref_base);
							$table{$key1}{$key_SNPvalue}=$result[0];
						}else{
							die "Error: unknown --Method $method";
						}
					}
				}
			}
		}
	}
	close IN;
	open IN,"gunzip -c $hash{$sample}{1.1} |" or die $!;
	while(<IN>){
		chomp;
		my @A=split /\s+/;
		my $ref_base=uc $A[2];
		foreach my $strand ("+","."){
			my $key1="$A[0],$A[1],$strand";
			if(exists $table{$key1}){
				my $key_RNA="$sample,RNA";
				if(!exists $table{$key1}{$key_RNA}){
					my $key_Editvalue="$sample,Editvalue";
					if($A[5]==0){
						$table{$key1}{$key_RNA}="0,0,0,0";
						$table{$key1}{$key_Editvalue}="1";
					}else{
						my ($rna_info,$seq,$qua) = &dealCoverage(@A,$phred_RNA);
						$table{$key1}{$key_RNA}=$rna_info;
						$table{$key1}{$key_Editvalue}=&edit_Pvalue($ref_base,$rna_info,$errorRate);
					}
				}
			}
		}
	}
	close IN;
	if(exists $hash{$sample}{2.1}){
		open IN,"gunzip -c $hash{$sample}{2.1} |" or die $!;
		while(<IN>){
			chomp;
			my @A=split /\s+/;
			my $ref_base=uc $A[2];
			my $strand="-";
			my $key1="$A[0],$A[1],$strand";
			if(exists $table{$key1}){
				my $key_RNA="$sample,RNA";
				if(!exists $table{$key1}{$key_RNA}){
					my $key_Editvalue="$sample,Editvalue";
					if($A[5]==0){
						$table{$key1}{$key_RNA}="0,0,0,0";
						$table{$key1}{$key_Editvalue}="1";
					}else{
						my ($rna_info,$seq,$qua) = &dealCoverage(@A,$phred_RNA);
						$table{$key1}{$key_RNA}=$rna_info;
						$table{$key1}{$key_Editvalue}=&edit_Pvalue($ref_base,$rna_info,$errorRate);
					}
				}
			}
		}
		close IN;
	}
}
### 
foreach my $key1 ( keys %table){  
	foreach my $sample (@order){
		my ($key_DNA,$key_SNPvalue,$key_RNA,$key_Editvalue)=("$sample,DNA","$sample,SNPvalue","$sample,RNA","$sample,Editvalue");
		if(!exists $table{$key1}{$key_DNA}){
			$table{$key1}{$key_DNA}="0,0,0,0";
			$table{$key1}{$key_SNPvalue}="NA";
		}
		if(!exists $table{$key1}{$key_RNA}){
			$table{$key1}{$key_RNA}="0,0,0,0";
			$table{$key1}{$key_Editvalue}="1";
		}
	}
}
##                   
foreach my $key1 (keys %table){      
	my $flag=0;
	foreach my $sample (@order){
		my $key_DNA="$sample,DNA";
		my @A=split /,/,$table{$key1}{$key_DNA};
		if($paralogous_D){
		my $dep=$A[0]+$A[1]+$A[2]+$A[3];
		if($dep>2*$peak_depth{$sample}){                 #Filter unnormal high coverage sites, it may due to CNV.
			$flag=1;
			last;
		}
		}
		my $key_RNA="$sample,RNA";
		my @R=split /,/,$table{$key1}{$key_RNA};
		my %rna_count=("A",$R[0],"C",$R[1],"G",$R[2],"T",$R[3]);
		delete $rna_count{$table{$key1}{gbase}};
		my @rna_edit = sort {$b<=>$a} values %rna_count;
		if(@rna_edit!=3){die "$key1 error!";}
		if($rna_edit[0]>0 && $rna_edit[1]/$rna_edit[0]>0.01){     #filter mulitple RNA editing sites.
			$flag=1;
			last;
		}
	}
	if($flag){delete $table{$key1};}
}
foreach my $key1 (keys %table){       
	my @K=split /,/,$key1;
	my ($key_plus,$key_minus)=("$K[0],$K[1],+","$K[0],$K[1],-");
	if(exists $table{$key_plus} && exists $table{$key_minus}){ # Filter antisence transcript sites at double-strand transcription regions
		my ($plus_dep,$minus_dep)=(0,0);
		foreach my $sample (@order){
			my $key_RNA="$sample,RNA";
			my @plus=split /,/,$table{$key_plus}{$key_RNA};
			my @minus=split /,/,$table{$key_minus}{$key_RNA};
			$plus_dep+=$plus[0]+$plus[1]+$plus[2]+$plus[3];
			$minus_dep+=$minus[0]+$minus[1]+$minus[2]+$minus[3];
		}
		if($plus_dep>$minus_dep){
			delete $table{$key_minus};
		}elsif($plus_dep<$minus_dep){
			delete $table{$key_plus};
		}
	}
}
####
my (%Fdr_array,%Binomial_array);
foreach my $key1 (keys %table){
	foreach my $sample (@order){

		my ($chr,$pos,$strand)=split /,/,$key1;
		my $key_Editvalue="$sample,Editvalue";

		push @{$Fdr_array{$sample}},[$chr,$pos,$strand,$table{$key1}{$key_Editvalue}];
		if($method eq "Binomial"){
			my $key_SNPvalue="$sample,SNPvalue";
			push @{$Binomial_array{$sample}},[$chr,$pos,$strand,$table{$key1}{$key_SNPvalue}];
		}
	}
}
foreach my $sample (keys %Fdr_array){
	my @A=@{$Fdr_array{$sample}};
	@A=sort {$b->[3]<=>$a->[3]} @A;
	my @B;
	for(my $i=0;$i<@A;$i++){
		push @B,$A[$i][3];	
	}
	my @B_Fdr= &Fdr(\@B);
	for(my $i=0;$i<@A;$i++){
		my $key1="$A[$i][0],$A[$i][1],$A[$i][2]";
		my $key_Editvalue="$sample,Editvalue";
		$table{$key1}{$key_Editvalue}=&formatP($B_Fdr[$i]);
	}
}

my $least_depth=1;
if($method eq "Binomial"){
	foreach my $sample (keys %Binomial_array){
		my @A=@{$Binomial_array{$sample}};
		@A=sort {$b->[3]<=>$a->[3]} @A;
		my @B;
		for(my $i=0;$i<@A;$i++){
			push @B,$A[$i][3];	
		}
		my @B_Fdr= &Fdr(\@B);
		for(my $i=0;$i<@A;$i++){
			my ($key1,$key_SNPvalue)=("$A[$i][0],$A[$i][1],$A[$i][2]","$sample,SNPvalue");
			$table{$key1}{$key_SNPvalue}=&formatP($B_Fdr[$i]);
		}
	}

	my @p=sort {$a<=>$b} values %peak_depth;         #evalueate the least coverage.
	while($least_depth<$p[0]){
		my $ratio=1/$ploidy;
		my $lower_tail_p = &pbinom(0, $least_depth, $ratio, "lower");	
		if($lower_tail_p<0.05){last;}
		$least_depth++;
	}
}
##
my @title=("Chromosome","Coordinate","Strand","Gbase","EditType");
foreach my $sample (@order){
	push @title,"$sample.DNA_baseCount[A,C,G,T]","$sample.RNA_baseCount[A,C,G,T];P_value";
}
print join("\t",@title),"\n";
foreach my $key1 (sort keys %table){
	my $base=$table{$key1}{gbase};
	my @info;
	my $flag=1;
	foreach my $sample (@order){
		my ($key_DNA,$key_RNA,$key_SNPvalue,$key_Editvalue)=("$sample,DNA","$sample,RNA","$sample,SNPvalue","$sample,Editvalue");
		my @dna_coverage=split /,/,$table{$key1}{$key_DNA};
		my $total_dna_dep=$dna_coverage[0]+$dna_coverage[1]+$dna_coverage[2]+$dna_coverage[3];
		my @rna_coverage=split /,/,$table{$key1}{$key_RNA};
		my $total_rna_dep=$rna_coverage[0]+$rna_coverage[1]+$rna_coverage[2]+$rna_coverage[3];
		if($table{$key1}{$key_SNPvalue} ne "NA"){
			if($method eq "Bayesian"){
				my ($str,$p)=split /,/,$table{$key1}{$key_SNPvalue};
				$str=~s/$base//g;
				if($str || $p<$Bayesian_Posterior_Probability){$flag=0;}
			}elsif($method eq "Binomial"){
				if($total_dna_dep<$least_depth){
					my $n=0;
					foreach my $d (@dna_coverage){
						if($d==0){$n++;}
					}
					if($n<3){$flag=0;}
				}else{
					if($table{$key1}{$key_SNPvalue}<$FDR_DNA_Heterozygosis){     
						$flag=0;
					}
				}
			}elsif($method eq "Frequency"){
				if($table{$key1}{$key_SNPvalue}>$Non_Ref_BaseCount){
					$flag=0;
				}
			}else{
				die "Error: unknown --Method $method";
			}
		}
		if($total_dna_dep>=$DNAdepth && $total_rna_dep>=$RNAdepth && $table{$key1}{$key_Editvalue}<0.05){
        	my ($alt_base)=$table{$key1}{type}=~/\w->(\w)/;
			my $strand=(split /,/,$key1)[2];
			if($strand eq "-"){
				$alt_base=~tr/ATGC/TACG/;
			}
			my @dep=split /\,/,$table{$key1}{$key_RNA};
			my %baseDepth=("A",$dep[0],"C",$dep[1],"G",$dep[2],"T",$dep[3]);			
			if($baseDepth{$alt_base}>0){
				push @info,$table{$key1}{$key_DNA},"$table{$key1}{$key_RNA};$table{$key1}{$key_Editvalue}*";
			}else{
				push @info,$table{$key1}{$key_DNA},"$table{$key1}{$key_RNA};$table{$key1}{$key_Editvalue}";
			}
		}else{
			push @info,$table{$key1}{$key_DNA},"$table{$key1}{$key_RNA};$table{$key1}{$key_Editvalue}";
		}
	}
	if($flag){
		my ($chr,$pos,$strand)=split /,/,$key1;
		print "$chr\t$pos\t$strand\t$base\t$table{$key1}{type}\t",join("\t",@info),"\n";
	}
}
###################################################################################################
sub edit_Pvalue {
	my ($ref,$info,$errorRate)=@_;
	my @I=split /,/,$info;
	my %h=('A'=>$I[0],'C'=>$I[1],'G'=>$I[2],'T'=>$I[3]);
	my ($editDep,$totalDep)=0;
	foreach my $base (keys %h){
		$totalDep+=$h{$base};
		if($base ne $ref && $h{$base}>$editDep){
			$editDep=$h{$base};
		}
	}
	my $upper_tail_p = &pbinom($editDep, $totalDep, $errorRate, "upper");
	return $upper_tail_p;
}
sub dealCoverage {
	my @array = @_;
	die "These bases($array[0]\t$array[1]) doesn't match!\n" unless ($array[3] =~ /^\[M\]:([ATCGNatcgn0]+)(;\[I\]:[ATCGNatcgn,]+)?(;\[D\]:[ATCGNatcgn,]+)?/);
	my @bases = split //,$1;
	my $number = $#bases+1;
	my @qualities;
	my @aa = split //,$array[4];
	for (my $k=4;$k<$number+4;$k++) {
		push @qualities,$aa[$k];
	}
	die "Something was wrong at the location[$array[0]\t$array[1]]!\n" unless @bases == @qualities;
	my $phred=$array[-1];
	my (@b_temp,@q_temp);
	for (my $i=0;$i<@bases;$i++) {
		my $score = ord($qualities[$i]) - $phred;
		if ($score >= $qual_cutoff) {
			push @b_temp,$bases[$i];
			push @q_temp,$qualities[$i];
		}
	}
	my $b_temp = join "",@b_temp;
	my $q_temp = join "",@q_temp;

	my %seq_hash;
	my @seq_array=split //,$b_temp;
	foreach my $seq_base (@seq_array){
		$seq_hash{$seq_base}++;
	}
	my @seq_info;
	if(!exists $seq_hash{'A'}){$seq_hash{'A'}=0;}
	if(!exists $seq_hash{'C'}){$seq_hash{'C'}=0;}
	if(!exists $seq_hash{'G'}){$seq_hash{'G'}=0;}
	if(!exists $seq_hash{'T'}){$seq_hash{'T'}=0;}
	push @seq_info,$seq_hash{'A'},$seq_hash{'C'},$seq_hash{'G'},$seq_hash{'T'};
	my $seq_info=join(",",@seq_info);
	return ($seq_info,$b_temp,$q_temp);
}
###############################
sub SNPvalue_Bayesian{
#       my $HomoPrior = 0.99;
#       my $HetePrior = 0.01;
#       my $rate = 4; #the transition-transversion rate ratio
	my $weight = 0.5; #the weight of A/C and G/T;
	my $weight_other = (1-$weight)/2; #the weight of the other substitutions
		my @pileup = split //,$_[0];
	my @quals = split //,$_[1];
	my $ploid = $_[2];
	my $BaseContent=$_[3];
	my @Type = ('A','T','G','C');
	my %FixedError = (
			'A'=>{
			'C'=>$weight,
			'T'=>$weight_other,
			'G'=>$weight_other,
			},
			'C'=>{
			'A'=>$weight,
			'T'=>$weight_other,
			'G'=>$weight_other,
			},
			'G'=>{
			'T'=>$weight,
			'A'=>$weight_other,
			'C'=>$weight_other,
			},
			'T'=>{
			'G'=>$weight,
			'A'=>$weight_other,
			'C'=>$weight_other,
			},
			);

	my %P_allele;   #The probability of a given allele
		foreach my $i (0..$#pileup) {
			next unless ($pileup[$i] =~ /[ATGC]/);
			my $score = &QualityProbability($quals[$i]);
			foreach my $base (@Type) {
				if ($base eq $pileup[$i]) {
					push @{$P_allele{$base}}, 1-$score
				} else {
					push @{$P_allele{$base}}, $score * $FixedError{$base}{$pileup[$i]}};
			}
		}
	my (%posterior, %prior, %TypePrior);
	foreach my $base (@Type) {
		my $genotype = $base x $ploid;
		my $value=1;
		foreach my $score (@{$P_allele{$base}}) {
			$value += log($score);
		}
		$prior{$genotype} = $value;     #The probability of pileup data given the genotype (homozygous): p(D|G);
		$TypePrior{$genotype} = $HomoPrior * $$BaseContent{$base};      #the prior probability of each genotype (homozygous): p(G);
	}
	my $HeteNum = ($ploid-1)*4+($ploid-1)*2*$rate;
	my %control;
	foreach my $key (keys %FixedError){
		foreach my $ke (keys %{$FixedError{$key}}){                                # $key: original base; $ke: observed base;
			next if (exists $control{$key}{$ke});
			next if ($ploid == 1);
			$control{$key}{$ke} ++;
			$control{$ke}{$key} ++;
			for(my $i=1;$i<$ploid;$i++){
				my ($genotype,$value) = ("",1);
				$genotype .= $key x $i;
				$genotype .= $ke x ($ploid-$i);

				for (my $j=0;$j<=$#{$P_allele{$key}};$j++){
					$value += log( ($i*$P_allele{$key}[$j] + ($ploid-$i)*$P_allele{$ke}[$j]) / $ploid);
				}
				$prior{$genotype} = $value;     #The probability of pileup data given the genotype (heterozygous): p(D|G);

				if (($key eq "G" && $ke eq "A") || ($key eq "A" && $ke eq "G") || ($key eq "T" && $ke eq "C") || ($key eq "C" && $ke eq "T")) {
					$TypePrior{$genotype} = $HetePrior/$HeteNum*$rate; #the prior probability of each genotype (heterozygous): p(G) [genotypes including the allele pair of A/G or T/C];
				} else {
					$TypePrior{$genotype} = $HetePrior/$HeteNum; #the prior probability of each genotype (heterozygous): p(G);
				}
			}
		}
	}

	my $P_pileup = 0;
	my @greatest = sort {$prior{$b}<=>$prior{$a}} keys %prior;
	my $max = $prior{$greatest[0]};
	foreach my $genotype (keys %prior){
		$prior{$genotype} -= $max;
		my $P_genotype; #the prior probability of corresponding genotype;
		if (exists $TypePrior{$genotype}){
			$P_genotype = $TypePrior{$genotype};
		} else {
			die "There is not the prior probability for $genotype\n";
		}
		$P_pileup += $P_genotype * exp($prior{$genotype});      #The probability of pileup data: p(D)
	}

	foreach my $genotype (keys %prior){
		my $P_genotype;
		if (exists $TypePrior{$genotype}){
			$P_genotype = $TypePrior{$genotype};
		} else {
			die "There is not the prior probability for $genotype\n";
		}
		$posterior{$genotype} = $P_genotype * exp($prior{$genotype}) / $P_pileup;
	}

	my @sort = sort{ $posterior{$b}<=>$posterior{$a} }keys %posterior;
	my $Type = &CheckHomo($sort[0]);

#       return($Type,$posterior{$sort[0]},$sort[0]);
	return($sort[0],$posterior{$sort[0]});
}

sub QualityProbability{
	my$score=shift;
	my$value=ord($score)-$phred_DNA;
	my$p=10**(0-$value/10);
#       $p=$p/(1+$p);  # for solexa
	return($p);
}

sub CheckHomo{
	my $string = $_[0];
	my $base = (split //, $string)[0];
	$string =~ s/$base//g;
	my $zygous;
	if ($string eq "") {
		$zygous = "HOMO";
	} else {
		$zygous = "HETE";
	}
	return $zygous;
}

sub BasePercent {
	my ($file,$ref) = @_;
	if ($file =~ /\.gz$/) {
		open (INN,"gunzip -c $file |") or die $!;
	} else {
		open (INN,$file) or die $!;
	}
	$/ = ">";
	<INN>;
	my ($total,$A_num,$T_num,$C_num,$G_num);
	while (<INN>) {
		chomp;
		my $scaf;
		if (/(\S+)/) {
			$scaf = $1;
		}
		s/.+\n//;
		s/\s+//g;
		$_ =~ s/N//ig;
		$total += length($_);
		$A_num += $_ =~ s/A//ig;
		$T_num += $_ =~ s/T//ig;
		$C_num += $_ =~ s/C//ig;
		$G_num += $_ =~ s/G//ig;
	}
	close INN;
	$/ = "\n";
	my $A_ratio = sprintf "%.4f",$A_num/$total;
	my $T_ratio = sprintf "%.4f",$T_num/$total;
	my $C_ratio = sprintf "%.4f",$C_num/$total;
	my $G_ratio = sprintf "%.4f",$G_num/$total;
	$ref->{"A"} = $A_ratio;
	$ref->{"T"} = $T_ratio;
	$ref->{"C"} = $C_ratio;
	$ref->{"G"} = $G_ratio;
}

sub SNPvalue_Binomial{
	my ($info,$r,$p)=@_;
	my $ratio=1/$p;
	my @I=split /,/,$info;
	my %info_hash=("A",$I[0],"C",$I[1],"G",$I[2],"T",$I[3]);
	my $totalDep=$I[0]+$I[1]+$I[2]+$I[3];
	my $diffDep=0;
	foreach my $base (keys %info_hash){
		if($base ne $r){
			$diffDep+=$info_hash{$base};
		}
	}
	my $lower_tail_p = &pbinom($diffDep, $totalDep, $ratio, "lower");
	return ($lower_tail_p);
}

sub SNPvalue_Frequency{
	my ($s,$r)=@_;
	my @S=split //,$s;
	my $totalDep=@S;
	my $diffDep=0;
	foreach my $base (@S){
		if($base ne $r){
			$diffDep++;
		}
	}
	my $ratio=sprintf "%.4f",$diffDep/$totalDep;
	return ($diffDep,$ratio);
}


sub Fdr{
	my $p=shift;
	my $n=@$p;
	my $temp;
	my $i=0;  my $pre;
	my @Fdr_result;
	foreach my $pvalue (@$p){
		if($i==0){$temp=$n;}elsif($pvalue<$pre){$temp=$n-$i}
		$pre=$pvalue;
		my $cor=$pvalue*$n/$temp;
		push @Fdr_result,$cor;
		$i++;
	}
	return @Fdr_result;
}

#################################################################################################
## Calculate the p for a given number of successes in a Bernoulli experiment.
## Example:
## my $p = dbinom(2, 10, 0.3);
sub dbinom {
	my ($x, $n, $pie) = @_; # $x: number of successes; $n: number of trials; $pie: probability of success on each trial.
		die if $x > $n || $x < 0 || $n < 0;
	die unless $pie > 0 && $pie < 1;
## calculate the factorial part of the binomial formula, use logarithmic operation to save computing time.
	my ($a, $b);
	if ($x > $n/2) {
		$a = $x;
		$b = $n - $x;
	} else {
		$a = $n - $x;
		$b = $x;
	}
	my $log_fac1 = 0;
	for (my $i = $n; $i > $a; $i --) {
		$log_fac1 += log($i);
	}
	my $log_fac2 = 0;
	for (my $i = $b; $i >= 1; $i --) {
		$log_fac2 += log($i);
	}
	my $log_fac = $log_fac1 - $log_fac2;

## calculate the p for a given number of successes.
	my $log_p = $log_fac + $x * log($pie) + ($n - $x) * log(1 - $pie);
	my $p = exp($log_p);
	return $p;
}

sub pbinom {
## $x: number of successes; $n: number of trials; $pie: probability of success on each trial.
## $lower_tail: if TRUE, probabilities are P[X <= x], otherwise, P[X > x].
	my ($x, $n, $pie, $tail) = @_;
	my $tail_p = 0;
	if ($tail eq "lower") {
		if ($x <= $n/2) {
			for (my $i = 0; $i <= $x; $i ++) {
				$tail_p += dbinom($i, $n, $pie);
			}
		} else {
			for (my $i = $x + 1; $i <= $n; $i ++) {
				$tail_p += dbinom($i, $n, $pie);
			}
			$tail_p = 1 - $tail_p;
		}
	} elsif ($tail eq "upper") {
		if ($x <= $n/2) {
			for (my $i = 0; $i <= $x - 1; $i ++) {
				$tail_p += dbinom($i, $n, $pie);
			}
			$tail_p = 1 - $tail_p;
		} else {
			for (my $i = $x; $i <= $n; $i ++) {
				$tail_p += dbinom($i, $n, $pie);
			}
		}
	} else {
		die "tail is not set correctly, set it to lower or upper\n";
	}
	$tail_p = dbinom($x, $n, $pie) if $tail_p <= 0;
	return $tail_p;
}

sub formatP {
	my ($p) = @_;
	if ($p > 0.0001) {
		$p = sprintf "%.8f", $p;
	} else {
		$p = sprintf "%.8e", $p;
	}
	return $p;
}



