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
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my ($outdir, $ref, $help, $index, $config, $split, $bwa, $junction, $run, $readlen);
$outdir ||= "./";
$index ||= 1;

GetOptions(
		"outDir:s"=>\$outdir,
		"ref:s"=>\$ref,
		"junction:s"=>\$junction,
		"readlen:s"=>\$readlen,
		"index:s"=>\$index,
		"config:s"=>\$config,
		"split"=>\$split,
		"bwa:s"=>\$bwa,
		"run"=>\$run,
		"help"=>\$help,
		);
if (!defined $ref || !defined $outdir || !defined $bwa || !defined $config || $help) {
	print <<"Usage End."; 

Description:
	Part1 of RES-Scanner alignment pipeline. 
	Step1: Index the reference genome for BWA alignment (optional).
	Step2: Preprocess the fastq files before BWA alignment. As the memory requirement is positively correlated with the size of fastq files, to limit the memory usage and to speed up the alignment step, the fastq file can be split into several files with no more than 3 Gb per file (optional).

options:
	--ref       FILE    The reference genome in FASTA format.
	--outDir    STR     The output directory, default="./"
	--index	    NUM     Index the reference genome for BWA or not, '1' for yes, '0' for no, default yes. [1]
	--junction  FILE    The file of junction information with POS format. Force --index. (Note: The length of all RNA-seq reads should keep the same.) [null]
	--readlen   INT     The length (bp) of RNA-seq reads for --junction option.[null]
	--bwa       FILE    The absolute path of BWA software pre-installed in local machine.
	--config    FILE    The configuration file that contains the information of DNA-Seq and RNA-Seq data.
	                    Format of configuration file (tab-delimited table): 
	                       1) Tag (DNA/RNA)
	                       2) sample ID
	                       3) lane ID
	                       4) Insert size
	                       5) 1.fastq file
	                       6) 2.fastq file (Note: the sixth column can be absent for single-end data).
	--split             Split fastq files for multithreading. [null]
	--run               Run the jobs directly with serial working mode.
	--help              Show the help information.

Usage:
	perl $0 --outDir ./outdir/ --ref reference.fa --bwa /path_to_bwa/bwa --index 1 --config config.file --split

Note:
After executing the command line of part1, one (--index 0) or two (--index 1) shell scripts will be generated in the directory './outdir/'. Then run the scripts with command 'sh ./outdir/step0.sh' and 'sh ./outdir/step1.sh' orderly or submit the scripts as jobs to computing servers. In addition, if one step contains multiple shell scripts, users can run these scripts in parallel. Before going ahead to the next step, please make sure that all the current jobs are completed successfully.

Usage End.

		exit;
}
#################################################################################################################
my $Split_Fq="$Bin/bin/Split_Fq.pl";
my $junction2flankSequenceRegion="$Bin/bin/junction2flankSequenceRegion.pl";
my $flankSequenceRegion2fa="$Bin/bin/flankSequenceRegion2fa.pl";
die "$Split_Fq is not existent!\n" unless -e $Split_Fq;
die "$junction2flankSequenceRegion is not existent!\n" unless -e $junction2flankSequenceRegion;

mkdir $outdir unless -e $outdir;
$outdir=abs_path $outdir;
$ref=abs_path $ref;

foreach my $script ($bwa) {
	die "$script not exist!" unless -e $script;
}

###check configure file

my %checkLane;
my %pair;
open IN,"$config" or die $!;
while(<IN>){
	chomp;
	my @A=split /\s+/,$_;
	if($A[0] ne "RNA" && $A[0] ne "DNA"){
		die "Warning: tag '$A[0]' should be 'DNA' or 'RNA' in configure file $config\n$_\n";
	}
	$checkLane{$A[2]}++;
	$pair{$A[1]}{$A[0]}=1;
	for(my $i=4;$i<@A;$i++){
		die "Error: $A[$i] is not existent in $config !\n" unless -e $A[$i];
	}
}
close IN;

foreach my $lane (keys %checkLane){
	if($checkLane{$lane}>1){
		die "Error: lane ID: $lane should be unique in file $config!\n";
	}
}

my $pair_flag=0;
foreach my $sample (sort keys %pair){
	if(!exists $pair{$sample}{'DNA'}){
		print STDERR "Warning: There is no matched DNA-Seq data for sample '$sample' in configure file $config!\n";
		$pair_flag=1;
	}
	if(!exists $pair{$sample}{'RNA'}){
		print STDERR "Warning: There is no matched RNA-Seq data for sample '$sample' in configure file $config!\n";
		$pair_flag=1;
	}
}
#if($pair_flag){
#	die "Error: The configure file $config should be revised with matched DNA-Seq and RNA-Seq data for each sample name.\n";
#}


## index reference
if($index){
	my $index_outdir= "$outdir/index";
	mkdir $index_outdir unless -e $index_outdir;
	my $refname=basename $ref;
	my $junctionCount=0;
	my $readLength=0;
	open OUT,">$outdir/step0.sh" or die $!;
	## add junction sequence
	if($junction){
		die "Error: $junction is not existent!\n" unless -e $junction;
		$junction=abs_path $junction;
		
		open JUN,"$junction" or die $!;
		while(<JUN>){
			$junctionCount++;
		}
		close JUN;

		if(defined $readlen){
			$readLength=$readlen;
		}else{
		my $readCount=0;
		my $totalLength=0;
		open IN,"$config" or die $!;
		while(<IN>){
			chomp;
			my @A=split /\s+/,$_;
			if($A[0] eq "RNA"){
				if($A[4]=~/\.gz$/){
					open INN,"gunzip -c $A[4] |" or die $!;
				}else{
					open INN,"$A[4]" or die $!;
				}
				while(<INN>){
					my $r=<INN>;
					chomp $r;
					$totalLength+=length($r);
					<INN>;
					<INN>;
					$readCount++;
					last if $readCount>=100;
				}
				close INN;
				last;
			}
		}	
		close IN;
		$readLength=$totalLength/$readCount;
		}
		print OUT "perl $junction2flankSequenceRegion $ref $junction $readLength > $index_outdir/junctionFlankSequenceRegion.txt\n";
		print OUT "perl $flankSequenceRegion2fa $ref $index_outdir/junctionFlankSequenceRegion.txt > $index_outdir/junctionFlankSequenceRegion.fa\n";
		print OUT "cat $ref $index_outdir/junctionFlankSequenceRegion.fa > $index_outdir/$refname\n";
	}else{
		print OUT "ln -s $ref $index_outdir/$refname\n";
	}

	##
	my $refLength=$junctionCount*($readLength-1)*2;
	if($ref=~/\.gz$/){
		open IN,"gunzip -c $ref |" or die $!;
	}else{
		open IN,"$ref" or die $!;
	}
	while(<IN>){
		chomp;
		next if $_=~/^>/;
		$refLength+=length($_);
	}
	close IN;
	my $a_parameter;
	if($refLength>=2000000000){
		$a_parameter="bwtsw";
	}else{
		$a_parameter="is";
	}
	print OUT "$bwa index -a $a_parameter $index_outdir/$refname 2>$index_outdir/$refname.log\n";
	print OUT "echo 'indexing of reference (step0) is completed!' > $outdir/step0.log\n";
	print OUT "echo 'indexing of reference (step0) is completed!'\n";
	if($junction){
		print OUT "echo '$index_outdir/$refname' is the file of assembly plus exonic sequences surrounding all provided splicing junctions\n";
	}
	close OUT;

}else{
	if(! -e "$ref.bwt" || ! -e "$ref.rbwt"){
		die "Warnning: $ref is not indexed, please call pipeline with '--index 1'.\n";
	}
}

########  split fastq files to limit memory  ###########
open IN,"$config" or die $!;
while(<IN>){
	chomp;
	my @A=split /\s+/,$_;
	my ($tag,$sample,$lane)=($A[0],$A[1],$A[2]);
	mkdir "$outdir/$tag" unless -e "$outdir/$tag";
	mkdir "$outdir/$tag/$sample" unless -e "$outdir/$tag/$sample";
	mkdir "$outdir/$tag/$sample/$lane" unless -e "$outdir/$tag/$sample/$lane";
	open OUT,">$outdir/$tag/$sample/$lane/prepareFq.sh" or die $!;
	if(@A==6){
		my $fq1= abs_path $A[4];
		my $fq2= abs_path $A[5];
		if($split){
			print OUT "perl $Split_Fq --fq $fq1 --fq $fq2 --outdir $outdir/$tag/$sample/$lane\n";
			print OUT "echo 'Preprocessing of the fastq files of $lane ($tag step1) is completed!' >> $outdir/$tag/$sample/$lane/step1.log\n";
			print OUT "echo 'Preprocessing of the fastq files of $lane ($tag step1) is completed!'\n";
		}else{
			print OUT "mkdir -p $outdir/$tag/$sample/$lane/fq1/\n";
			if($fq1=~/\.gz$/){
				print OUT "ln -s $fq1 $outdir/$tag/$sample/$lane/fq1/1.1.fq.gz\n";
			}else{
				print OUT "gzip -c $fq1 > $outdir/$tag/$sample/$lane/fq1/1.1.fq.gz\n";
			}
			if($fq2=~/\.gz$/){
				print OUT "ln -s $fq2 $outdir/$tag/$sample/$lane/fq1/1.2.fq.gz\n";
			}else{
				print OUT "gzip -c $fq2 > $outdir/$tag/$sample/$lane/fq1/1.2.fq.gz\n";
			}
			print OUT "echo 'Preprocessing of the fastq files of $lane ($tag step1) is completed!' >> $outdir/$tag/$sample/$lane/step1.log\n";
			print OUT "echo 'Preprocessing of the fastq files of $lane ($tag step1) is completed!'\n";
		}
	}elsif(@A==5){
		my ($fq1)=(abs_path $A[4]);	
		if($split){
			print OUT "perl $Split_Fq --fq $fq1 --outdir $outdir/$tag/$sample/$lane\n";
			print OUT "echo 'Preprocessing of the fastq files of $lane ($tag step1) is completed!' >> $outdir/$tag/$sample/$lane/step1.log\n";
			print OUT "echo 'Preprocessing of the fastq files of $lane ($tag step1) is completed!'";
		}else{
			print OUT "mkdir -p $outdir/$tag/$sample/$lane/fq1/\n";
			if($fq1=~/\.gz$/){
				print OUT "ln -s $fq1 $outdir/$tag/$sample/$lane/fq1/1.1.fq.gz\n";
			}else{
				print OUT "gzip -c $fq1 > $outdir/$tag/$sample/$lane/fq1/1.1.fq.gz\n";
			}
			print OUT "echo 'Preprocessing of the fastq files of $lane ($tag step1) is completed!' >> $outdir/$tag/$sample/$lane/step1.log\n";
			print OUT "echo 'Preprocessing of the fastq files of $lane ($tag step1) is completed!'\n";
		}
	}else{
		die "Error: Unrecognized config file: $config !\nConfigure file must be 5 columns or 6 columns\n";
	}
	close OUT;
	system "echo sh $outdir/$tag/$sample/$lane/prepareFq.sh >> $outdir/step1.sh";
}
close IN;

############################ help information  ########################
if($run){
	if($index){
		system "sh $outdir/step0.sh";
	}
	system "sh $outdir/step1.sh";
}else{
	if($index){
		print STDERR "*****************************************************************************************************************************\n";
		print STDERR "There are two shell scripts to be finished, users have to run the jobs step by step, before going ahead to the next step, please make sure the current job is finished successfully. If one step contains multiple shell scripts, users can run these scripts in parallel.\n";
		print STDERR "Please run step0 with command 'sh $outdir/step0.sh' or by submitting it to computing servers.\n";
		print STDERR "Please run step1 with command 'sh $outdir/step1.sh' or by submitting it to computing servers.\n";
		print STDERR "*****************************************************************************************************************************\n";
	}else{
		print STDERR "*****************************************************************************************************************************\n";
		print STDERR "There are one shell script to be finished, before going ahead to the next part, please make sure the current job is finished successfully. If the step contains multiple shell scripts, users can run these scripts in parallel.\n";
		print STDERR "Please run step1 with command 'sh $outdir/step1.sh' or by submitting it to computing servers.\n";
		print STDERR "*****************************************************************************************************************************\n";
	}
}
