#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my ($outdir, $ref, $help, $index, $config, $split, $bwa);
$outdir ||= "./";
$index ||= 1;

GetOptions(
		"outDir:s"=>\$outdir,
		"ref:s"=>\$ref,
		"index:s"=>\$index,
		"config:s"=>\$config,
		"split"=>\$split,
		"bwa:s"=>\$bwa,
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
	--bwa       FILE    The absolute path of BWA software pre-installed in local machine.
	--config    FILE    The configuration file that contains the information of DNA-Seq and RNA-Seq data.
	                    Format of configuration file (tab-delimited table): 
	                       1) Tag (DNA/RNA)
	                       2) sample ID
	                       3) lane ID
	                       4) Insert size
	                       5) 1.fastq file
	                       6) 2.fastq file (Note: the sixth column can be absent for single-end data).
	--split            Split fastq files for multithreading. [null]
	--help             Show the help information.

Usage:
	perl $0 --outDir ./outdir/ --ref reference.fa --bwa /path_to_bwa/bwa --index 1 --config config.file --split

Note:
After executing the command line of part1, one (--index 0) or two (--index 1) shell scripts will be generated in the directory './outdir/'. Then run the scripts with command 'sh ./outdir/step0.sh' and 'sh ./outdir/step1.sh' orderly or submit the scripts as jobs to computing servers. In addition, if one step contains multiple shell scripts, users can run these scripts in parallel. Before going ahead to the next step, please make sure that all the current jobs are completed successfully.

Usage End.

		exit;
}
#################################################################################################################
my $Split_Fq="$Bin/Split_Fq.pl";
die "$Split_Fq not exist!\n" unless -e $Split_Fq;

mkdir $outdir unless -e $outdir;
$outdir=abs_path $outdir;
$ref=abs_path $ref;

foreach my $script ($bwa) {
	die "$script not exist!" unless -e $script;
}


## index reference
if($index){
	my $index_outdir= "$outdir/index";
	mkdir $index_outdir unless -e $index_outdir;

	my $refLength=0;
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
	my $refname=basename $ref;
	open OUT,">$outdir/step0.sh" or die $!;
	print OUT "ln -s $ref $index_outdir/$refname\n";
	print OUT "$bwa index -a $a_parameter $index_outdir/$refname 2>$index_outdir/$refname.log\n";
	print OUT "echo 'indexing of reference (step0) is completed!' > $outdir/step0.log\n";
	print OUT "echo 'indexing of reference (step0) is completed!'\n";
	close OUT;

}else{
	if(! -e "$ref.bwt" || ! -e "$ref.rbwt"){
		die "Warnning: $ref is not indexed, please call pipeline with '--index 1'.\n";
	}
}

########  split fastq files to limit memory  ###########
my %checkLane;
open IN,"$config" or die $!;
while(<IN>){
	chomp;
	my @A=split /\s+/,$_;
	if($A[0] ne "RNA" && $A[0] ne "DNA"){
		die "Warning: tag '$A[0]' should be 'DNA' or 'RNA' in configure file $config\n$_\n";
	}
	$checkLane{$A[2]}++;
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

