#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd 'abs_path';
my ($outdir, $ref, $help, $config, $n, $bwa, $samtools, $k, $l, $t, $index, $run, $junction);
$outdir ||= "./";
$n ||= 0.08;
$k ||= 3;
$t ||= 1;
$l ||= 32;

GetOptions(
		"outDir:s"=>\$outdir,
		"ref:s"=>\$ref,
		"config:s"=>\$config,
		"index:s"=>\$index,
		"n:s"=>\$n,
		"k:s"=>\$k,
		"t:s"=>\$t,
		"l:s"=>\$l,
		"bwa:s"=>\$bwa,
		"samtools:s"=>\$samtools,
		"junction:s"=>\$junction,
		"run"=>\$run,
		"help"=>\$help,
		);
if (!defined $ref || !defined $outdir || !defined $bwa || !defined $config || $help || !defined $samtools) {
	print <<"Usage End."; 

Description:
	Part2 of RES-Scanner alignment pipeline.
	Step3:	Call BWA for read mapping.
	Step4:	Merge bam files for each sample.

Options:
	--ref        FILE   Index of reference, please input the absolute path.
	--outDir     STR    Set the output directory, must be the same with that of part1. default="./"
	--n          NUM    Maximum mismatch number (int) or rate (float) for read alignment, required by BWA. [0.08]
	--k          INT    Maximum mismatch number in the seed sequence for read alignment, required by BWA [3]
	--l          INT    Length of seed sequence for read alignment, required by BWA. [32]
	--t          INT    Number of threads for BWA mapping [1]
	--config     FILE   The configuration file that contains the information of DNA-Seq or RNA-Seq data.
	--bwa        FILE   The absolute path of BWA software pre-installed in local machine.
	--samtools   FILE   The absolute path of SAMtools software pre-installed in local machine.
	--index      NUM    The reference genome is indexed in part1 or not, '1' for yes, '0' for not, default yes. [1]
	--junction   FILE   The file of junction information with POS format. Force --index. [null]
	--run               Run the jobs directly with serial working mode.
	--help              Show the help information.

Usage:
	perl $0 --config config.file --ref reference.fa --outDir ./outdir --bwa /path_to_bwa/bwa --samtools /path_to_samtools/samtools --index 1

Note: 
Before going ahead to the part2, please make sure that jobs in part1 are all finished successfully. After executing command line of part2, two shell scripts named 'step2.sh' and 'step3.sh' will be generated in the directory './outdir/'. Then run the scripts with command 'sh ./outdir/step2.sh' and 'sh ./outdir/step3.sh' orderly or submit the scripts as jobs to computing servers. If one step contains multiple shell scripts, users can run these scripts in parallel. Before going ahead to the next step, please make sure that all the current jobs are completed successfully.
	
Usage End.

		exit;
}
#################################################################################################################

die "Error: $outdir is not existent" unless -e $outdir;
$outdir=abs_path $outdir;
$ref=abs_path $ref;

if($index){
	my $refname=basename $ref;
	$ref="$outdir/index/$refname";
	die "Warning: work of part1 is not completed! as $ref is not existent.\n" unless -e $ref;
}else{
	if(! -e "$ref.bwt" || ! -e "$ref.rbwt"){
		die "Warnning: $ref is not indexed, please call pipeline with '--index 1'.\n";
	}
}

my $convert_BAM_coordinate4junction="$Bin/bin/convert_BAM_coordinate4junction.pl";
foreach my $script ($bwa, $samtools, $convert_BAM_coordinate4junction) {
	die "$script is not existent!" unless -e $script;
}

############################ call bwa mapping  ########################
if(-e "$outdir/step2.sh"){
	system "rm $outdir/step2.sh";
}
my %bam_hash;
open IN,"$config" or die $!;
while(<IN>){
	chomp;
	my @A=split /\s+/;
	my ($tag,$lib,$lane,$insertSize)=@A[0..3];
	die "Error: $outdir/$tag/$lib/$lane is not existent!\n" unless -d "$outdir/$tag/$lib/$lane";
	my $j;
	for ($j=1;-d "$outdir/$tag/$lib/$lane/fq$j/";$j++){}
	my $file_num=$j-1;
	my $x= $A[3] * 1.3;
	if(@A==6){
		for (my $i=1;$i<=$file_num;$i++){
			my $a_fq="$outdir/$tag/$lib/$lane/fq$i/$i.1.fq.gz";
			my $b_fq="$outdir/$tag/$lib/$lane/fq$i/$i.2.fq.gz";
			die "Error: $a_fq is not existent!\n" unless -e $a_fq;
			die "Error: $b_fq is not existent!\n" unless -e $b_fq;
			my $Qual64=checkQuality($a_fq);
			open OUT,">$outdir/$tag/$lib/$lane/fq$i/run.sh";
			print OUT "if [ ! -f \"$outdir/$tag/$lib/$lane/step1.log\" ];then echo \"Warning: $outdir/step1.sh work is not completed! as $outdir/$tag/$lib/$lane/step1.log is not existent\"\nexit 0\nfi\n";
			if($Qual64){
				print OUT "$bwa aln -I -k $k -l $l -n $n -t $t $ref $a_fq > $outdir/$tag/$lib/$lane/fq$i/$lane.1.sai 2>$outdir/$tag/$lib/$lane/fq$i/$lane.1.sai.log\n";
				print OUT "$bwa aln -I -k $k -l $l -n $n -t $t $ref $b_fq > $outdir/$tag/$lib/$lane/fq$i/$lane.2.sai 2>$outdir/$tag/$lib/$lane/fq$i/$lane.2.sai.log\n";	
			}else{
				print OUT "$bwa aln -k $k -l $l -n $n -t $t $ref $a_fq > $outdir/$tag/$lib/$lane/fq$i/$lane.1.sai 2>$outdir/$tag/$lib/$lane/fq$i/$lane.1.sai.log\n";
				print OUT "$bwa aln -k $k -l $l -n $n -t $t $ref $b_fq > $outdir/$tag/$lib/$lane/fq$i/$lane.2.sai 2>$outdir/$tag/$lib/$lane/fq$i/$lane.2.sai.log\n";
			}
			print OUT "$bwa sampe -a $x $ref $outdir/$tag/$lib/$lane/fq$i/$lane.1.sai $outdir/$tag/$lib/$lane/fq$i/$lane.2.sai $a_fq $b_fq 2>$outdir/$tag/$lib/$lane/fq$i/$lane.bwa.bam.log | $samtools view -hbS - > $outdir/$tag/$lib/$lane/fq$i/$lane.bwa.bam 2>>$outdir/$tag/$lib/$lane/fq$i/$lane.bwa.bam.log\n";
			print OUT "rm $outdir/$tag/$lib/$lane/fq$i/$lane.1.sai $outdir/$tag/$lib/$lane/fq$i/$lane.2.sai\n";
			print OUT "echo $tag $lane BWA mapping work-complete! > $outdir/$tag/$lib/$lane/fq$i/$lane.log\n";
			print OUT "echo $tag $lane BWA mapping work-complete!\n";
			close OUT;
			system("echo sh $outdir/$tag/$lib/$lane/fq$i/run.sh >> $outdir/step2.sh");
			push @{$bam_hash{$tag}{$lib}},"$outdir/$tag/$lib/$lane/fq$i/$lane.bwa.bam";
		}
	}elsif(@A==5){
		for (my $i=1;$i<=$file_num;$i++){
			my $a_fq="$outdir/$tag/$lib/$lane/fq$i/$i.1.fq.gz";
			die "Error: $a_fq is not existent!\n" unless -e $a_fq;
			my $Qual64=checkQuality($a_fq);
			open OUT,">$outdir/$tag/$lib/$lane/fq$i/run.sh";
			print OUT "if [ ! -f \"$outdir/$tag/$lib/$lane/step1.log\" ];then echo \"Warning: $outdir/step1.sh work is not completed! as $outdir/$tag/$lib/$lane/step1.log is not existent\"\nexit 0\nfi\n";
			if($Qual64){
				print OUT "$bwa aln -I -k $k -l $l -n $n -t $t $ref $a_fq > $outdir/$tag/$lib/$lane/fq$i/$lane.1.sai 2>$outdir/$tag/$lib/$lane/fq$i/$lane.1.sai.log\n";
			}else{
				print OUT "$bwa aln -k $k -l $l -n $n -t $t $ref $a_fq > $outdir/$tag/$lib/$lane/fq$i/$lane.1.sai 2>$outdir/$tag/$lib/$lane/fq$i/$lane.1.sai.log\n";
			}
			print OUT "$bwa samse $ref $outdir/$tag/$lib/$lane/fq$i/$lane.1.sai $a_fq 2>$outdir/$tag/$lib/$lane/fq$i/$lane.bwa.bam.log | $samtools view -hbS - > $outdir/$tag/$lib/$lane/fq$i/$lane.bwa.bam 2>>$outdir/$tag/$lib/$lane/fq$i/$lane.bwa.bam.log\n";
			print OUT "rm $outdir/$tag/$lib/$lane/fq$i/$lane.1.sai\n";
			print OUT "echo $tag $lane BWA mapping work-complete! > $outdir/$tag/$lib/$lane/fq$i/$lane.log\n";
			print OUT "echo $tag $lane BWA mapping work-complete!\n";
			close OUT;
			system("echo sh $outdir/$tag/$lib/$lane/fq$i/run.sh >> $outdir/step2.sh");
			push @{$bam_hash{$tag}{$lib}},"$outdir/$tag/$lib/$lane/fq$i/$lane.bwa.bam";
		}
	}else{
		die "Error: Unrecognized config file: $config !\nConfigure file must be 5 columns or 6 columns\n\n";
	}
}
close IN;

######################################   merge bam files for each samples  ######################################

my %config_hash;
if(-e "$outdir/step3.sh"){
	system "rm $outdir/step3.sh";
}
foreach my $tag (keys %bam_hash){
	foreach my $lib (keys %{$bam_hash{$tag}}){
	my $dir="$outdir/$tag/$lib";
	die "Error: $dir is not existent" unless -e $dir;
	my @With= @{$bam_hash{$tag}{$lib}};
	open OUT,">$dir/mergeBam.$lib.sh" or die $!;
	if(@With>1){
		foreach my $file(@With){
			my $logfile=$file;
			$logfile=~s/bwa.bam$/log/;
			print OUT "if [ ! -f \"$logfile\" ];then echo \"Warning: work is not completed! as $logfile is not existent\"\nexit 0\nfi\n";	
		}
		print OUT "$samtools merge $dir/$lib.merge.bam";
		foreach my $file(@With){
			print OUT " $file";
		}
		print OUT "\n";
	}elsif(@With==1){
		my $file=$With[0];
		my $logfile=$file;
		$logfile=~s/bwa.bam$/log/;
		print OUT "if [ ! -f \"$logfile\" ];then echo \"Warning: work is not completed! as $logfile is not existent\"\nexit 0\nfi\n";
		print OUT "ln -sf $With[0] $dir/$lib.merge.bam\n";
	}else{
		die "Error: no bam files to be processed";
	}
	if(defined $junction){
		print OUT "if [ ! -f \"$outdir/index/junctionFlankSequenceRegion.txt\" ];then echo \"Warning: $outdir/index/junctionFlankSequenceRegion.txt is not existent when --junction is activated\"\nexit 0\nfi\n";
		print OUT "perl $convert_BAM_coordinate4junction $outdir/index/junctionFlankSequenceRegion.txt $dir/$lib.merge.bam $dir/$lib.merge.converted.bam $samtools\n";
		print OUT "$samtools sort $dir/$lib.merge.converted.bam $dir/$lib.merge.converted.sort\n";
		print OUT "mv -f $dir/$lib.merge.bam $dir/$lib.merge.unconverted.bam\n";
		print OUT "rm -f $dir/$lib.merge.converted.bam\n";
		print OUT "mv -f $dir/$lib.merge.converted.sort.bam $dir/$lib.merge.bam\n";
	}
	print OUT "echo Merging BAM files of $tag $lib is completed ! > $dir/mergeBam.$lib.log\n";
	print OUT "echo '$dir/$lib.merge.bam' is the final mapping result of $tag $lib. > $dir/mergeBam.$lib.log\n";
	print OUT "echo Merging BAM files of $tag $lib is completed !\n";
	print OUT "echo '$dir/$lib.merge.bam' is the final mapping result of $tag $lib.\n";
	close OUT;
	$config_hash{$lib}{$tag}="$dir/$lib.merge.bam";
	system "echo sh $dir/mergeBam.$lib.sh >> $outdir/step3.sh";
}
}

## generate configuration file automatically for RES_Scanner_indentification pipeline

open OUT,">$outdir/RES_Scanner_indentification_config.txt" or die $!;
foreach my $sample (sort keys %config_hash){
	my $flag=1;
	if(!exists $config_hash{$sample}{DNA}){
		print STDERR "Warning: DNA-Seq data for $sample is not existent!\n";
		$flag=0;
	}
	if(!exists $config_hash{$sample}{RNA}){
		print STDERR "Warning: RNA-Seq data for $sample is not existent!\n";
		$flag=0;
	}
	if($flag){
		print OUT "$sample\t$config_hash{$sample}{DNA}\t$config_hash{$sample}{RNA}\n";
	}
}
close OUT;

###  Usage notice  ###
if($run){
	system "sh $outdir/step2.sh";
	system "sh $outdir/step3.sh";
}else{
	print STDERR "*****************************************************************************************************************************\n";
	print STDERR "There are two shell scripts to be finished, users have to run the jobs step by step, before going ahead to the next step, please make sure the current job is finished successfully. If one step contains multiple shell scripts, users can run these scripts in parallel.\n";
	print STDERR "Please run step2 with command 'sh $outdir/step2.sh' or by submitting it to computing servers.\n";
	print STDERR "Please run step3 with command 'sh $outdir/step3.sh' or by submitting it to computing servers.\n";
	print STDERR "After all the jobs are completed, the file named 'RES_Scanner_indentification_config.txt' in the output directory can be the configuration file for RES_Scanner_indentification pipeline\n";
	print STDERR "*****************************************************************************************************************************\n";
}
################ subrouting ###############
sub checkQuality {
	my $file=shift;
	if($file=~/\.gz$/){
		open INN,"gunzip -c $file |" or die $!;
	}else{
		open INN,"$file" or die $!;
	}
	my $count=0;
	my %hash;
	while(<INN>){
		my $id=$_;
		my $seq=<INN>;
		my $tag=<INN>;
		my $qual=<INN>;
		chomp ($id,$seq,$tag,$qual);
		$count++;
		last if $count>10000;
		my @Array=split //,$qual;
		foreach my $q (@Array){
			my $value=ord($q);
			$hash{$value}=1;
		}
	}
	close INN;
	my @array=(sort {$a<=>$b} keys %hash);
	my $value=$array[0];
	if($value<64){
		return 0;
	}else{
		return 1;
	}
}

