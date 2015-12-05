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
use Cwd 'abs_path';
use File::Basename qw(basename);
use Getopt::Long;
die "<bam> <outdir> <samtools> --ss [1|0] --mq [20] --rmdup [1|0] --DNA|--RNA  --uniqTag [0|1]\n" unless @ARGV>=4;
my ($ss,$mq,$rmdup,$DNA,$RNA,$uniqTag);
$ss ||= 1;
$mq ||= 20;
$uniqTag ||= 1;
$rmdup ||= 1;
my $samtools=$ARGV[2];
die "$samtools does not exists!\n" unless -e $samtools;

GetOptions(
		"ss:s"=>\$ss,
		"mq:s"=>\$mq,
		"rmdup:s"=>\$rmdup,
		"DNA"=>\$DNA,
		"RNA"=>\$RNA,
		"uniqTag:s"=>\$uniqTag,
		);

if($ARGV[0]=~/\.sam.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}elsif($ARGV[0]=~/\.sam$/){
	open IN,"$ARGV[0]" or die $!;
}elsif($ARGV[0]=~/\.bam$/){
	open IN,"$samtools view -h $ARGV[0] |" or die $!;
}else{
	die "the suffix of input file $ARGV[0] should be .bam or .sam.gz/.sam!";
}
my $outdir = abs_path($ARGV[1]);
my $filename= basename($ARGV[0]);
if($DNA){
	$ss=0;
}

if($ss){
	open OUT1,"| $samtools view -hbS - >$outdir/$filename.negative.bam" or die $!;
	open OUT2,"| $samtools view -hbS - >$outdir/$filename.positive.bam" or die $!;
}else{
	open OUT3,"| $samtools view -hbS - >$outdir/$filename.best.bam" or die $!;
}

while(<IN>){
	chomp;
	my @array4tile=split /\t/;
	if(@array4tile<9){
		if($ss){
			print OUT1 "$_\n";
			print OUT2 "$_\n";
		}else{
			print OUT3 "$_\n";
		}
		next;
	}

	my @A=split /\t/;
	if($A[4]<$mq){next;}
	my $flag=$A[1];
	my %hash;
	&index($flag,\%hash);
	if(exists $hash{1024} && $rmdup){next;}
	if($RNA){
        if( $uniqTag && ($_!~/XT:A:U\s+/ || $_!~/X0:i:1\s+/ || $_!~/X1:i:0\s+/) ){
			next;
		}
		if(exists $hash{256}){next;}
		if($flag eq 67 || $flag eq 131 || $flag eq 115 || $flag eq 179){next;}
		if(exists $hash{32} && exists $hash{16}){next;}
		if($ss){
			if(  exists $hash{64}  && exists $hash{32} && $A[8]>=0  ){
				print OUT1 "$_\n";
			}elsif( exists $hash{128}  && exists $hash{16} && $A[8]<=0 ){
				print OUT1 "$_\n";
			}elsif( exists $hash{64}  && exists $hash{16} && $A[8]<=0 ){
				print OUT2 "$_\n";
			}elsif( exists $hash{128}  && exists $hash{32} && $A[8]>=0){
				print OUT2 "$_\n";
			}elsif(exists $hash{64} && !exists $hash{16} && !exists $hash{32}){
				print OUT1 "$_\n";
			}elsif( exists $hash{128} && !exists $hash{16} && !exists $hash{32}){
				print OUT2 "$_\n";
			}elsif ($flag==16){
				print OUT2 "$_\n";
			}elsif ($flag==0){
				print OUT1 "$_\n";
			}
		}else{
			print OUT3 "$_\n";
		}
	}elsif($DNA){
		print OUT3 "$_\n";
	}else{
		die "Error: option --DNA or --RNA is not activated!\n";
	}
}
close IN;
if($ss){
	close OUT1;
	close OUT2;
}else{
	close OUT3;
}



sub index{
	my $num=shift;
	my $hash=shift;
	my $temp=1;
	if($num>=1){
		while($temp<=$num){
			$temp=$temp*2;
		}
		$temp=$temp/2;
		$$hash{$temp}=1;
		my $diff=$num-$temp;
		if($diff>0){
			&index($diff,\%$hash);
		}
	}
}

