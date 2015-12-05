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
die "perl $0 <edit.site> <knownSNPs.gff>\n" unless @ARGV==2;

my %hash;
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,"$ARGV[0]" or die $!;
}
while(<IN>){
	chomp;
	next if $_=~/^1.Chromosome_ID/;
	my @A=split /\t/;
	$hash{$A[0]}{$A[1]}=1;
}
close IN;

if($ARGV[1]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[1] |" or die $!;
}else{
	open IN,"$ARGV[1]" or die $!;
}
while(<IN>){
	chomp;
	next if $_=~/^#/;
	my @A=split /\t/;
	if(exists $hash{$A[0]}{$A[4]}){
		delete $hash{$A[0]}{$A[4]};
	}
}
close IN;

if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,"$ARGV[0]" or die $!;
}
while(<IN>){
	chomp;
	if($_=~/^1.Chromosome_ID/){
		print "$_\n";
		next;
	}
	my @A=split /\t/;
	if(exists $hash{$A[0]}{$A[1]}){
		print "$_\n";
	}
}
close IN;



