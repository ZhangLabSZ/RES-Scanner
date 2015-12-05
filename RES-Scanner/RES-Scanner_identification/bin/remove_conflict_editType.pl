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
die "usage: perl $0 <RNAedit.merged.table>\n" unless @ARGV==1;
my @array;
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,$ARGV[0] or die $!;
}
while(<IN>){
	chomp;
	if($_=~/^Chromosome|^#/){
		print "$_\n";
		next;
	}
	my @A=split /\t/;
	push @array,[$A[0],$A[1],$_];
}
close IN;

@array=sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @array;

my ($chr_last,$pos_last,$info_last);
foreach my $value (@array){
	my ($chr,$pos,$info)=@$value;
	if(!defined $chr_last){
		($chr_last,$pos_last,$info_last)=($chr,$pos,$info);
	}else{
		if($chr eq $chr_last && $pos eq $pos_last){
			($chr_last,$pos_last,$info_last)=();
		}else{
			print "$info_last\n";
			($chr_last,$pos_last,$info_last)=($chr,$pos,$info);
		}
	}
}

if(defined $info_last){
	print "$info_last\n";
}
