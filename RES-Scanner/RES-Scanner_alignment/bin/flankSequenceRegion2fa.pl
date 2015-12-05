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
die "perl $0 <genome.fa> <cdsPOS2junction.coordinate>\n" unless @ARGV==2;

my %genome;
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,"$ARGV[0]" or die $!;
}
$/=">";
<IN>;
while(<IN>){
	chomp;
	$_=~/(.+?)\n/;
	my $id=(split /\s+/,$1)[0];
	$_=~s/.+?\n//;
	$_=~s/\s+//g;
	$genome{$id}=uc $_;
}
$/="\n";
close IN;


open IN,"$ARGV[1]" or die $!;

while(<IN>){
	chomp;
	next if $_=~/^#/;
	my @A=split /\t/;
	my @B=split /,/,$A[2];
	my $string;
	foreach my $region (@B){
		my @R=split /\-/,$region;
		my $len=$R[1]-$R[0]+1;
		my $str=substr($genome{$A[1]},$R[0]-1,$len);
		$string.=$str;
	}
	print ">$A[0]\t$A[1]\t$A[2]\n$string\n";
}
close IN;
