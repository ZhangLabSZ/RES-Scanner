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
die "Usage:perl $0 <RNAedtingType>\n" unless @ARGV == 1;

my %count;
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open (IN,$ARGV[0]) or die $!;
}
while (<IN>) {
	next if $_=~/^#|^1.Chromosome_ID|^Chromosome/;
	chomp;
	my @info = split /\t/;
	my $id = join "#",@info;
	$count{$id} ++;
	die "Something($info[0]\t$id) was wrong!\n" if ( $count{$id} >= 2);
	my $strand = $info[2];
	my $line = join "\t",$id,$info[0],$strand,$info[1],$info[1];
	print "$line\n";
}
close IN;


