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
die "<genome.fa> <sam2base>\n" unless @ARGV==2;
my %chrHash;
my %hash;
my $covered_base=0;
my $total_depth=0;
if($ARGV[1]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[1] |" or die $!;
}else{
	open IN,"$ARGV[1]" or die $!;
}
while(<IN>){
	chomp;
	my @A=split /\t/;
	$chrHash{$A[0]}=1;
	if($A[5]>0){
		$covered_base++;
	}
	$hash{$A[5]}++;
	$total_depth+=$A[5];
}
close IN;


my $len;
if ($ARGV[0] =~ /\.gz$/) {
	open IN, "gunzip -c $ARGV[0] |" or die $!;;
} else {
	open IN, "$ARGV[0]" or die $!;;
}
$/ = ">";
<IN>;
while (<IN>) {
	chomp;
	$_=~/(.+?)\n/;
	my $id=(split /\s+/,$1)[0];
	s/.+\n//;
	s/\s+//g;
	if(exists $chrHash{$id}){
		$len += length($_);
	}
}
$/="\n";
close IN;

my @array= sort {$hash{$b}<=>$hash{$a}} keys %hash;
my $peak_depth;
if(@array){
	$peak_depth=$array[0];
}else{
	$peak_depth=0;
}
my $average_depth= $covered_base==0 ? 0 : sprintf "%.2f",$total_depth/$covered_base;
my $coverage=sprintf "%.2f",$covered_base/$len*100;
print "#total_bases\t#covered_bases\tcoveage(%)\tcoverage(X)\tpeak_coverage(X)\n";
print "$len\t$covered_base\t$coverage\t$average_depth\t$peak_depth\n";






