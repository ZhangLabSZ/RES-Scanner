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
die "Usage: perl $0 <blat_output.psl> <RES_final_result.txt>\n" unless @ARGV >= 2;
my $alignRate ||= 0.5;
GetOptions(
		"alignRate:s"=>\$alignRate,
		);

my (%hash,%paralogous);
if ($ARGV[0] =~ /\S+\.gz$/) {
	open (IN,"gunzip -c $ARGV[0] |") or die $!;
} else {
	open (IN,$ARGV[0]) or die $!;
}
for (my $i=0;$i<5;$i++) {
	<IN>;
}
while (<IN>) {
	my @info = split /\s+/;
	my @query = split /;/,$info[9];
	if(@query!=4){die "Error: chromosome/scaffold ID contatins semicolon ';' !\n";}
	my ($scaffold,$site,$st,$ed) = @query;
	if($info[0]<($ed-$st+1)*$alignRate || $info[5]>=5 || $info[7]>=5 || $info[1]/$info[0]>0.05){   #align rate>=0.5 & gapSize<5bp & -minIdentity=95
		next;
	}
	my ($match_st,$match_ed,$site_offset)=($info[11]+1,$info[12],$site-$st+1);
	if($site_offset<$match_st || $site_offset>$match_ed){   # edit site should be located in the matching regions.
		next;
	}
	if(exists $hash{$info[9]}){
		$paralogous{$scaffold}{$site}=1;
	}else{
		$hash{$info[9]}=1;
	}
}
close IN;



open IN,"$ARGV[1]" or die $!;
while(<IN>){
	chomp;
	my @A=split /\s+/;
	if($_=~/^#|^Chromosome/ || !exists $paralogous{$A[0]}{$A[1]}){
		print "$_\n";
	}
}
close IN;
