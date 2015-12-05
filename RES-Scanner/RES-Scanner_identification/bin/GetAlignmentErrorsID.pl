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
die "Usage: perl $0 <*.fq.psl.gz> [outfile *.list]\n" unless @ARGV == 1;

my (%reads,%not);
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
	if(@query!=7){die "Error: ID contatins semicolon ';' !\n";}
	my ($scaffold,$st,$ed,$cigar,$strand,$mismatch) = @query[1..6];
	next unless $info[0]>=($ed-$st+1)*0.9;
	my $head = $st-1;
	my $tail = $ed;
	if ($info[1] <= $mismatch) {
#		if (($info[8] ne $strand) || ($info[13] ne $scaffold)) {
		if ($info[13] ne $scaffold) {
			if (! exists $not{$query[0]}) {
				print "\@$info[9]\n";
				$not{$query[0]} = 1;
				next;
			}
		} else {
			if (($info[15] >= $tail) || ($info[16] <= $head)) {
				if (! exists $not{$query[0]}) {
					print "\@$info[9]\n";
					$not{$query[0]} = 1;
					next;
				}
			}
		}
	}	
}
close IN;

