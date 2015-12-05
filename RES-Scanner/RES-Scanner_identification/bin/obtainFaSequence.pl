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
die "Usage: perl $0 <*.list> <*1.fq> [<*2.fq>]\n" unless @ARGV >= 2;

my %hash;
if ($ARGV[0] =~ /\.gz$/) {
	open (IN,"gunzip -c $ARGV[0] |") or die $!;
} else {
	open (IN,$ARGV[0]) or die $!;
}
while (<IN>) {
	chomp;
	my ($id,$tail) = split /;/;
	$hash{$id} = $_;
}
close IN;

&dealFq($ARGV[1],\%hash,1);

if(@ARGV==3){
	&dealFq($ARGV[2],\%hash,2);
}



############################################	SUBROUTINE	###################################################
sub dealFq {
	my ($file,$ref,$tag) = @_;
	if ($file =~ /\S+\.gz/) {
		open (IN,"gunzip -c $file |") or die $!;
	} else {
		open (IN,$file) or die $!;
	}
	while (<IN>) {
		chomp;
		my $fqID=(split /\s+/)[0];
		$fqID=~s/\/1$|\/2$//;
		$fqID=$fqID."/$tag";
		if (exists $ref->{$fqID}) {
			my $id = $ref->{$fqID};
			$id =~ s/@/>/;
			print "$id\n";
			my $line = <IN>;
			print "$line";
			<IN>;
			<IN>;
		} else {
			<IN>;
			<IN>;
			<IN>;
		}
	}
	close IN;
}
