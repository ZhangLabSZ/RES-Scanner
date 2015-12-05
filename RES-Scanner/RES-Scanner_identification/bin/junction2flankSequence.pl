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
die "perl $0 <genome.fa> <junction.pos> [readLength]\n" unless @ARGV==3;
my $readLen=$ARGV[2];
my $id;
my @Array;

my %length;
my $fa_file = $ARGV[0];

if ($fa_file =~ /\.gz$/) {
	open IN, "gunzip -c $fa_file |";
} else {
	open IN, $fa_file;
}
$/ = ">";
<IN>;
while (<IN>) {
	chomp;
	my $scaf;
	if (/(\S+)/) {
		$scaf = $1;
	}
	s/.+\n//;
	s/\s+//g;
	my $len = length($_);
	$length{$scaf}=$len;
}
$/="\n";
close IN;

if($ARGV[1]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[1] |" or die $!;
}else{
	open IN,"$ARGV[1]" or die $!;
}
while(<IN>){
	chomp;
	my @A=split /\s+/;
	my ($ID)=$A[1];
	if(!defined $id){
		$id=$ID;
	}else{
		if($id ne $ID){
			&junction2flankSequence(\@Array);
			@Array=();
			$id=$ID;
		}
	}
	push @Array,[@A];
}
close IN;

if(@Array){
	&junction2flankSequence(\@Array);
}

#############################################################################################
sub junction2flankSequence{
	my $arr=shift;
	my @array=@$arr;
	@array=sort{$a->[3]<=>$b->[3]} @array;
	for(my $i=0;$i<@array;$i++){
		my $up_index=$i-1;
		my $down_index=$i+1;
		my ($up_size,$down_size)=(0,0);
		my (@up_array,@down_array);
		while($up_index>=0){
			if($up_size>=$readLen-1){last;}
			my ($beg,$end);
			if($array[$up_index+1][3]<=$array[$up_index][4]){
				$up_index--;
				next;
			}
			if($array[$up_index+1][3]-1-($array[$up_index][4]+1)+1+$up_size >= $readLen-1){
				$beg=$array[$up_index+1][3]-($readLen-1-$up_size);
				$end=$array[$up_index+1][3]-1;
			}else{
				$beg=$array[$up_index][4]+1;
				$end=$array[$up_index+1][3]-1;
			}
			unshift @up_array,"$beg-$end";
			$up_size+=$end-$beg+1;
			$up_index--;
		}
		if($up_size>$readLen-1){die;}
		if($up_size<$readLen-1 && $up_index==-1){
			my ($beg,$end);
			if($array[$up_index+1][3]-1>=$readLen-1-$up_size){
				$beg=$array[$up_index+1][3]-($readLen-1-$up_size);
				$end=$array[$up_index+1][3]-1;
			}else{
				$beg=1;
				$end=$array[$up_index+1][3]-1;
			}
			unshift @up_array,"$beg-$end";
		}

##
		while($down_index<=@array-1){
			if($down_size>=$readLen-1){last;}
			my ($beg,$end);
			if($array[$down_index-1][4]>=$array[$down_index][3]){
				$down_index++;
				next;
			}
			if($array[$down_index][3]-1-($array[$down_index-1][4]+1)+1+$down_size >= $readLen-1){
				$beg=$array[$down_index-1][4]+1;
				$end=$array[$down_index-1][4]+($readLen-1-$down_size);
			}else{
				$beg=$array[$down_index-1][4]+1;
				$end=$array[$down_index][3]-1;
			}
			push @down_array,"$beg-$end";
			$down_size+=$end-$beg+1;
			$down_index++;
		}
		if($down_size>$readLen-1){die;}
		if($down_size<$readLen-1 && $down_index==@array){
			my ($beg,$end);
			if($length{$array[$down_index-1][1]}-$array[$down_index-1][4]>=$readLen-1-$down_size){
				$beg=$array[$down_index-1][4]+1;
				$end=$array[$down_index-1][4]+$readLen-1-$down_size;
			}else{
				$beg=$array[$down_index-1][4]+1;
				$end=$length{$array[$down_index-1][1]};
			}
			push @down_array,"$beg-$end";
		}
##
		my $region=join(",",@up_array,@down_array);
		my $regionID="$array[$i][0]";
		print "$regionID\t$array[$i][1]\t$region\n";
	}
}
