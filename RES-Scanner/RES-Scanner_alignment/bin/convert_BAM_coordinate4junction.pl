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

die "perl $0 <junction2flankSequence.coordinate> <input.bam> <output.bam> <samtools>\n" unless @ARGV==4;
my %hash;
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,"$ARGV[0]" or die $!;
}
while(<IN>){
	chomp;
	my @A=split /\s+/;
	$hash{$A[0]}{chr}=$A[1];
	$hash{$A[0]}{junction}=$A[2];
}
close IN;


my $samtools=$ARGV[3];
die "Error: $samtools is not existent!\n" unless -e $samtools;

if($ARGV[1]=~/\.bam$/){
	open IN,"$samtools view -h $ARGV[1] |" or die $!;
}elsif($ARGV[1]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[1] |" or die $!;
}else{
	open IN,"$ARGV[1]" or die $!;
}

open OUT,"| $samtools view -bhS - > $ARGV[2]" or die $!;
#open OUT,"> $ARGV[2]" or die $!;
while(<IN>){
	chomp;
	$_=~s/\s+XA:Z:\S+;$//;
	my @A=split /\t/;
	if(@A<10){
		print OUT "$_\n";
		next;
	}
	next if $A[5] eq "*";

	if($A[6] ne "=" && $A[6] ne "." && exists $hash{$A[6]}){   # for the other pair 
		my $chrPE=$hash{$A[6]}{chr};
		my @B=split /,/,$hash{$A[6]}{junction};
		my $offset=$A[7];
		my @array;
		foreach my $region (@B){
			my @R=split /\-/,$region;
			for(my $i=$R[0];$i<=$R[1];$i++){
				push @array,$i;
			}
		}
		my $coordinate=$array[$offset-1];
		($A[6],$A[7])=($chrPE,$coordinate);
	}

	if(exists $hash{$A[2]}){
		my $chr=$hash{$A[2]}{chr};
		my @B=split /,/,$hash{$A[2]}{junction};
		my $offset=$A[3];
		my @array;
		my $flankSeqLength;
		foreach my $region (@B){
			my @R=split /\-/,$region;
			$flankSeqLength+=$R[1]-$R[0]+1;
			for(my $i=$R[0];$i<=$R[1];$i++){
				push @array,$i;
			}
		}
		my $end=&offset($A[3],$A[5]);
		if($flankSeqLength<$end){next;}
		my $coordinate=$array[$offset-1];

		if($A[6] eq "="){             # for the other pair
			$A[7]=$array[$A[7]-1];		
		}
#
		if($A[5] ne "*"){
			my @CIAGR;
			my @M=$A[5]=~/(\d+[A-Z])/g;
			my $OFFSET=$offset;
			for(my $i=0;$i<@M;$i++){
				if($M[$i]=~/M|D|N/){
					$M[$i]=~s/([MDN])//;
					my $tag=$1;
					my $n=0;
					for(my $j=0;$j<$M[$i];$j++){
						if($j>0){
							my $span=$array[$OFFSET+$j-1]-$array[$OFFSET+$j-2];
							if($span>1){
								my $Nlen=$span-1;
								push @CIAGR,"${n}$tag","${Nlen}N";
								$n=0;
							}
						}elsif($j==0){
							if($OFFSET>1 && $i>0 && $M[$i-1]!~/S/ && $array[$OFFSET-1]-$array[$OFFSET-2]>1){
								my $Nlen=$array[$OFFSET-1]-$array[$OFFSET-2]-1;
								push @CIAGR,"${Nlen}N";
							}
						}
						$n++;
					}
					if($n>0){
						push @CIAGR,"${n}$tag";
					}
					$OFFSET+=$M[$i];
				}elsif($M[$i]=~/S|I/){
					push @CIAGR,$M[$i];
				}else{
					die "the CIAGR is '$A[5]', sam format error!\n";
				}
			}
			$A[5]=join("",@CIAGR);
		}
#
		($A[2],$A[3])=($chr,$coordinate);
		print OUT join("\t",@A),"\n";
	}else{
		print OUT join("\t",@A),"\n";
	}
}
close IN;
close OUT;


sub offset{
	my $beg=shift;
	my $matchInfo=shift;
	my @M=$matchInfo=~/(\d+[A-Z])/g;
	my $offset=$beg;
	for(my $i=0;$i<@M;$i++){
		if($M[$i]=~/M/){
			$M[$i]=~s/M//;
			$offset=$offset+$M[$i]-1;
		}elsif($M[$i]=~/D/){
			$M[$i]=~s/D//;
			$offset=$offset+$M[$i]+1;
		}elsif($M[$i]=~/I/){
			$offset=$offset+1;
		}elsif($M[$i]=~/S/){
		}elsif($M[$i]=~/N/){
			$M[$i]=~s/N//;
			$offset=$offset+$M[$i]+1;
		}else{
			die "$matchInfo Error!";
		}
	}
	return $offset;
}




