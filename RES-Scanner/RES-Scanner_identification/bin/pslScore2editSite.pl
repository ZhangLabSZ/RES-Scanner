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
my ($pslScore,$junctionCoordinate,$editTable,$bestHitRatio);
$bestHitRatio ||= 0.5;

GetOptions(
		"junctionCoordinate:s"=>\$junctionCoordinate,
		"pslScore:s"=>\$pslScore,
		"editTable:s"=>\$editTable,
		"bestHitRatio:s"=>\$bestHitRatio,
		);

if(!$pslScore || !$editTable){
	die "Usage: perl $0 --pslScore <psl.score> --junctionCoordinate [junctionCoordinate.file] --editTable [editTable] --bestHitRatio [0.5]\n";
}

die "Error: $pslScore is not existent!" unless -e $pslScore;
die "Error: $editTable is not existent!" unless -e $editTable;

my %junc;
if(defined $junctionCoordinate){
	die "Error: $junctionCoordinate is not existent!" unless -e $junctionCoordinate;
	if($junctionCoordinate=~/\.gz$/){
		open IN,"gunzip -c $junctionCoordinate |" or die $!;
	}else{
		open IN,"$junctionCoordinate" or die $!;
	}
	while(<IN>){
		chomp;
		my @A=split /\s+/;
		$junc{$A[0]}{chr}=$A[1];
		$junc{$A[0]}{junction}=$A[2];
	}
	close IN;
}

my %hash;
my @array;
my $lastID;
my $lastKey;
open IN,"$pslScore" or die $!;
while(<IN>){
	chomp;
	my @A=split /\s+/;
	if(defined $junctionCoordinate && exists $junc{$A[0]}){
		my $chrID=$junc{$A[0]}{chr};
		my @B=split /,/,$junc{$A[0]}{junction};
		my @AR;
		foreach my $region (@B){
			my @R=split /\-/,$region;
			for(my $i=$R[0];$i<=$R[1];$i++){
				push @AR,$i;
			}
		}
		my $beg=$AR[$A[1]];   #zero based half open of *.psl
			my $end=$AR[$A[2]-1];
		($A[0],$A[1],$A[2])=($chrID,$beg,$end);
	}

	my ($queryID)=$A[3]=~/^(\S+):\d+\-\d+$/;
	my ($editChr,$editPos,$bamBeg)=$queryID=~/^\S+;(\S+):(\d+):(\d+)$/;
	if(!defined $lastID){
		$lastID=$queryID;
		$lastKey="$editChr,$editPos";
	}else{
		if($lastID ne $queryID || $lastKey ne "$editChr,$editPos"){
			my $editFlag= &assess_pslScore(\@array);
			$hash{$lastKey}{edit}+=$editFlag;
			$lastID=$queryID;
			$lastKey="$editChr,$editPos";
			@array=();
		}
	}
	push @array,[$A[0],$A[1],$A[2],$queryID,$A[4],$editChr,$editPos,$bamBeg];
}
close IN;

if(@array){
	my ($editChr,$editPos,$bamBeg)=$lastID=~/^\S+;(\S+):(\d+):(\d+)$/;
	my $editFlag= &assess_pslScore(\@array);
	$hash{$lastKey}{edit}+=$editFlag;
}

if($editTable=~/\.gz$/){
	open IN,"gunzip -c $editTable |" or die $!;
}else{
	open IN,"$editTable" or die $!;
}
while(<IN>){
	chomp;
	if($_=~/^Chromosome|^1.Chromosome_ID/){
		print "$_\t18.Extracted_Read_Num(Blat)\t19.BestHit_Read_Num(Blat)\t20.BestHit_Read_Ratio(Blat)\n";
		next;
	}
	my @A=split /\t/;
	my @B=split /,/,$A[9];
	my %baseHash=("A",$B[0],"C",$B[1],"G",$B[2],"T",$B[3]);
	my ($editBase)=$A[10]=~/\w->(\w)/;
	if($A[2] eq "-"){
		$editBase=~tr/ATCG/TAGC/;
	}
	my $totalEdit=$baseHash{$editBase};

	my $key="$A[0],$A[1]";
	if(!exists $hash{$key}){
		$hash{$key}{edit}=0;
	}
	my $ratio=sprintf "%.4f",$hash{$key}{edit}/$totalEdit;
	if($ratio>=$bestHitRatio ){
		print "$_\t$totalEdit\t$hash{$key}{edit}\t$ratio\n";
	}
}
close IN;

############################################################################################
sub assess_pslScore{
	my $array_ref=shift;
	my @array=@$array_ref;
	@array=sort{$b->[4] <=> $a->[4]} @array;
	my $bestScore=$array[0][4];
	if($bestScore==0){return 0;}
	foreach my $a (@array){
		my ($tID,$tbeg,$tend,$qID,$score,$editChr,$editPos,$bamBeg)=@$a;
#		if( $tID eq $editChr && $bamBeg==$tbeg+1  && $editPos>=$tbeg && $editPos<=$tend){
		if( $tID eq $editChr && $editPos>=$tbeg && $editPos<=$tend){	
		}else{
			if($score/$bestScore>=0.95){
				return 0;
			}
		}
	}
	return 1;
}


