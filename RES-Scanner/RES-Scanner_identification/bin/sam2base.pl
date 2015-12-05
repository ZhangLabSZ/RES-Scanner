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
use Cwd 'abs_path';
use File::Basename qw(dirname);
use Getopt::Long;
die "
################################## Usage ##################################
Usage: $0 <genome.fa> <sorted.sam> <outfile.gz> --trim [0,0] --samtools <abs_path samtools>
Attention: The file.sam should be sorted by chromosome ID and leftmost mapping position!
################################## Usage ##################################
\n" unless @ARGV>=5;

my $cut ||= "0,0";
my $samtools;
GetOptions(
		"trim:s"=>\$cut,
		"samtools:s"=>\$samtools,
		);

die "Error: $samtools not exists!\n" unless -e $samtools;
if($cut!~/^\d+,\d+$/){
	die "Error: --trim option should be two numbers separated by a comma, i.e. '6,6'.\n";
}

if($ARGV[2]!~/\.gz$/){die "$ARGV[2] should be the file.gz format!\n";}
if($ARGV[0]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[0] |" or die $!;
}else{
	open IN,"$ARGV[0]" or die $!;
}
$/=">";
<IN>;
my %genome;
my %hash;
while(<IN>){
	chomp;
	$_=~/(.+?)\n/;
	my $id=(split /\s+/,$1)[0];
	$_=~s/.+?\n//;
	$_=~s/\s+//g;
	$genome{$id}=$_;
	$hash{$id}=1;
}
$/="\n";
close IN;

my $outfile=abs_path($ARGV[2]);
my $outdir=dirname($outfile);

my $genomeID;
if($ARGV[1]=~/\.gz$/){
	open IN,"gunzip -c $ARGV[1] |" or die $!;
}elsif($ARGV[1]=~/\.bam$/){
	open IN,"$samtools view  $ARGV[1] |" or die $!;
}else{
	open IN,"$ARGV[1]" or die $!;
}

my %storage;
my (%insert,%delete);
my $offset=1;
my %sortedORnot;

open OUT,"| gzip >$outfile" or die $!;
while(<IN>){
	chomp;
	my @A=split /\t/;
	if(exists $sortedORnot{$A[2]}){die "sam file is not sorted!";}
	if(!exists $genome{$A[2]}){next;}
	if(!defined $genomeID){$genomeID=$A[2];delete $hash{$A[2]};}
##edit on 20120304
	next unless $A[5] ne "*";
	my @M=$A[5]=~/(\d+[A-Z])/g;
	my $OFFSET=$A[3];
	for(my $i=0;$i<@M;$i++){
		if($M[$i]=~/M/){
			$M[$i]=~s/M//;
			$OFFSET+=$M[$i]-1;
		}elsif($M[$i]=~/D/){
			$M[$i]=~s/D//;
			$OFFSET+=$M[$i]+1;
		}elsif($M[$i]=~/I/){
			$M[$i]=~s/I//;
			$OFFSET++;
		}elsif($M[$i]=~/S/){
		}elsif($M[$i]=~/N/){
			$M[$i]=~s/N//;
			$OFFSET+=$M[$i]+1;
		}else{
			die "the sixth column is '$A[5]', sam format error!\n";			
		}
	}

##
	if(length($genome{$A[2]})<$OFFSET){
		next;
	}
	while($genomeID ne $A[2]){
		$sortedORnot{$genomeID}=1;
		for(my $i=$offset;$i<=length($genome{$genomeID});$i++){ 
			my $base=substr($genome{$genomeID},$i-1,1);
			if(exists $storage{$genomeID}{$i}){
				my $dep=length($storage{$genomeID}{$i}{seq});
				my $seqString="[M]:$storage{$genomeID}{$i}{seq}";
				my $qualString="[M]:$storage{$genomeID}{$i}{qual}";
				if(exists $insert{$genomeID}{$i}){
					$seqString .= ";[I]:".join(",",@{$insert{$genomeID}{$i}{Iseq}});
					$qualString .= ";[I]:".join(",",@{$insert{$genomeID}{$i}{Iqual}});
					delete $insert{$genomeID}{$i};
				}
				if(exists $delete{$genomeID}{$i}){
					$seqString .= ";[D]:".join(",",@{$delete{$genomeID}{$i}{Dseq}});
					$qualString .= ";[D]:".join(",",@{$delete{$genomeID}{$i}{Dqual}});
					delete $delete{$genomeID}{$i};
				}			
				print OUT "$genomeID\t$i\t$base\t$seqString\t$qualString\t$dep\n";
				delete $storage{$genomeID}{$i};
			}
		}
		$genomeID=$A[2];
		delete $hash{$A[2]};
		$offset=1;
	}
	if($A[3]>$offset){
		for(my $i=$offset;$i<$A[3];$i++){                   # delete the hash key which will not use again to release memory.
			my $base=substr($genome{$genomeID},$i-1,1);
			if(exists $storage{$genomeID}{$i}){
				my $dep=length($storage{$genomeID}{$i}{seq});
				my $seqString="[M]:$storage{$genomeID}{$i}{seq}";
				my $qualString="[M]:$storage{$genomeID}{$i}{qual}";
				if(exists $insert{$genomeID}{$i}){
					$seqString .= ";[I]:".join(",",@{$insert{$genomeID}{$i}{Iseq}});
					$qualString .= ";[I]:".join(",",@{$insert{$genomeID}{$i}{Iqual}});
					delete $insert{$genomeID}{$i};
				}
				if(exists $delete{$genomeID}{$i}){
					$seqString .= ";[D]:".join(",",@{$delete{$genomeID}{$i}{Dseq}});
					$qualString .= ";[D]:".join(",",@{$delete{$genomeID}{$i}{Dqual}});
					delete $delete{$genomeID}{$i};
				}			
				print OUT "$genomeID\t$i\t$base\t$seqString\t$qualString\t$dep\n";
				delete $storage{$genomeID}{$i};			
			}
		}
	}

	$offset=$A[3];

#	add for trim reads begin
	my $readLen=length($A[9]);
	my $qualLen=length($A[10]);
	if($readLen != $qualLen){next;}
    my ($newStart,$newMatchInfo);

	my $cutinfo=$cut;
	if( ($A[1]&16) == 16){
		my @cut_array=split /,/,$cutinfo;
		$cutinfo="$cut_array[1],$cut_array[0]";
	}

	($newStart,$newMatchInfo)=&start_MatchInfo($A[3],$A[5],$cutinfo);

	my ($cut5end,$cut3end)=split /,/,$cutinfo;
    my $read=substr($A[9],$cut5end,$readLen-$cut5end-$cut3end);
	my $qual=substr($A[10],$cut5end,$readLen-$cut5end-$cut3end);

	my @seq=split //,$read;
	my @qual=split //,$qual;

#   add for trim reads end
	@M=();
	@M=$newMatchInfo=~/(\d+[A-Z])/g;
	my $refPos=$newStart;
	my $seqPos=0;
	for(my $i=0;$i<@M;$i++){
		if($M[$i]=~/M/){
			$M[$i]=~s/M//;
			for(my $j=0;$j<$M[$i];$j++){
				$storage{$A[2]}{$refPos+$j}{seq} .= $seq[$seqPos];
				$storage{$A[2]}{$refPos+$j}{qual} .= $qual[$seqPos];
				$seqPos++;
			}
			$refPos+=$M[$i];
		}elsif($M[$i]=~/D/){
			$M[$i]=~s/D//;
			my $Dqual;
			for(my $j=0;$j<$M[$i];$j++){$Dqual .= 0;}
			my $Dseq=substr($genome{$genomeID},$refPos-1,$M[$i]);
			if($M[$i]>0){
				push ( @{$delete{$A[2]}{$refPos}{Dseq}},$Dseq );
				push ( @{$delete{$A[2]}{$refPos}{Dqual}},$Dqual );
			}
			$refPos+=$M[$i];
		}elsif($M[$i]=~/I/){
			$M[$i]=~s/I//;
			my ($Iseq,$Iqual);
			for(my $j=0;$j<$M[$i];$j++){
				$Iseq .= $seq[$seqPos];
				$Iqual .= $qual[$seqPos];
				$seqPos++;
			}
			if($M[$i]>0){
				push ( @{$insert{$A[2]}{$refPos-1}{Iseq}},$Iseq);
				push ( @{$insert{$A[2]}{$refPos-1}{Iqual}},$Iqual);
			}
		}elsif($M[$i]=~/S/){
			$M[$i]=~s/S//;
			$seqPos+=$M[$i];
		}elsif($M[$i]=~/N/){
			$M[$i]=~s/N//;
			$refPos+=$M[$i];
		}else{
			die "the sixth column is '$A[5]', sam format error!\n";			
		}
	}
}

if(%storage){							 
	for(my $i=$offset;$i<=length($genome{$genomeID});$i++){ 
		my $base=substr($genome{$genomeID},$i-1,1);
		if(exists $storage{$genomeID}{$i}){
			my $dep=length($storage{$genomeID}{$i}{seq});
			my $seqString="[M]:$storage{$genomeID}{$i}{seq}";
			my $qualString="[M]:$storage{$genomeID}{$i}{qual}";
			if(exists $insert{$genomeID}{$i}){
				$seqString .= ";[I]:".join(",",@{$insert{$genomeID}{$i}{Iseq}});
				$qualString .= ";[I]:".join(",",@{$insert{$genomeID}{$i}{Iqual}});
				delete $insert{$genomeID}{$i};
			}
			if(exists $delete{$genomeID}{$i}){
				$seqString .= ";[D]:".join(",",@{$delete{$genomeID}{$i}{Dseq}});
				$qualString .= ";[D]:".join(",",@{$delete{$genomeID}{$i}{Dqual}});
				delete $delete{$genomeID}{$i};
			}			
			print OUT "$genomeID\t$i\t$base\t$seqString\t$qualString\t$dep\n";
			delete $storage{$genomeID}{$i};			
		}
	}
	$offset=1;
}
close IN;
close OUT;


sub start_MatchInfo{
	my $start=shift;
	my $matchInfo=shift;
	my $cutlen=shift;
	my ($cutlen5,$cutlen3)=split /,/,$cutlen;
	my @M=$matchInfo=~/(\d+[A-Z])/g;
##for head
	my $headcutlen=0;
	my $offset=0;
	my $tag;
	while($headcutlen<$cutlen5){
		my $info=shift @M;
		if($info=~/M/){
			$info=~s/M//;
			if($headcutlen+$info>=$cutlen5){
				my $diff=$info-($cutlen5-$headcutlen);
				if($diff>0){$tag="${diff}M";}
				$offset+=($cutlen5-$headcutlen);
				$headcutlen=$cutlen5;
			}else{
				$headcutlen+=$info;
				$offset+=$info;
			}
		}elsif($info=~/D|N/){
			$info=~s/D|N//;
			$offset+=$info;
		}elsif($info=~/I|S/){
			$info=~s/I|S//;
			if($headcutlen+$info>=$cutlen5){
				my $diff=$info-($cutlen5-$headcutlen);
				if($diff>0){$tag="${diff}S";}   # "S" for edge
					$headcutlen=$cutlen5;
			}else{
				$headcutlen+=$info;
			}
		}else{
			die "Error matchInfo: $matchInfo\n";
		}
	}
	if($tag){
		unshift (@M,$tag);
	}
	my $newStart=$start+$offset;
	my $firstTag=$M[0];
	if($firstTag=~/D|N/){
		$firstTag=~s/D|N//;
		$newStart+=$firstTag;
		shift @M;
	}elsif($firstTag=~/I/){
		$firstTag=~s/I/S/;
		$M[0]=$firstTag;
	}
#### for tail
	my $tailcutlen=0;
	$tag="";
	while($tailcutlen<$cutlen3){
		my $info=pop @M;
		if($info=~/M/){
			$info=~s/M//;
			if($tailcutlen+$info>=$cutlen3){
				my $diff=$info-($cutlen3-$tailcutlen);
				if($diff>0){$tag="${diff}M";}
				$tailcutlen=$cutlen3;
			}else{
				$tailcutlen+=$info;
			}
		}elsif($info=~/D|N/){
		}elsif($info=~/I|S/){
			$info=~s/I|S//;
			if($tailcutlen+$info>=$cutlen3){
				my $diff=$info-($cutlen3-$tailcutlen);
				if($diff>0){$tag="${diff}S";}   # "S" for edge
					$tailcutlen=$cutlen3;
			}else{
				$tailcutlen+=$info;
			}
		}else{
			die "Error matchInfo: $matchInfo\n";
		}
	}
	if($tag){
		push (@M,$tag);
	}
	my $lastTag=$M[-1];
	if($lastTag=~/D|N/){
		pop @M;
	}elsif($lastTag=~/I/){
		$lastTag=~s/I/S/;
		$M[-1]=$lastTag;
	}

######
	my $newMatchInfo=join("",@M);
	return ($newStart,$newMatchInfo);
}

