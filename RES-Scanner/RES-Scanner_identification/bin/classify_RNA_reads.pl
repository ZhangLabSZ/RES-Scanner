#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'abs_path';
use File::Basename qw(basename);

my ($RNA_bam,$editTable,$samtools,$outdir);
my $goodRead ||= 1;
my $readType ||= 3;
my $phred ||= "33,33";
my $qual_cutoff ||= 30;
my $cut ||= "6,6";
my $paralogous_R ||= 1;
my $refined ||= 1;

GetOptions (
		"samtools:s"=>\$samtools,
		"RNA_bam:s"=>\$RNA_bam,
		"editTable:s"=>\$editTable,
		"readType:s"=>\$readType,
		"refined:s"=>\$refined,
		"refinedDepth:s"=>\$goodRead,
		"phred:s"=>\$phred,
		"qual_cutoff:s"=>\$qual_cutoff,
		"trim:s"=>\$cut,
		"outdir:s"=>\$outdir,
		"paralogous_R:s"=>\$paralogous_R,
		);

if( !$editTable || !$RNA_bam || !$samtools || !$outdir){
	print <<"Usage End.";

Options:
		--editTable     FILE    Table of candidate RNA editing sites.
		--readType      INT     The minimum number of unique RNA reads supporting editing for a candidate editing site. [3]
		--refined       NUM     Whether refined the number of RNA reads supporting candidate editing sites. 1 for yes, 0 for no. [1];
		--refinedDepth  INT     The minimum number of RNA reads in the middle of its length supporting editing for a candidate editing site.
		                        (e.g. from positions 23~68 of a 90-bp read). [1]
		--RNA_bam       FILE    The file.bam of RNA-seq data.
		--samtools      FILE    The absolute path of pre-installed SAMtools package.
		--phred         NUM     The Phred base quality for query QUALity of DNA.bam and RNA.bam files respectively, default DNA_ASCII-33,RNA_ASCII-33. [33,33]
		--qual_cutoff   NUM     Quality cutoff for BWA alignment. [30]
		--trim          INT     The number of bases solf-clipped at "5'end,3'end" of a read. [6,6]
		--paralogous_R  NUM     Remove candidate editing sites from those regions that are similar to other parts of the genome
		                        by BLAT alignment. 1 for yes, 0 for not. Note: force --blat. [1]
		--outdir        STR     Output directory.

		Example:
		perl $0 --editTable editTable.txt --RNA_bam RNA.bam --readType 3 --refinedDepth 1 --samtools samtools --phred 33,33 --qual_cutoff 30 --outdir ./outdir/ --paralogous_R 0 --refined 1

Usage End.
		exit;
}
if(!$refined && !$paralogous_R){
	die "Error: Do nothing for --paralogous_R 0 && --refined 0";
}

die "Error: RNA bam file $RNA_bam is not existent!" unless -e $RNA_bam;
mkdir $outdir unless -e $outdir;
$outdir=abs_path $outdir;
if($refined){
	my $editTable_filename=basename $editTable;
	if($editTable_filename=~/\.gz$/){
		open OUT1,"| gzip >$outdir/temp.$editTable_filename" or die $!;
	}else{
		open OUT1,">$outdir/temp.$editTable_filename" or die $!;
	}
}
if($paralogous_R){
	my $RNA_bam_filename=basename $RNA_bam;
	open OUT2,">$outdir/$RNA_bam_filename.editRead.fa" or die $!;
}
my ($phred_DNA,$phred_RNA)=split /,/,$phred;
my %site;
if($editTable=~/\.gz$/){
	open IN,"gunzip -c $editTable |" or die $!;
}else{
	open IN,"$editTable" or die $!;
}
while(<IN>){
	chomp;
	if ($_=~/^1.Chromosome_ID/ && $refined){
		print  OUT1 "$_\t15.Quadrant1/4_Read_Num(RNA)\t16.Quadrant2/3_Read_Num(RNA)\t17.Read_Type_Num(RNA)\n";
		next;
	}
	my @A=split /\t/;
	$site{$A[0]}{$A[1]}{info}=$_;
}
close IN;

###
my %site_sort;
my %log_total;
foreach my $chr (keys %site){
	my @array=sort {$a<=>$b} keys %{$site{$chr}};
	@{$site_sort{$chr}}=@array;
	my $n=@array;
	$log_total{$chr}=$n;
}


my $time= &getTime();
print STDERR "$time\tCalculating number of RNA reads that supporting RNA editing sites...\n";


##
if($RNA_bam=~/\.bam$/){
	open IN,"$samtools view $RNA_bam |" or die $!;
}elsif($RNA_bam=~/\.gz$/){
	open IN,"gunzip -c $RNA_bam |" or die $!;
}else{
	open IN,"$RNA_bam" or die $!;
}
while(<IN>){
	chomp;
	my @A=split /\t/;

	if(exists $site{$A[2]} && @{$site_sort{$A[2]}}>0){
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
		my $cutMaxLength=$cut5end>$cut3end ? $cut5end : $cut3end;
		my $read=substr($A[9],$cut5end,$readLen-$cut5end-$cut3end);
		my $qual=substr($A[10],$cut5end,$readLen-$cut5end-$cut3end);
		my $beg=$newStart;
		my $end= &offset($beg,$newMatchInfo);

		my @array= @{$site_sort{$A[2]}};
		my @next_array=@array;
		my $len=length($read);
		my $tailEnd_len=int $len/4;
		foreach my $pos (@array){
			if($pos+$cutMaxLength+1<$beg && $pos+$cutMaxLength+1<$A[3]){
				shift @next_array;
			}
			if($pos>$end || $array[-1]<$beg){last;}
			if($pos>=$beg && $pos<=$end){
				my @info_array=split /\t/,$site{$A[2]}{$pos}{info};
				my ($editType)=$info_array[10]=~/\w\-\>(\w)/;
				$editType=uc $editType;
				if($info_array[2] eq "-"){
					$editType=~tr/ATCGN/TAGCN/;
				}
				my $posOnRead=&locateOnRead($pos-$beg+1,$newMatchInfo);
				if($posOnRead==0){next;}
				my $readBase=uc (substr($read, $posOnRead-1, 1));
				my $readQual=substr($qual,$posOnRead-1,1);
				if($readBase ne $editType){next;}
				my $q=ord($readQual)-$phred_RNA;
				if($q<$qual_cutoff){next;}
				if($paralogous_R){
					my $bamStrand= ($A[1]&16)==16 ? "-" : "+";
					my $bamFq= ($A[1]&128)==128 ? "2" : "1";
					my $bamRead=uc $A[9];
					my $bamID="$A[0]/$bamFq;$A[2]:$pos";
					if($bamStrand eq "-"){
						$bamRead=reverse $bamRead;
						$bamRead=~tr/ATCGN/TAGCN/;
					}
					print OUT2 ">$bamID\n$bamRead\n";
				}

				if($refined){
					my $flag=1;
					if(!exists $site{$A[2]}{$pos}{lastBeg} && !exists $site{$A[2]}{$pos}{lastEnd} && !exists $site{$A[2]}{$pos}{lastScaf}){
					}elsif($site{$A[2]}{$pos}{lastBeg} eq $beg && $site{$A[2]}{$pos}{lastEnd} eq $end && $site{$A[2]}{$pos}{lastScaf} eq $A[2]){
						$flag=0;
					}
					my $beg_len=$posOnRead;
					my $end_len=$len-$posOnRead;
					if($beg_len< $tailEnd_len || $end_len< $tailEnd_len){
						$site{$A[2]}{$pos}{end}++;
					}else{
						$site{$A[2]}{$pos}{mid}++;
					}
					if($flag){$site{$A[2]}{$pos}{readNum}++;}
					($site{$A[2]}{$pos}{lastBeg},$site{$A[2]}{$pos}{lastEnd},$site{$A[2]}{$pos}{lastScaf})=($beg,$end,$A[2]);
				}
			}
		}
		@{$site_sort{$A[2]}}=@next_array;
		my $n=@array;
		my $next_n=@next_array;
		my $ratio=sprintf "%.2f",100-$next_n/$log_total{$A[2]}*100;
		if($next_n<$n){
			print STDERR "\t$A[2]: ${ratio}% sites finish.\n";
		}
	}
}
close IN;

###
undef %site_sort;

if($refined){
	foreach my $scaf (keys %site){
		foreach my $pos (keys  %{$site{$scaf}}){
			if(!exists $site{$scaf}{$pos}{end}){
				$site{$scaf}{$pos}{end}=0;
			}
			if(!exists $site{$scaf}{$pos}{mid}){
				$site{$scaf}{$pos}{mid}=0;
			}
			if(!exists $site{$scaf}{$pos}{readNum}){
				$site{$scaf}{$pos}{readNum}=0;
			}
			next if ($site{$scaf}{$pos}{mid}<$goodRead || $site{$scaf}{$pos}{readNum}<$readType);
			print OUT1 "$site{$scaf}{$pos}{info}\t$site{$scaf}{$pos}{end}\t$site{$scaf}{$pos}{mid}\t$site{$scaf}{$pos}{readNum}\n";
		}
	}
	close OUT1;
}
if($paralogous_R){
	close OUT2;
}
#########################################

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




sub locateOnRead{
	my $span2beg=shift;
	my $matchInfo=shift;
	my @M=$matchInfo=~/(\d+[A-Z])/g;
	my $offset=0;
	my $count=0;
	while($span2beg>0){
		if(!@M){return 0;}
		my $info=shift @M;
		$count++;
		if($info=~/M/){
			$info=~s/M//;
			if($info>=$span2beg){
				$offset+=$span2beg;
				$span2beg=0;
			}else{
				$offset+=$info;
				$span2beg-=$info;
			}
		}elsif($info=~/D/){
			$info=~s/D//;
			if($info>=$span2beg){
				return 0;
			}else{
				$span2beg-=$info;
			}
		}elsif($info=~/I/){
			$info=~s/I//;
			$offset+=$info;
		}elsif($info=~/S/){
			$info=~s/S//;
			if($count==1){$offset+=$info;}
			else{
				return 0;
			}
		}elsif($info=~/N/){
			$info=~s/N//;
			if($info>=$span2beg){
				return 0;
			}else{
				$span2beg-=$info;
			}
		}else{
			die "$matchInfo Error!";
		}
	}
	return $offset;
}


sub getTime{
	my $timestamp=time;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($timestamp);
	$sec=($sec<10) ? "0$sec" : $sec;
	$min=($min<10) ? "0$min" : $min;
	$hour=($hour<10) ? "0$hour" : $hour;
	$mday=($mday<10) ? "0$mday" : $mday;
	$mon=($mon<9) ? "0".($mon+1) : $mon+1;
	$year = $year + 1900;
	my $time=$year."-".$mon."-".$mday." ".$hour.":".$min.":"."$sec";
	return $time;
}



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

