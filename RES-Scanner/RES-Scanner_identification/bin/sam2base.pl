#!/usr/bin/perl -w
use strict;
use Cwd 'abs_path';
use File::Basename qw(dirname);
use Getopt::Long;
die "
################################## Usage ##################################
Usage: $0 <genome.fa> <sorted.sam> <outfile.gz> --trim [0] --samtools <abs_path samtools>
Attention: The file.sam should be sorted by chromosome ID and leftmost mapping position!
################################## Usage ##################################
\n" unless @ARGV>=5;

my $cut ||= 0;
my $samtools;
GetOptions(
		"trim:s"=>\$cut,
		"samtools:s"=>\$samtools,
		);

die "Error: $samtools not exists!\n" unless -e $samtools;

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

#	add for cutoff reads begin
	my $readLen=length($A[9]);
	my $qualLen=length($A[10]);
	if($readLen != $qualLen){next;}
    my ($newStart,$newMatchInfo);
	($newStart,$newMatchInfo)=&start_MatchInfo($A[3],$A[5],$cut);

    my $read=substr($A[9],$cut,$readLen-2*$cut);
	my $qual=substr($A[10],$cut,$readLen-2*$cut);

	my @seq=split //,$read;
	my @qual=split //,$qual;

#   add for cutoff reads end
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
			push ( @{$delete{$A[2]}{$refPos}{Dseq}},$Dseq );
			push ( @{$delete{$A[2]}{$refPos}{Dqual}},$Dqual );
			$refPos+=$M[$i];
		}elsif($M[$i]=~/I/){
			$M[$i]=~s/I//;
			my ($Iseq,$Iqual);
			for(my $j=0;$j<$M[$i];$j++){
				$Iseq .= $seq[$seqPos];
				$Iqual .= $qual[$seqPos];
				$seqPos++;
			}
			push ( @{$insert{$A[2]}{$refPos-1}{Iseq}},$Iseq);
			push ( @{$insert{$A[2]}{$refPos-1}{Iqual}},$Iqual);
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
	my @M=$matchInfo=~/(\d+[A-Z])/g;
##for head
	my $headcutlen=0;
	my $offset=0;
	my $tag;
	while($headcutlen<$cutlen){
		my $info=shift @M;
		if($info=~/M/){
			$info=~s/M//;
			if($headcutlen+$info>=$cutlen){
				my $diff=$info-($cutlen-$headcutlen);
				if($diff>0){$tag="${diff}M";}
				$offset+=($cutlen-$headcutlen);
				$headcutlen=$cutlen;
			}else{
				$headcutlen+=$info;
				$offset+=$info;
			}
		}elsif($info=~/D|N/){
			$info=~s/D|N//;
			$offset+=$info;
		}elsif($info=~/I|S/){
			$info=~s/I|S//;
			if($headcutlen+$info>=$cutlen){
				my $diff=$info-($cutlen-$headcutlen);
				if($diff>0){$tag="${diff}S";}   # "S" for edge
					$headcutlen=$cutlen;
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
	while($tailcutlen<$cutlen){
		my $info=pop @M;
		if($info=~/M/){
			$info=~s/M//;
			if($tailcutlen+$info>=$cutlen){
				my $diff=$info-($cutlen-$tailcutlen);
				if($diff>0){$tag="${diff}M";}
				$tailcutlen=$cutlen;
			}else{
				$tailcutlen+=$info;
			}
		}elsif($info=~/D|N/){
		}elsif($info=~/I|S/){
			$info=~s/I|S//;
			if($tailcutlen+$info>=$cutlen){
				my $diff=$info-($cutlen-$tailcutlen);
				if($diff>0){$tag="${diff}S";}   # "S" for edge
					$tailcutlen=$cutlen;
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

