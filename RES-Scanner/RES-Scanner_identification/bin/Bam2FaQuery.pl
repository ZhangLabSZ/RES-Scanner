#!/usr/bin/perl -w
use strict;
use Getopt::Long;
die "Usage:perl $0  --sam <input.file.bam/sam[.gz]> --out <outfile> --ss [1|0] [--positive|--negative] --samtools <absolute_path_to_samtools>\n" unless @ARGV >= 7;
my (@sam,$ss,$outfile,$positive,$negative,$samtools);
$ss ||= 1;
GetOptions(
		"sam:s"=>\@sam,
		"ss:s"=>\$ss,
		"out:s"=>\$outfile,
		"positive"=>\$positive,
		"negative"=>\$negative,
		"samtools:s"=>\$samtools,
		);


die "Error: $samtools not exists!\n" unless -e $samtools;


open (OT,"> $outfile") or die $!;



if($ss){
	my $tag;
	if($positive && $negative){die;}
	if($positive){$tag="+";}
	elsif($negative){$tag="-";}
	&dealSam($sam[0],$tag,$ss);
}else{
	&dealSam($sam[0],"+",$ss);
}
close OT;

################################################    SUBROUTINE  ######################################################

sub dealSam {
	my ($file,$strand,$ss) = @_;
	if ($file =~ /\S+\.bam$/) {
		open (IN,"$samtools view $file |") or die $!;
	} elsif ($file =~ /\S+\.gz$/) {
		open (IN,"gunzip -c $file |") or die $!;
	}else {
		open (IN,$file) or die $!;
	}
	while (<IN>) {
		my @info = split /\s+/;
		next unless $_=~/NM:i:(\d+)/;
		my $mismath=$1;

		my $info_fq = &dealFq($info[0],$info[1]);
		my $info_strand = &dealStrand($info[0],$info[1]);
		my $start = $info[3];
		my $end = &getAlignment($start,$info[5]);
		my $tail = join ";",$info[2],$start,$end,$info[5],$info_strand,$mismath;

		if($ss){
		if ($strand eq "-") {
			unless ( ( ($info_fq == 1) && ($info_strand eq "+") )  || ( ($info_fq == 2) && ($info_strand eq "-") ) ) {
				die "Error: This($info[0]\t$info[1]\t) is wrong!";
			}
		} else {
			unless ( ( ($info_fq == 1) && ($info_strand eq "-") )  || ( ($info_fq == 2) && ($info_strand eq "+") ) ) {
				die "Error: This($info[0]\t$info[1]\t) is wrong!";
			}
		}
		}

		my $li = join "\/",$info[0],$info_fq;
		my $lin = join ";",$li,$tail;
		my $line = join "",">",$lin;
		print OT "$line\n";
		my $sequnce=$info[9];
		if($info_strand eq "-"){
			$sequnce=reverse $info[9];
			$sequnce=~tr/ATCGNatcgn/TAGCNtagcn/;
		}
		print OT "$sequnce\n";
	}
	close IN;
}
#####################################
sub dealFq {
	my ($id,$flag) = @_;
	my $fq;
	my $temp = $flag & 128;
	if ($temp == 0) {
		$fq = 1;
	} elsif ($temp == 128) {
		$fq = 2;
	} else {
		die "$id\t$flag\t doesn't exist!\n";
	}
	return $fq;
}
#####################################
sub dealStrand {
	my ($id,$flag) = @_;
	my $strand;
	my $temp = $flag & 16;
	if ($temp == 0) {
		$strand = "+";
	} elsif ($temp == 16) {
		$strand = "-";
	} else {
		die "$id\t$flag\t doesn't exist strand!\n";
	}
	return $strand;
}
######################################
sub getAlignment {
	my ($be,$cigar) = @_;
	my @M = $cigar =~ /(\d+[A-Z])/g;
	my $OFFSET=$be;
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
			die "the sixth column is '$cigar', sam format error!\n";
		}
	}
		return $OFFSET;
}
