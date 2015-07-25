#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use Cwd 'abs_path';

my ($RNA_singleBase,$DNA_singleBase,$DNAdepth,$phred,$qual_cutoff,$RNAdepth,$editLevel,$editDepth,$strand,$RNA_bam,$ploidy,$samtools,$genome,$method);

my $HomoPrior ||= 0.99;
my $HetePrior ||= 0.01;
my $rate ||= 4; #the rate of transition over transversion
#$DNAdepth ||= 7;
#$RNAdepth ||= 3;
$phred ||= 33;
$qual_cutoff ||= 30;
#$editLevel ||= 0.05;
#$editDepth ||= 3;
$ploidy ||= 2;
$method ||= "Bayesian";

GetOptions (
		"RNA_singleBase:s"=>\$RNA_singleBase,
		"DNA_singleBase:s"=>\$DNA_singleBase,
		"RNA_bam:s"=>\$RNA_bam,
		"genome:s"=>\$genome,
#		"DNAdepth:s"=>\$DNAdepth,
		"phred:s"=>\$phred,
		"qual_cutoff:s"=>\$qual_cutoff,
#		"RNAdepth:s"=>\$RNAdepth,
#		"editLevel:s"=>\$editLevel,
#		"editDepth:s"=>\$editDepth,
		"strand:s"=>\$strand,
		"ploidy:s"=>\$ploidy,
		"samtools:s"=>\$samtools,
		"HomoPrior:s"=>\$HomoPrior,
		"HetePrior:s"=>\$HetePrior,
		"rate:s"=>\$rate,
		"method:s"=>\$method,
		);

if (!$RNA_singleBase || !$DNA_singleBase || !$RNA_bam || !$strand || !$samtools ) {
	print <<"Usage End.";
Description:
	This script is used to obtain the homo-site for RNA editing set.
		--RNA_singleBase    The raw file of RNAedinting_sam2base
		--DNA_singleBase    The single base file for whole genome.
		--RNA_bam           The file.bam of RNA-seq data.
		--genome            The refference genome (Required when --method Bayesian).
		--phred             ASCII-33 gives the Phred base quality for query QUALity of BAM file. [33]
		--qual_cutoff       Quality cutoff for BWA alignment. [30]
		--strand            The strand of data, '+' or '-' for strand specific data, 'unknown' for non-strand specific data.
		--ploidy            The ploidy level, 1 for monoploid , 2 for diploid, 3 for triploid, 4 for tetraploid, and so on. [2].
		--samtools          The absolute path of pre-installed SAMtools package;
		--method            Method for detecting SNPs.'Bayesian' or 'Binomial' or 'Stringent'. [Bayesian]
		--HomoPrior         The prior probability of homozygous genomic positions. (force --method Bayesian) [0.99]
		--HetePrior         The prior probability of heterozygous genomic positions. (force --method Bayesian) [0.01]
		--rate              The rate of transitions over transversions of the genome. (force --method Bayesian) [4]


Example:
	perl $0 --RNA_singleBase RNA.sam2base --DNA_singleBase DNA.sam2base --RNA_bam RNA.bam --phred 33 --qual_cutoff 30 --strand unknown --ploidy 2 --samtools /abs_path/samtools --method Bayesian --genome genome.fa

Usage End.
		exit;
}

die "$RNA_singleBase does not exist!\n" unless -e $RNA_singleBase;
die "$DNA_singleBase doesn't exist!\n" unless -e $DNA_singleBase;
die "$samtools doesn't exist!\n" unless -e $samtools;
$samtools=abs_path $samtools;

if($strand ne "+" && $strand ne "-" && $strand ne "unknown"){
	die "Input Error: --strand should be '+' or '-' or 'unknown'!\n";
}

if($method ne "Bayesian" && $method ne "Binomial" && $method ne "Stringent"){
	die "Error: unknown options value '$method' for --Method [Bayesian|Binomial|Stringent]\n";
}elsif($method eq "Bayesian" && !$genome){
	die "Error: options --genome is undefined in Bayesian mode\n";
}


if ($RNA_singleBase =~ /\S+\.gz$/) {
	open (IN,"gunzip -c $RNA_singleBase |") or die $!;
} else {
	open (IN,$RNA_singleBase) or die $!;
}
my $time= &getTime();
print STDERR "$time\tProccessing candidate RNA edit site ...\n";

my %RNAediting;
while (<IN>) {
	chomp;
	my @info = split /\s+/;
#	next unless $info[5] >= $RNAdepth;
	next unless $info[5] >= 3;   #optional
	die "These bases($info[0]\t$info[1]) doesn't match!\n" unless ($info[3] =~ /^\[M\]:([ATCGNatcgn]+)(;\[I\]:[ATCGNatcgn,]+)?(;\[D\]:[ATCGNatcgn,]+)?/);
	my @bases = split //,$1;
	my $number = $#bases+1;
	my @qualities;
	my @aa = split //,$info[4];
	for (my $k=4;$k<$number+4;$k++) {
		push @qualities,$aa[$k];
	}
	die "Something was wrong at the location[$info[0]\t$info[1]]!\n" unless @bases == @qualities;
	my (@b_temp,@q_temp);
	for (my $i=0;$i<@bases;$i++) {
		my $score = ord($qualities[$i]) - $phred;
		if ($score >= $qual_cutoff) {
			push @b_temp,$bases[$i];
			push @q_temp,$qualities[$i];
		}
	}
	my $coverage = @b_temp;
	next unless $coverage >= 3;    #optional
	my $Bases = join "","[M]:",@b_temp;
	my $Quals = join "","[M]:",@q_temp;
	my %editHash;
	my $ref_base = uc($info[2]);
	foreach my $base (@b_temp) {
		$base = uc($base);
		next if ($base eq $ref_base);
		$editHash{$base} += 1;
	}
	my $flag = 0;
	foreach my $key (sort keys %editHash) {
		$flag = 1 if ($editHash{$key}/$coverage > 0);   #optional
	}
	if ($flag == 1) {
		my $tail = "$ref_base\t$Bases\t$Quals\t$coverage";
		$RNAediting{$info[0]}{$info[1]} = $tail;
	}
##	
}
close IN;

$time= &getTime();
print STDERR "$time\tProccessing homozygous DNA sites ...\n";

my %BaseContent;
if($method eq "Bayesian"){
	&BasePercent($genome,\%BaseContent);
}

if($DNA_singleBase =~ /\S+\.gz$/) {
	open (IN,"gunzip -c $DNA_singleBase |") or die $!;
} else {
	open (IN,$DNA_singleBase) or die $!;
}

my %homo;
while (<IN>) {
	my $line1 = $_;
	chomp $line1;
	my @info = split /\s+/,$line1;

	next unless (exists $RNAediting{$info[0]}{$info[1]});

	my ($ref_base,$m_bas,$m_qua,$C) = (split /\s+/,$RNAediting{$info[0]}{$info[1]});
	$ref_base=uc $ref_base;
	$info[2]=uc $info[2];
	die "Wrong! RNA ref_base no equal DNA ref_base! $ref_base ne $info[2] at $info[0]:$info[1]" unless ($ref_base eq $info[2]);

#	my $coverage = $info[5]; 
#	next unless $coverage >= $DNAdepth;
	my $coverage = 0;
	my ($seq,$qua) = &dealBasequailty(@info);
	$seq=uc $seq;
	if ($seq) {
		$coverage += length($seq);
	}else{
		next;
	}
#	next unless $coverage >= $DNAdepth;
###
	my %dna_hash;
	my @dna_array=split //,$seq;
	foreach my $dna_base (@dna_array){
		$dna_hash{$dna_base}++;
	}
	my @dna_info;
	if(!exists $dna_hash{'A'}){$dna_hash{'A'}=0;}
	if(!exists $dna_hash{'C'}){$dna_hash{'C'}=0;}
	if(!exists $dna_hash{'G'}){$dna_hash{'G'}=0;}
	if(!exists $dna_hash{'T'}){$dna_hash{'T'}=0;}
	push @dna_info,$dna_hash{'A'},$dna_hash{'C'},$dna_hash{'G'},$dna_hash{'T'};
	my $dna_info=join(",",@dna_info);
###
	my @result;
	if($method eq "Bayesian"){
		@result = &SNPvalue_Bayesian($seq,$qua, $ploidy,\%BaseContent);
		$result[1]=sprintf "%.4f",$result[1];
		$homo{$info[0]}{$info[1]} = [$coverage,$ref_base,$dna_info,$result[0],$result[1]];
	}elsif($method eq "Binomial"){
		@result = &SNPvalue_Binomial($dna_info,$ref_base,$ploidy);
		$homo{$info[0]}{$info[1]} = [$coverage,$ref_base,$dna_info,$result[0]];
	}elsif($method eq "Stringent"){
		@result = &SNPvalue_Stringent($seq,$ref_base);
		$homo{$info[0]}{$info[1]} = [$coverage,$ref_base,$dna_info,$result[0],$result[1]];
	}else{
		die "Error: unknown --Method $method";
	}

###
}
close IN;

if($method eq "Binomial"){
	my @Barray;
	my @Parray;
	foreach my $chr (keys %homo){
		foreach my $site (keys %{$homo{$chr}}){
			push @Barray,[@{$homo{$chr}{$site}},$chr,$site];
			push @Parray,$homo{$chr}{$site}->[3];
		}
	}
	@Barray=sort{$b->[3]<=>$a->[3]} @Barray;
	@Parray=sort{$b<=>$a} @Parray;
	my @Parray_cor=&Fdr(\@Parray);
	for(my $i=0;$i<@Barray;$i++){
		$Barray[$i][3] = &formatP($Barray[$i][3]);
		$Parray_cor[$i] = $Parray_cor[$i]>1 ? &formatP(1) : &formatP($Parray_cor[$i]);
		$homo{$Barray[$i][4]}{$Barray[$i][5]}= [$Barray[$i][0],$Barray[$i][1],$Barray[$i][2],$Barray[$i][3],$Parray_cor[$i]]
	}
}

my (%base_complementrity) = (
		'T' =>'A',
		'A' =>'T',
		'C' =>'G',
		'G' =>'C',
		);
my %site;

$time= &getTime();
print STDERR "$time\tProccessing detail information of candidate RNA editing sites ...\n";

foreach my $chr (sort keys %RNAediting) {
	for my $pos (sort {$a <=> $b} keys %{$RNAediting{$chr}}) {
		if (exists $homo{$chr}{$pos}) {
			my @temp = split /\t/,$RNAediting{$chr}{$pos};
			my ($cov,$homobase,$dna_info,$homo_tag,$posterior) = @{$homo{$chr}{$pos}};
			my $line = join "\t",$chr,$pos,$homobase,$cov,@temp[1..$#temp];
			my @info = split /\s+/,$line;
			next if ($info[2] eq "N");
			die "These bases($info[0]\t$info[1]) doesn't match!\n" unless ($info[4] =~ /^\[M\]:([ATCGNatcgn]+)$/);
			my @bases = split //,$1;
			die "These qualities($info[0]\t$info[1]) doesn't match!\n" unless ($info[5] =~ /^\[M\]:(\S+)/);
			my @qualities = split //,$1;
			die "Something was wrong at the location[$info[0]\t$info[1]]!\n" unless @bases == @qualities;
			my %editHash;
			my $ref_count = 0;

			my %rna_hash;
			foreach my $base (@bases) {
				$base = uc($base);
				if ($base eq "$info[2]") {
					$ref_count += 1;
				} else {
					$editHash{$base} += 1;
				}
				$rna_hash{$base}++;
			}
			if(!exists $rna_hash{'A'}){$rna_hash{'A'}=0;}
			if(!exists $rna_hash{'C'}){$rna_hash{'C'}=0;}
			if(!exists $rna_hash{'G'}){$rna_hash{'G'}=0;}
			if(!exists $rna_hash{'T'}){$rna_hash{'T'}=0;}
			my $rna_info="$rna_hash{'A'},$rna_hash{'C'},$rna_hash{'G'},$rna_hash{'T'}";

			my (@best,@tag,$RNA_degree);
			my $mark = 0;
			foreach my $K (sort { $editHash{$b} <=> $editHash{$a} } keys %editHash) {
				my $stat = join ":",$K,$editHash{$K};
				if ($mark == 0) {
					push @best,$K,$editHash{$K};
					$RNA_degree = sprintf "%.4f",$editHash{$K}/$info[6];
				}
				push @tag,$stat;
				$mark ++;
			}

			if(@tag>2){             # filter multiple editing type
				my @num_array;
				foreach my $tag_ele (@tag){
					my ($num)=(split /:/,$tag_ele)[1];
					push @num_array,$num;
				}
				@num_array=sort{$b<=>$a} @num_array;
				next if $num_array[1]/$num_array[0]>0.01;
			}

			my $type;
			if ($strand eq "+") {
				$type = join "->",$info[2],$best[0];
				my $a=join ("",sort($info[2],$best[0]));
			} elsif ($strand eq "-") {
				my $ref = $base_complementrity{$info[2]};
				my $con = $base_complementrity{$best[0]};
				$type = join "->",$ref,$con;
				my $a=join ("",sort($ref,$con));
			} elsif ($strand eq "unknown") {
				my $ref = $base_complementrity{$info[2]};
				my $con = $base_complementrity{$best[0]};
				$type = "${info[2]}->${best[0]}|${ref}->${con}";
			}
			my $outline;
#			next if ($best[1] < $editDepth);
			next if ($best[1] < 3);
			if($strand eq "+" || $strand eq '-'){
				$outline = join "\t",@info[0,1],$strand,@info[2,3],$dna_info,$homo_tag,$posterior,$type,$info[6],$rna_info,$best[0],$best[1],$RNA_degree;
			}else{
				$outline = join "\t",@info[0,1],".",@info[2,3],$dna_info,$homo_tag,$posterior,$type,$info[6],$rna_info,$best[0],$best[1],$RNA_degree;
			}
			$site{$info[0]}{$info[1]}{info}=$outline;
		}
	}
}
undef %RNAediting;

my %site_sort;
my %log_total;
foreach my $chr (keys %site){
	my @array=sort {$a<=>$b} keys %{$site{$chr}};
	@{$site_sort{$chr}}=@array;
	my $n=@array;
	$log_total{$chr}=$n;
}


$time= &getTime();
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
	my $beg=$A[3];
	my $end= &offset($A[3],$A[5]);
	if(exists $site{$A[2]}){
		my @array= @{$site_sort{$A[2]}};
		my @next_array=@array;
		my $len=length($A[9]);
		my $tailEnd_len=int $len/4;
		foreach my $pos (@array){
			if($pos<$beg){
				shift @next_array;
			}
			if($pos>$end || $array[-1]<$beg){last;}
			if($pos>=$beg && $pos<=$end){
				my @info_array=split /\t/,$site{$A[2]}{$pos}{info};
				my $editType=$info_array[11];
				$editType=uc $editType;
				my $posOnRead=&locateOnRead($pos-$beg+1,$A[5]);
				if($posOnRead==0){next;}
				my $readBase=uc (substr($A[9], $posOnRead-1, 1));
				my $readQual=substr($A[10],$posOnRead-1,1);
				if($readBase ne $editType){next;}
				my $q=ord($readQual)-$phred;
				if($q<$qual_cutoff){next;}
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
		@{$site_sort{$A[2]}}=@next_array;
		my $n=@array;
		my $next_n=@next_array;
		my $ratio=sprintf "%.2f",100-$next_n/$log_total{$A[2]}*100;
		if($next_n<$n){
#			print STDERR "\t$A[2]: ${ratio}% sites finish. (Total: $log_total{$A[2]})\n";
			print STDERR "\t$A[2]: ${ratio}% sites finish.\n";
		}
	}
}
close IN;

my @readInfo_result;
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
		my @line=split /\t/,"$site{$scaf}{$pos}{info}\t$site{$scaf}{$pos}{end}\t$site{$scaf}{$pos}{mid}\t$site{$scaf}{$pos}{readNum}";
		push @readInfo_result,[@line];
	}
}
undef %site_sort;
undef %site;
my @Parray;

$time= &getTime();
print STDERR "$time\tCalculating P-value of each RNA editing site ...\n";


foreach my $value (@readInfo_result){
	my $line=join("\t",@$value);
	my $editDep=$value->[12];
	my $totalDep=$value->[9];

	my $errorRate=10**(-1*$qual_cutoff/10);       # Q=-10log(10)P

	my $upper_tail_p = &pbinom($editDep, $totalDep, $errorRate, "upper");
	push @Parray,[($line,$upper_tail_p)];

}
@readInfo_result=();
##
@Parray=sort {$b->[1]<=>$a->[1]} @Parray;
my @p;
foreach my $a (@Parray){
	push @p,$a->[1];
}
my @p_cor = Fdr(\@p);

my @pvalue_result;
for(my $i=0;$i<@Parray;$i++){
	$p_cor[$i]= &formatP($p_cor[$i]);
	$Parray[$i][1]= &formatP($Parray[$i][1]);
	my @A=split /\t/,$Parray[$i][0];
	push @pvalue_result, [@A,$Parray[$i][1],$p_cor[$i]];
}

@pvalue_result= sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @pvalue_result;

my ($title1,$title2);
if($method eq "Bayesian"){
	($title1,$title2)=("Bayesian_Genotype(DNA)","Bayesian_Posterior_Probability(DNA)");
}elsif($method eq "Binomial"){
	($title1,$title2)=("P_value(DNA_Heterozygosis)","FDR(DNA_Heterozygosis)");
}elsif($method eq "Stringent"){
	($title1,$title2)=("Non_Ref_BaseCount(DNA)","Non_Ref_BaseRatio(DNA)");
}else{
	die "Error: unknown --Method $method\n";
}

print "1.Chromosome_ID\t2.Coordinate\t3.strand\t4.Ref_base\t5.coverage(DNA)\t6.DNA_BaseCount[A,C,G,T]\t7.$title1\t8.$title2\t9.coverage(RNA)\t10.RNA_BaseCount[A,C,G,T]\t11.RNA_Editing_Type\t12.RNA_Editing_Level\t13.Quadrant1/4_Read_Num(RNA)\t14.Quadrant2/3_Read_Num(RNA)\t15.Read_Type_Num(RNA)\t16.P_value(RNA_editing)\t17.FDR(RNA_editing)\n";

foreach my $value (@pvalue_result){
	my ($Chromosome_ID,$Coordinate,$strand,$gbase,$coverage_DNA,$DNA_BaseCount,$label1,$label2,$editType,$coverage_RNA,$RNA_BaseCount,$BestEditBase_RNA,$coverage_BestEditBase,$editing_degree,$badRead,$goodRead,$ReadTypeNum,$Pvalue,$FDR)=@$value;
	my $line=join("\t",$Chromosome_ID,$Coordinate,$strand,$gbase,$coverage_DNA,$DNA_BaseCount,$label1,$label2,$coverage_RNA,$RNA_BaseCount,$editType,$editing_degree,$badRead,$goodRead,$ReadTypeNum,$Pvalue,$FDR);
	print "$line\n";
}

$time= &getTime();
print STDERR "$time\tRaw candidate RNA editing sites is generated.!\n";

##

#####################################	subroutine	##############################################
sub dealBasequailty {
	my @array = @_;
	die "These bases($array[0]\t$array[1]) doesn't match!\n" unless ($array[3] =~ /^\[M\]:([ATCGNatcgn0]+)(;\[I\]:[ATCGNatcgn,]+)?(;\[D\]:[ATCGNatcgn,]+)?/);
	my @bases = split //,$1;
	my $number = $#bases+1;
	my @qualities;
	my @aa = split //,$array[4];
	for (my $k=4;$k<$number+4;$k++) {
		push @qualities,$aa[$k];
	}
	die "Something was wrong at the location[$array[0]\t$array[1]]!\n" unless @bases == @qualities;
	my (@b_temp,@q_temp);
	for (my $i=0;$i<@bases;$i++) {
		my $score = ord($qualities[$i]) - $phred;
		if ($score >= $qual_cutoff) {
			push @b_temp,$bases[$i];
			push @q_temp,$qualities[$i];
		}
	}
	my $b_temp = join "",@b_temp;
	my $q_temp = join "",@q_temp;
	return ($b_temp,$q_temp);
}


sub QualityProbability{
	my$score=shift;
	my$value=ord($score)-$phred;
	my$p=10**(0-$value/10);
#	$p=$p/(1+$p);  # for solexa
	return($p);
}


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
			$M[$i]=~s/S//;
			if($i>0){$offset=$offset+$M[$i];}
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

sub Fdr{
	my $p=shift;
	my $n=@$p;
	my $temp;
	my $i=0;  my $pre;
	my @Fdr_result;
	foreach my $pvalue (@$p){
		if($i==0){$temp=$n;}elsif($pvalue<$pre){$temp=$n-$i}
		$pre=$pvalue;
		my $cor=$pvalue*$n/$temp;
		push @Fdr_result,$cor;
		$i++;
	}
	return @Fdr_result;
}

#################################################################################################
## Calculate the p for a given number of successes in a Bernoulli experiment.
## Example:
## my $p = dbinom(2, 10, 0.3);
sub dbinom {
	my ($x, $n, $pie) = @_; # $x: number of successes; $n: number of trials; $pie: probability of success on each trial.
		die if $x > $n || $x < 0 || $n < 0;
	die unless $pie > 0 && $pie < 1;
## calculate the factorial part of the binomial formula, use logarithmic operation to save computing time.
	my ($a, $b);
	if ($x > $n/2) {
		$a = $x;
		$b = $n - $x;
	} else {
		$a = $n - $x;
		$b = $x;
	}
	my $log_fac1 = 0;
	for (my $i = $n; $i > $a; $i --) {
		$log_fac1 += log($i);
	}
	my $log_fac2 = 0;
	for (my $i = $b; $i >= 1; $i --) {
		$log_fac2 += log($i);
	}
	my $log_fac = $log_fac1 - $log_fac2;

## calculate the p for a given number of successes.
	my $log_p = $log_fac + $x * log($pie) + ($n - $x) * log(1 - $pie);
	my $p = exp($log_p);
	return $p;
}


sub pbinom {
## $x: number of successes; $n: number of trials; $pie: probability of success on each trial.
## $lower_tail: if TRUE, probabilities are P[X <= x], otherwise, P[X > x].
	my ($x, $n, $pie, $tail) = @_;
	my $tail_p = 0;
	if ($tail eq "lower") {
		if ($x <= $n/2) {
			for (my $i = 0; $i <= $x; $i ++) {
				$tail_p += dbinom($i, $n, $pie);
			}
		} else {
			for (my $i = $x + 1; $i <= $n; $i ++) {
				$tail_p += dbinom($i, $n, $pie);
			}
			$tail_p = 1 - $tail_p;
		}
	} elsif ($tail eq "upper") {
		if ($x <= $n/2) {
			for (my $i = 0; $i <= $x - 1; $i ++) {
				$tail_p += dbinom($i, $n, $pie);
			}
			$tail_p = 1 - $tail_p;
		} else {
			for (my $i = $x; $i <= $n; $i ++) {
				$tail_p += dbinom($i, $n, $pie);
			}
		}
	} else {
		die "tail is not set correctly, set it to lower or upper\n";
	}
	$tail_p = dbinom($x, $n, $pie) if $tail_p <= 0;
	return $tail_p;
}



sub formatP {
	my ($p) = @_;
	if ($p > 0.0001) {
		$p = sprintf "%.8f", $p;
	} else {
		$p = sprintf "%.8e", $p;
	}
	return $p;
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


##

sub SNPvalue_Binomial{
	my ($info,$r,$p)=@_;
	my $ratio=1/$p;	
	my @I=split /,/,$info;
	my %info_hash=("A",$I[0],"C",$I[1],"G",$I[2],"T",$I[3]);
	my $totalDep=$I[0]+$I[1]+$I[2]+$I[3];
	my $diffDep=0;
	foreach my $base (keys %info_hash){
		if($base ne $r){
			$diffDep+=$info_hash{$base};
		}
	}
	my $lower_tail_p = &pbinom($diffDep, $totalDep, $ratio, "lower");
	return ($lower_tail_p);
}


sub SNPvalue_Stringent{
	my ($s,$r)=@_;
	my @S=split //,$s;
	my $totalDep=@S;
	my $diffDep=0;
	foreach my $base (@S){
		if($base ne $r){
			$diffDep++;
		}
	}
	my $ratio=sprintf "%.4f",$diffDep/$totalDep;
	return ($diffDep,$ratio);
}
##

sub SNPvalue_Bayesian{
#	my $HomoPrior = 0.99;
#	my $HetePrior = 0.01;
#	my $rate = 4; #the transition-transversion rate ratio
	my $weight = 0.5; #the weight of A/C and G/T;
	my $weight_other = (1-$weight)/2; #the weight of the other substitutions 
	my @pileup = split //,$_[0];
	my @quals = split //,$_[1]; 
	my $ploid = $_[2];
	my $BaseContent=$_[3];
	my @Type = ('A','T','G','C');
	my %FixedError = (
		'A'=>{
			'C'=>$weight,
			'T'=>$weight_other,
			'G'=>$weight_other,
		},
		'C'=>{
			'A'=>$weight,
			'T'=>$weight_other,
			'G'=>$weight_other,
		},
		'G'=>{
			'T'=>$weight,
			'A'=>$weight_other,
			'C'=>$weight_other,
		},
		'T'=>{
			'G'=>$weight,
			'A'=>$weight_other,
			'C'=>$weight_other,
		},
	);

	my %P_allele;	#The probability of a given allele
	foreach my $i (0..$#pileup) {
		next unless ($pileup[$i] =~ /[ATGC]/);
		my $score = &QualityProbability($quals[$i]);
		foreach my $base (@Type) {
			if ($base eq $pileup[$i]) { 
				push @{$P_allele{$base}}, 1-$score
			} else {
				push @{$P_allele{$base}}, $score * $FixedError{$base}{$pileup[$i]}};
		}
	}

	my (%posterior, %prior, %TypePrior);
	foreach my $base (@Type) {
		my $genotype = $base x $ploid;
		my $value=1;
		foreach my $score (@{$P_allele{$base}}) {
			$value += log($score);
		}
		$prior{$genotype} = $value;	#The probability of pileup data given the genotype (homozygous): p(D|G);
		$TypePrior{$genotype} = $HomoPrior * $$BaseContent{$base};	#the prior probability of each genotype (homozygous): p(G);
	}
	my $HeteNum = ($ploid-1)*4+($ploid-1)*2*$rate;	
	my %control;
	foreach my $key (keys %FixedError){
		foreach my $ke (keys %{$FixedError{$key}}){                                # $key: original base; $ke: observed base;
			next if (exists $control{$key}{$ke});
			next if ($ploid == 1);
			$control{$key}{$ke} ++;
			$control{$ke}{$key} ++;
			for(my $i=1;$i<$ploid;$i++){
				my ($genotype,$value) = ("",1);
				$genotype .= $key x $i;
				$genotype .= $ke x ($ploid-$i); 
				
				for (my $j=0;$j<=$#{$P_allele{$key}};$j++){
					$value += log( ($i*$P_allele{$key}[$j] + ($ploid-$i)*$P_allele{$ke}[$j]) / $ploid);
				}
				$prior{$genotype} = $value;	#The probability of pileup data given the genotype (heterozygous): p(D|G);
				
				if (($key eq "G" && $ke eq "A") || ($key eq "A" && $ke eq "G") || ($key eq "T" && $ke eq "C") || ($key eq "C" && $ke eq "T")) {
					$TypePrior{$genotype} = $HetePrior/$HeteNum*$rate; #the prior probability of each genotype (heterozygous): p(G) [genotypes including the allele pair of A/G or T/C];
				} else {
					$TypePrior{$genotype} = $HetePrior/$HeteNum; #the prior probability of each genotype (heterozygous): p(G);
				}
			}
		}
	}

	my $P_pileup = 0;
	my @greatest = sort {$prior{$b}<=>$prior{$a}} keys %prior;
	my $max = $prior{$greatest[0]};
	foreach my $genotype (keys %prior){
		$prior{$genotype} -= $max;
		my $P_genotype;	#the prior probability of corresponding genotype;
		if (exists $TypePrior{$genotype}){
			$P_genotype = $TypePrior{$genotype};
		} else {
			die "There is not the prior probability for $genotype\n";
		}
		$P_pileup += $P_genotype * exp($prior{$genotype});	#The probability of pileup data: p(D)
	}
	
	foreach my $genotype (keys %prior){
		my $P_genotype;
		if (exists $TypePrior{$genotype}){
			$P_genotype = $TypePrior{$genotype};
		} else {
			die "There is not the prior probability for $genotype\n";
		}
		$posterior{$genotype} = $P_genotype * exp($prior{$genotype}) / $P_pileup;
	}
	
	my @sort = sort{ $posterior{$b}<=>$posterior{$a} }keys %posterior;
	my $Type = &CheckHomo($sort[0]);

#	return($Type,$posterior{$sort[0]},$sort[0]);
	return($sort[0],$posterior{$sort[0]});
}

sub CheckHomo{
	my $string = $_[0];
	my $base = (split //, $string)[0];
	$string =~ s/$base//g;
	my $zygous;
	if ($string eq "") {
		$zygous = "HOMO";
	} else {
		$zygous = "HETE";
	}
	return $zygous;
}


sub BasePercent {
	my ($file,$ref) = @_;
	if ($file =~ /\.gz$/) {
		open (INN,"gunzip -c $file |") or die $!;
	} else {
		open (INN,$file) or die $!;
	}
	$/ = ">";
	<INN>;
	my ($total,$A_num,$T_num,$C_num,$G_num);
	while (<INN>) {
		chomp;
		my $scaf;
		if (/(\S+)/) {
			$scaf = $1;
		}
		s/.+\n//;
		s/\s+//g;
		$_ =~ s/N//ig;
		$total += length($_);
		$A_num += $_ =~ s/A//ig;
		$T_num += $_ =~ s/T//ig;
		$C_num += $_ =~ s/C//ig;
		$G_num += $_ =~ s/G//ig;
	}
	close INN;
	$/ = "\n";
	my $A_ratio = sprintf "%.4f",$A_num/$total;
	my $T_ratio = sprintf "%.4f",$T_num/$total;
	my $C_ratio = sprintf "%.4f",$C_num/$total;
	my $G_ratio = sprintf "%.4f",$G_num/$total;
	$ref->{"A"} = $A_ratio;
	$ref->{"T"} = $T_ratio;
	$ref->{"C"} = $C_ratio;
	$ref->{"G"} = $G_ratio;
}
