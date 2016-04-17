#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <rawSite.file1> <rawSite.file2> ... <editSite.file1> <editSite.file2> ...\n" unless @ARGV>=2;
my @array;
my $filenum=@ARGV;
die "Error: number of input files is not even!\n" if $filenum%2;

for(my $i=0;$i<$filenum/2;$i++){
	if($ARGV[$i]=~/\.gz$/){
		open IN,"gunzip -c $ARGV[$i] |" or die $!;
	}else{
		open IN,"$ARGV[$i]" or die $!;
	}
	while(<IN>){
		chomp;
		if($_=~/^1.Chromosome_ID/){
			next;
		}
		my @A=split /\t/;
		push @array,[$A[0],$A[1],$A[2],$A[12]];
	}
	close IN;
}


my $title;
my %hash;
for(my $i=$filenum/2;$i<$filenum;$i++){
	if($ARGV[$i]=~/\.gz$/){
		open IN,"gunzip -c $ARGV[$i] |" or die $!;
	}else{
		open IN,"$ARGV[$i]" or die $!;
	}
	while(<IN>){
		chomp;
		if($_=~/^1.Chromosome_ID/){
			$title=$_;
			next;
		}
		my @A=split /\t/;
		my $key=$A[0]."_".$A[1]."_".$A[2];
		$hash{$key}=$_;
	}
	close IN;
}


@array= sort {$b->[3] <=> $a->[3]} @array;
my @p;
foreach my $a (@array){
	push @p,$a->[3];
}
my @p_cor = Fdr(\@p);

for(my $i=0;$i<@array;$i++){
	my $key=$array[$i][0]."_".$array[$i][1]."_".$array[$i][2];
	if(exists $hash{$key}){
		my @A=split /\t/,$hash{$key};
		$p_cor[$i]= &formatP($p_cor[$i]);
		$A[13]=$p_cor[$i];
		$hash{$key}=join("\t",@A);
	}

}
print "$title\n";
foreach my $key (sort keys %hash){
	print "$hash{$key}\n";
}

#########################################
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



sub formatP {
	my ($p) = @_;
	if ($p > 0.0001) {
		$p = sprintf "%.8f", $p;
	} else {
		$p = sprintf "%.8e", $p;
	}
	return $p;
}






