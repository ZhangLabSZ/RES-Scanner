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
