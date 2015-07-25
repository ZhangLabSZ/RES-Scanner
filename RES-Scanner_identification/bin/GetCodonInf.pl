#!/usr/bin/perl -w
use strict;
die "Usage: perl $0 <CDS.pos> <genome.fa>\n" unless @ARGV == 2;

my %CODE = (	
				'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',                               # Alanine
				'TGC' => 'C', 'TGT' => 'C',                                                           # Cysteine
				'GAC' => 'D', 'GAT' => 'D',                                                           # Aspartic Acid
				'GAA' => 'E', 'GAG' => 'E',                                                           # Glutamic Acid
				'TTC' => 'F', 'TTT' => 'F',                                                           # Phenylalanine
				'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                               # Glycine
				'CAC' => 'H', 'CAT' => 'H',                                                           # Histidine
				'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                             # Isoleucine
				'AAA' => 'K', 'AAG' => 'K',                                                           # Lysine
				'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',   # Leucine
				'ATG' => 'M',                                                                         # Methionine
				'AAC' => 'N', 'AAT' => 'N',                                                           # Asparagine
				'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                               # Proline
				'CAA' => 'Q', 'CAG' => 'Q',                                                           # Glutamine
				'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',   # Arginine
				'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',   # Serine
				'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                               # Threonine
				'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                               # Valine
				'TGG' => 'W',                                                                         # Tryptophan
				'TAC' => 'Y', 'TAT' => 'Y',                                                           # Tyrosine
				'TAA' => 'U', 'TAG' => 'U', 'TGA' => 'U'                                              # Stop
);

my %CDS;

#######################
my (%gene,%ID);
dealPOS($ARGV[0],\%gene);
######################
if ($ARGV[1] =~ /\.gz$/) {
	open (IN,"zcat $ARGV[1] |") or die $!;
} else {
	open (IN,$ARGV[1]) or die $!;
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
	my $seq = $_;
	my $length = length($_);
	my $chr_pp=$gene{$scaf};
	foreach  my $gene (sort keys %$chr_pp) {
		my $strand=$$chr_pp{$gene}{strand};
		my @exon = @{$chr_pp->{$gene}{exon}};
		my $Temp_seq;
		for (my $i=0; $i<@exon; $i++) {
			my $exon = substr($seq,$exon[$i][0]-1, $exon[$i][1] - $exon[$i][0] + 1);
			$exon = Complement_Reverse($exon) if($strand eq '-');
			$Temp_seq .= uc($exon);
		}
		next unless (length($Temp_seq)%3 == 0);
#		print STDERR "$gene\n" unless (length($Temp_seq)%3 == 0);
		$ID{$gene} = 1;
		my @bases = split //,$Temp_seq;
		my $Codon = join "",@bases[-3,-2,-1];
		$CODE{$Codon} = "*" unless (exists $CODE{$Codon});
		if ($CODE{$Codon} eq "U") {
		} else {
			my $stop = substr($seq,$$chr_pp{$gene}{codon}->[0] - 1,3);
			$stop = Complement_Reverse($stop) if($strand eq '-');
			$CODE{$stop} = "*" unless (exists $CODE{$stop});
			if ($CODE{$stop} eq "U") {
				$Temp_seq .= $stop;
			}
		}
#		Display_seq(\$Temp_seq);
#		print "$Temp_seq";
	
		$CDS{$gene} = $Temp_seq;
	}
}
close IN;
$/ = "\n";


my %pos;
open (IN,$ARGV[0]) or die $!;
while (<IN>) {
	next if (/^#/);
	chomp;
	my @info = split /\s+/;
	my @A = (split /_/,$info[0]);
	pop @A;
	my $id=join("_",@A);
	push @{$pos{$id}{$info[2]}},[$info[3],$info[4],$info[1]];
}
close IN;

foreach my $id (sort keys %pos) {
	foreach my $strand (sort keys %{$pos{$id}}) {
		my %p_to_bg;
		foreach my $p (@{$pos{$id}{$strand}}) {
			$p_to_bg{$p} = $p->[0];
		}
		@{$pos{$id}{$strand}} = sort {$p_to_bg{$a} <=> $p_to_bg{$b}} @{$pos{$id}{$strand}} if $strand eq "+";
		@{$pos{$id}{$strand}} = sort {$p_to_bg{$b} <=> $p_to_bg{$a}} @{$pos{$id}{$strand}} if $strand eq "-";
	}
}

my %locus;
foreach my $id (sort keys %ID) {
	foreach my $strand (sort keys %{$pos{$id}}) {
		my $k = 0;
		foreach my $p (@{$pos{$id}{$strand}}) {
			my ($aa,$bb,$sca) = @{$p};
			if ($strand eq "+") {
				for (my $i=$aa;$i<=$bb;$i++) {
					$k ++;
					my $num = int(($k-1)/3);
					my $re = ($k%3);
					my $str_a = 3*$num;
					my $str_b = 3*($num+1)-1;
					my $codon = substr($CDS{$id},$str_a,3);
					my $phase;
					$phase = 0 if ($re == 1);
					$phase = 1 if ($re == 2);
					$phase = 2 if ($re == 0);
					my $base = substr($CDS{$id},$str_a+$phase,1);
					#$locus{$sca}{$i}{$strand} = [$base,$phase,$codon];
					$CODE{$codon} = "*" unless (exists $CODE{$codon});
					my $lin = join "\t",$sca,$i,$strand,$base,$phase,$codon,$CODE{$codon},$num+1;
					print "$lin\n"
				}
			}  
			if ($strand eq "-") {
					for (my $i=$bb;$i>=$aa;$i--) {
					$k ++;
					my $num = int(($k-1)/3);
					my $re = ($k%3);
					my $str_a = 3*$num;
					my $str_b = 3*($num+1)-1;
					my $codon = substr($CDS{$id},$str_a,3);
					my $phase;
					$phase = 0 if ($re == 1);
					$phase = 1 if ($re == 2);
					$phase = 2 if ($re == 0);
					my $base = substr($CDS{$id},$str_a+$phase,1);
					#$locus{$sca}{$i}{$strand} = [$base,$phase,$codon];
					$CODE{$codon} = "*" unless (exists $CODE{$codon});
					my $lin = join "\t",$sca,$i,$strand,$base,$phase,$codon,$CODE{$codon},$num+1;
					print "$lin\n"
				}

			}
		}
	}
}

########################	Subroutine	#######################
sub dealPOS {
	my ($file,$ref) = @_;
	open (IN,$file) or die $!;
	while (<IN>) {
		next if (/^#/);
		my @info = split /\s+/;
		my @A = (split /_/,$info[0]);
		pop @A;
		my $id=join("_",@A);
		$ref->{$info[1]}{$id}{strand} = $info[2];
		push @{$ref->{$info[1]}{$id}{exon}},[$info[3],$info[4]];
	}
	close IN;
	foreach my $chr (keys %$ref) {
		my $chr_p = $ref->{$chr};
		foreach my $gene (keys %$chr_p) {
			my $gene_p = $chr_p->{$gene};
			my @exon;
			@exon = sort {$a->[0] <=> $b->[0]} @{$gene_p->{exon}} if ($gene_p->{strand} eq "+");     # for stop codon
			$gene_p->{codon} = [$exon[-1]->[1]+1, $exon[-1]->[1]+3] if ($gene_p->{strand} eq "+");
			@exon = sort {$b->[0] <=> $a->[0]} @{$gene_p->{exon}} if ($gene_p->{strand} eq "-");
			$gene_p->{codon} = [$exon[-1]->[0]-3, $exon[-1]->[0]-1] if ($gene_p->{strand} eq "-");
			$gene_p->{exon} = \@exon;
		}
	}
}

sub Complement_Reverse{
	my $seq=shift;
	$seq=~tr/AGCTagct/TCGAtcga/;
	$seq=reverse($seq);
	return $seq;
}

sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;
	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
########################	END	#######################
