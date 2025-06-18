#!/usr/bin/perl

# Xin Qiao, 25 Feb 2025
# Xin Qiao, 23 Mar 2025
# Xin Qiao, 27 Mar 2025
# Xin Qiao, 11 Jun 2025
# Xin Qiao, 13 Jun 2025

use warnings;
use strict;
use File::Basename;
use File::Spec;

#die "Usage: perl ",basename($0)," <lineage dataset> <species code>\n" if $#ARGV < 1;

use Bio::SeqIO;
use Bio::Seq;

use Getopt::Std;

# obtain current directory
my $path_curf = File::Spec->rel2abs(__FILE__);
#print "$path_curf\n";
my ($vol, $dir, $file) = File::Spec->splitpath($path_curf);
#print "$vol\ncurrent dir: $dir\n$file\n";

my %options=();
getopts("i:o:l:s:", \%options);
if(!exists $options{i} || !exists $options{o} || !exists $options{l} || !exists $options{s})
{
print "Usage: ulove.pl -i input_directory -o output_directory -l lineage_dataset -s species_code\n";
print "#####################\n";
print "-i the directory storing protein sequence file, please name the file like this Arath.pep (*species_code*.pep)\n";
print "-o the directory storing output data\n";
print "-l the lineages dataset used for completeness assessement, please choose an appropriate lineage for your species:
    viridiplantae
    chlorophyta
    streptophyta
    embryophyta
    tracheophyta
    spermatophyta
    angiosperms
    monocots
    eudicots\n";
print "-s species code, for example, Arath can be used as the species code of Arabidopsis thaliana\n";
exit;
}
$options{i}=~s/\/$//;
$options{o}=~s/\/$//;

my $hmm=$dir."hmm";
#system "mkdir $options{o}\/hmmsearch_results";
my $hmm_out="$options{o}\/hmmsearch_results";
unless(-d $hmm_out)
{
	system "mkdir $options{o}\/hmmsearch_results";
}
open IN, $dir."hmm/$options{l}/$options{l}_LC_OGs_list.txt" or die "Cannot open $options{l}_LC_OGs_list.txt!";
while(<IN>)
{
	chomp;
	system "hmmsearch -E 1e-5 --domtblout $hmm_out\/$_\_$options{s}\.domtblout -o $hmm_out\/$_\_$options{s}\.out $hmm\/$options{l}\/$_.hmm $options{i}\/$options{s}\.pep";
}
close IN;

my $viri=1163;
my $chlo=27;
my $stre=1136;
my $embr=1129;
my $trac=1104;
my $sper=1093;
my $angi=1081;
my $mono=155;
my $eudi=899;

# Loading the list of score_cutoff and length_cutoff for universal low-copy orthologs
my %score_cutoff;
open IN, $dir."hmm/$options{l}/cutoff/score_cutoff" or die "Cannot open score_cutoff!";
while(<IN>)
{
	chomp;
	my @a=split /\t/, $_;
	$score_cutoff{$a[0]}=$a[1];
}
close IN;

my %length_cutoff_mean;
my %length_cutoff_sigma;
open IN, $dir."hmm/$options{l}/cutoff/length_cutoff" or die "Cannot open length_cutoff!";
while(<IN>)
{
	chomp;
	my @a=split /\t/, $_;
	$length_cutoff_sigma{$a[0]}=$a[2];
	$length_cutoff_mean{$a[0]}=$a[3];
}
close IN;

# Determine the maximal alignment score for each of genes contained in hmmsearch results
my %max_score;
open IN, $dir."hmm/$options{l}/$options{l}_LC_OGs_list.txt" or die "Cannot open $options{l}_LC_OGs_list.txt!";
while(<IN>)
{
	chomp;
	my $og=$_;
	open INN, "$hmm_out\/$og\_$options{s}.domtblout" or die "Cannot open $og\_$options{s}.domtblout!";
	while(<INN>)
	{
		chomp;
		if($_ =~ /^\#/){next;}
		my @a=split /\s+/, $_;
		if(!exists $max_score{$a[0]})
		{
			$max_score{$a[0]}=$a[7];
		}
		else
		{
			if($max_score{$a[0]} > $a[7])
			{
				$max_score{$a[0]}=$a[7];
			}
		}
	}
	close INN;
}
close IN;

# Gene annotation completeness assessment
my $o=0;
my $s=0;
my $d=0;
my $f=0;
my $m=0;
open IN, $dir."hmm/$options{l}/$options{l}_LC_OGs_list.txt" or die "Cannot open $options{l}_LC_OGs_list.txt!";
while(<IN>)
{
	$o++;
	chomp;
	my $og=$_;
	my $i=0;
	my %hmm_ali;
	my %score;
	open INN, "$hmm_out\/$og\_$options{s}.domtblout" or die "Cannot open $og\_$options{s}.domtblout!";
	while(<INN>)
	{
		chomp;
		if($_ =~ /^\#/){next;}
		$i++;
		my @a=split /\s+/, $_;
		if($a[7] <= $score_cutoff{$og})
		{
			$i--;
			next;
		}
		my $ali=$a[16]-$a[15];
		if(!exists $hmm_ali{$a[0]})
		{
			$hmm_ali{$a[0]}=$ali;
		}
		else
		{
			$hmm_ali{$a[0]}=$hmm_ali{$a[0]}+$ali;
		}
		$score{$a[0]}=$a[7];
	}
	close INN;
	if($i == 0)
	{
		$m++;
		next;
	}
	else
	{
		my $j=0;
		my $k=0;
		my $l=0;
		my $target_num=0;
		my %complete;
		my %vlarge;
		my %fragment;
		foreach my $key (keys %hmm_ali)
		{
			$target_num++;
			my $zeta=($length_cutoff_mean{$og}-$hmm_ali{$key})/$length_cutoff_sigma{$og};
			if($zeta >= -2 && $zeta <= 2) #complete
			{
				$j++;
				#For any duplicate gene matches within the same rank for different BUSCOs, keep only the highest scoring gene match.
				if($score{$key} < $max_score{$key})
				{
					$j--;
				}
				else
				{
					$complete{$key}=$score{$key};
				}
			}
			elsif($zeta < -2) #vlarge
			{
				$k++;
				if($score{$key} < $max_score{$key})
				{
					$k--;
				}
				else
				{
					$vlarge{$key}=$score{$key};
				}
			}
			else #fragment
			{
				$l++;
				if($score{$key} < $max_score{$key})
				{
					$l--;
				}
				else
				{
					$fragment{$key}=$score{$key};
				}
			}
		}
		my $line=0;
		my $n85_percent_of_top_score;
		open INN, "$hmm_out\/$og\_$options{s}.domtblout" or die "Cannot open $og\_$options{s}.domtblout!";
		while(<INN>)
		{
			chomp;
			if($_ =~ /^\#/){next;}
			my @a=split /\s+/, $_;
			if($a[7] <= $score_cutoff{$og})
			{
				next;
			}
			if($a[7] < $max_score{$a[0]})
			{
				next;
			}
			$line++;
			if($line == 1)
			{
				if(exists $fragment{$a[0]})
				{
					$line--;
					next;
				}
				else
				{
					$n85_percent_of_top_score=$a[7]*0.85;
				}
			}
		}
		close INN;
		if($target_num == 1)
		{
			#Remove duplicate gene matches of lesser importance, i.e. keep the complete ones, then the very large ones and finally the fragments.
			if($j == 0)
			{
				if($k == 0)
				{
					$f++;
				}
				else
				{
					$s++;
				}
			}
			else
			{
				$s++;
			}
		}
		else
		{
			#Remove duplicate gene matches of lesser importance, i.e. keep the complete ones, then the very large ones and finally the fragments.
			if($j == 0)
			{
				if($k == 0)
				{
					$f++;
				}
				else
				{
					if($k == 1)
					{
						$s++;
					}
					else
					{
						#$d++;
						#print "type1: $og\n";
						my $n=0;
						foreach my $key (keys %vlarge)
						{
							#remove any gene matches that score less than 85% of the top gene match score for each BUSCO.
							if($vlarge{$key} >= $n85_percent_of_top_score)
							{
								$n++;
							}
						}
						if($n == 0)
						{
							#$f++;
							print "type1: $og\n";
						}
						elsif($n == 1)
						{
							$s++;
						}
						else
						{
							$d++;
						}
					}
				}
			}
			elsif($j == 1)
			{
				if($k == 0)
				{
					$s++;
				}
				else
				{
					#$d++;
					#print "type2: $og\n";
					$s++;
				}
			}
			else
			{
				#$d++;
				#print "type3: $og\n";
				my $n=0;
				foreach my $key (keys %complete)
				{
					#remove any gene matches that score less than 85% of the top gene match score for each BUSCO.
					if($complete{$key} >= $n85_percent_of_top_score)
					{
						$n++;
					}
				}
				if($n == 0)
				{
					#$f++;
					if($k == 0)
					{
						print "type2: $og\n";
					}
					elsif($k == 1)
					{
						$s++;
					}
					else
					{
						my $nn=0;
						foreach my $key (keys %vlarge)
						{
							#remove any gene matches that score less than 85% of the top gene match score for each BUSCO.
							if($vlarge{$key} >= $n85_percent_of_top_score)
							{
								$nn++;
							}
						}
						if($nn == 0)
						{
							#$f++;
							print "type3: $og\n";
						}
						elsif($nn == 1)
						{
							$s++;
						}
						else
						{
							$d++;
						}
					}
				}
				elsif($n == 1)
				{
					$s++;
				}
				else
				{
					$d++;
				}
			}
		}
	}
}
close IN;

my $c=$s+$d;
my $sp=($s/$o)*100;
my $spf = sprintf "%.1f", $sp;
my $dp=($d/$o)*100;
my $dpf = sprintf "%.1f", $dp;
my $cp=($c/$o)*100;
my $cpf = sprintf "%.1f", $cp;
my $fp=($f/$o)*100;
my $fpf = sprintf "%.1f", $fp;
my $mp=($m/$o)*100;
my $mpf = sprintf "%.1f", $mp;

open OUT, ">$options{o}\/short_summary.specific.$options{l}.$options{s}.ulove.txt";
print OUT "# ULOVE version is: 1.0.0\n";
if($options{l} eq "viridiplantae")
{
	print OUT "# The lineage dataset is: viridiplantae (Creation date: 2025-06-06, number of genomes: $viri, number of ULOVEs: $o)\n";
}
elsif($options{l} eq "chlorophyta")
{
	print OUT "# The lineage dataset is: chlorophyta (Creation date: 2025-06-06, number of genomes: $chlo, number of ULOVEs: $o)\n";
}
elsif($options{l} eq "streptophyta")
{
	print OUT "# The lineage dataset is: streptophyta (Creation date: 2025-06-06, number of genomes: $stre, number of ULOVEs: $o)\n";
}
elsif($options{l} eq "embryophyta")
{
	print OUT "# The lineage dataset is: embryophyta (Creation date: 2025-06-06, number of genomes: $embr, number of ULOVEs: $o)\n";
}
elsif($options{l} eq "tracheophyta")
{
	print OUT "# The lineage dataset is: tracheophyta (Creation date: 2025-06-06, number of genomes: $trac, number of ULOVEs: $o)\n";
}
elsif($options{l} eq "spermatophyta")
{
	print OUT "# The lineage dataset is: spermatophyta (Creation date: 2025-06-06, number of genomes: $sper, number of ULOVEs: $o)\n";
}
elsif($options{l} eq "angiosperms")
{
	print OUT "# The lineage dataset is: angiosperms (Creation date: 2025-06-06, number of genomes: $angi, number of ULOVEs: $o)\n";
}
elsif($options{l} eq "monocots")
{
	print OUT "# The lineage dataset is: monocots (Creation date: 2025-06-06, number of genomes: $mono, number of ULOVEs: $o)\n";
}
elsif($options{l} eq "eudicots")
{
	print OUT "# The lineage dataset is: eudicots (Creation date: 2025-06-06, number of genomes: $eudi, number of ULOVEs: $o)\n";
}
else
{
	die "Please check the lineage name!\n";
}
print OUT "# Summarized benchmarking in ULOVE notation for file $options{i}\/$options{s}.pep\n";
print OUT "# ULOVE was run in mode: protein\n\n";
print OUT "\t***** Results: *****\n\n";
print OUT "\tC:$cpf%[S:$spf%,D:$dpf%],F:$fpf%,M:$mpf%,n:$o\n";
print OUT "\t$c\tComplete ULOVEs (C)\n";
print OUT "\t$s\tComplete and single-copy ULOVEs (S)\n";
print OUT "\t$d\tComplete and duplicated ULOVEs (D)\n";
print OUT "\t$f\tFragmented BUSCOs (F)\n";
print OUT "\t$m\tMissing ULOVEs (M)\n";
print OUT "\t$o\tTotal ULOVE groups searched\n";
print OUT "\nDependencies and versions:\n";
print OUT "\thmmsearch: 3.3.2\n";
print OUT "\tulove: 1.0.0\n";
close OUT;

# generate a figure to present ULOVE assessment result
system "python $dir\/generate_figure.py $options{s} $options{o}\/short_summary.specific.$options{l}.$options{s}.ulove.txt $options{o}";

__END__