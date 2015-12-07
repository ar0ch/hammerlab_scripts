#!/usr/bin/perl -w
# Aroon Chande
# Extracted HMM profiled proteins from a proteins.faa file
# Input os output of  hmmscan --tblout
# ./protein_extractor.pl -hmm hmmdir -prot proteindir -o outdir
# For hmmscan batch jobs:
# for i in $(ls *.hmm);do; for j in $(ls ./Proteins/*.faa);do;\
# k=$(echo $j| cut -c 12-);l=$(echo $i| rev| cut -c 5- | rev);\
# echo $k;echo $l_out;hmmscan --tblout ./$l.out/$k.txt -o \
# $l.out/$l_$k.out $i $j;done;done

use strict;
use Getopt::Long;
use File::Basename;
my $prog = basename($0);
sub print_usage
{
    warn <<"EOF";

USAGE
  $prog -hmm <hmmdir> -prot <protdir> -o <outdir>

DESCRIPTION
  This program extracts the hits from an `hmmscan` file
  hmmscan must be run with `--tblout`. Other hmm* profiling
  output is also accepted. 

OPTIONS
  -h            Print this help message
  -hmm	dir		Directory containing HMMscan output in tab delim
				format
  -o    dir     Dumps extracted seqs to specifed directory
  -prot	dir		Directory containing protein files in .faa format
				Protein and hmm file must have same basename


EXAMPLES
  $prog -hmm hmmscan_vgrg -prot reference_prots -o hits_out
  $prog -h

EXIT STATUS
  0     Successful completion
  >0    An error occurred

EOF
}
my $h=0;
my ($hmmDir,$protDir,$outdir,$hmmFile,$protFile,$outfile,$prots,@prots,%prots,@hits);
$prots = '';
if (@ARGV < 1){print_usage();exit 1;}
GetOptions ('hmm=s' => \$hmmDir, 'prot=s' => \$protDir, 'o=s' => \$outdir, 'h' => \$h );
if (eval $h){ print_usage();exit 1;}
opendir DIR, $hmmDir or die "cannot open dir $hmmDir: $!";
my @file= readdir DIR;
closedir DIR;
foreach (@file){
	if ($_ =~ /.txt/){
		$hmmFile = join('/',$hmmDir,$_);
		$_ =~ s/\.txt//g;
		$protFile = join('/',$protDir,$_);
		$outfile = join('/',$outdir,$_);
		#Load protein file
		open PROT, $protFile or die "Cannot open $protFile: $!\n";
		foreach (<PROT>){
				$prots .= $_;
		}
		close PROT;
		@prots = split(/\>/,$prots);
		shift @prots;
		foreach (@prots){
			my($first, $last) = split(/\r?\n/,$_,2);
			$prots{$first} = $last;
		}

		open HMM, $hmmFile or die "Cannot open $hmmFile: $!\n";
		my $i = 0;
		foreach (<HMM>){
			if(($_ =~ /\#/)){next;}
			else{
				my @columns = split(/\s+/,$_);
				$hits[$i] = $columns[2];
				$i++;
				}
		}
		close HMM;
		open OUT, ">$outfile" or die "Cannot open $outfile: $!\n";
		foreach (@hits){
			print OUT ">$_\n$prots{$_}\n"
		}
	}
	else {next;}
}
