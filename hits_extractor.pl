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
use File::Temp;
use Parallel::ForkManager;
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
  -h    		Print this help message
  -hmm	dir		Directory with tab delimited HMMscan output
  -o    di		 Dumps extracted seqs to specifed directory
  -prot	dir		Directory containing protein files in .faa format
			Protein and hmm file must have same basename
  -t	int		Number of threads for rps-bpast stage. Default: 10
  -type	str		Type of effector to verify. Default: all, types:
			all,hydrolase,lipase,lysm,ntpase,transferase,unknown


EXAMPLES
  $prog -hmm hmmscan_vgrg -prot reference_prots -o hits_out
  $prog -h

EXIT STATUS
  0     Successful completion
  >0    An error occurred

EOF
}
my $h=0;
my $threads=10;
my $type='all';
my ($hmmDir,$protDir,$outdir,$hmmFile,$protFile,$outfile,$prots,@prots,%prots,@hits);
$prots = '';
if (@ARGV < 1){print_usage();exit 1;}
GetOptions ('hmm=s' => \$hmmDir, 'prot=s' => \$protDir, 'o=s' => \$outdir, 'h' => \$h, 't=i' => \$threads,'type=s' => \$type);
if (eval $h){ print_usage();exit 1;}
my $manager = Parallel::ForkManager -> new ( $threads );
opendir DIR, $hmmDir or die "cannot open dir $hmmDir: $!";
my @file= readdir DIR;
closedir DIR;
foreach (@file){
	@hits = ( );
	$prots = '';
	%prots = (	);
	if ($_ =~ /.txt$/){
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
			my($desc, $seq) = split(/\r?\n/,$_,2);
			$prots{$desc} = $seq;
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
		foreach my $hit (@hits){
			$manager->start and next;
			my $rpstemp=temp_filename();
			my $rpsresulttemp=temp_filename();
			open RPS, ">$rpstemp" or die "Cannot open $rpstemp: $!";
			print RPS ">$hit\n$prots{$hit}\n";
			close RPS;
			`rpsblast -db /home/achande3/bin/db/cdd/Cdd -max_target_seqs 1 -outfmt "6" -out $rpsresulttemp -query $rpstemp`;
			open RPSR, "$rpsresulttemp" or die "Cannot open $rpsresulttemp: $!";
			foreach (<RPSR>){	
				my @columns = split(/\s+/,$_);
				if($type eq "all"){
					if ($columns[1] =~ /226032|274542|261044|260130|250847|226199|223187|197609|182849|178040|273628|255682|254315|226200|212030|184001|177427/){
						print OUT ">$hit\n$prots{$hit}\n";
					}}	
				elsif ($type eq "lipase"){
					if ($columns[1] =~ /250847|254315|226200/){
					print OUT ">$hit\n$prots{$hit}\n";}
				}
				elsif($type eq "hydrolase"){
					if ($columns[1] =~ /255682|226199|178040|182849|177427/){
						print OUT ">$hit\n$prots{$hit}\n"
					}}
				elsif($type eq "lysm"){
					if ($columns[1] =~ /212030|197609/){
						print OUT ">$hit\n$prots{$hit}\n"
					}}
				elsif($type =~ /ntpase|transferase/){
					if ($columns[1] =~ /260130|227061|224093|223187|273628|224228|184001/){
						print OUT ">$hit\n$prots{$hit}\n"
					}}
				elsif($type eq "unknown"){
					print OUT ">$hit\n$prots{$hit}\n"
					}
				else{next;}}
				$manager->finish;
				}
			close RPSR;
		close OUT;
		
	}	
	else {next;}
}
$manager->wait_all_children;
exit 0;
##################
sub temp_filename{
	    my $file = File::Temp->new(
	        TEMPLATE => 'tempXXXXX',
	        DIR      => '/tmp/',
	    );
	}
###
##!/usr/bin/bash
#for i in $(ls *.hmm);do;
#l=$(echo $i|  cut -c 17- | rev | cut -c 5- | rev )
#echo "";echo $l
#for j in $(ls ./Proteins/*.faa);do;
#k=$(echo $j| cut -c 12-)
#mkdir -p $l.out
#hmmscan --cpu 12 --tblout ./$l.out/$k.txt -o ./$l.out/$k.out $i $j
#echo -n "."
#done
#echo ""
#for i in $(ls | grep ".out");do
#./extractor.pl -hmm $i -prot ./Proteins -o  ./out -a
#done
#done
