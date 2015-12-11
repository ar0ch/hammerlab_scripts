#!/usr/bin/perl -w
# Aroon Chande
# Extracted HMM profiled proteins from a proteins.faa file
# Input os output of  hmmscan --tblout
# ./protein_extractor.pl -hmm hmmdir -prot proteindir -o outdir
use strict;
use Getopt::Long;
use File::Basename;
use File::Temp;
use Parallel::ForkManager;
#use network.pm;
my $prog = basename($0);
my ($h,$append,$threads,$type,$scan,$summary)=('0','0','10','all','0','0');
my ($inDir,$hmmDir,$protDir,$outdir,$hmmFile,$protFile,$outfile,$prots,@prots,%prots,@hits);
$prots = '';
if (@ARGV < 1){print_usage();exit 1;}
GetOptions ('hmm=s' => \$hmmDir, 'in=s' => \$inDir 'prot=s' => \$protDir, 'o=s' => \$outdir, 'h' => \$h, 't=i' => \$threads,'type=s' => \$type, 	'a' => \$append, 'scan' => \$scan);
if (eval $h){ print_usage();exit 1;}
my $manager = Parallel::ForkManager -> new ( $threads );

if($scan){
	my @hmms = glob ( "$hmmDir/*.hmm" );
	my @prots = glob ( "$protDir/*.faa" );
	foreach my $hmm (@hmms){
		my $base = basename($hmm);
		print STDERR "$base\n[";
		my $tblout = join('.',$base,"out");
		system(`mkdir -p $outdir/$tblout`);
		foreach my $pr (@prots){
			my $protOut = basename($pr);
			my $out = temp_filename();
			system(`hmmscan --cpu $threads -o $out --tblout $outdir/$tblout/$protOut $hmm $pr`);
			print STDERR ".";
		}
		print STDERR "]\n";
	}
}
opendir DIR, $inDir or die "cannot open dir $inDir: $!";
my @file= readdir DIR;
closedir DIR;
foreach (@file){
	@hits = ( );
	$prots = '';
	%prots = (	);
	if ($_ =~ /.txt$/){
		$hmmFile = join('/',$inDir,$_);
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
			$seq =~ s/\*[^ACDEFGHIJKLMNPQRSTVWY]//g;
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
		if ($summary && $append){
			open OUT, ">>$outfile" or die "Cannot open $outfile: $!\n";
			open HY, ">>$outdir/hydrolase.summary,txt" or die "Cannot open $outfile: $!\n";
			open LI, ">>$outdir/lipase.summary.txt" or die "Cannot open $outfile: $!\n";
			open LY, ">>$outdir/lysm.summary.txt" or die "Cannot open $outfile: $!\n";
			open NT, ">>$outdir/ntpase.transferase.summary.txt" or die "Cannot open $outfile: $!\n";
			open UN, ">>$outdir/unknown.summary.txt" or die "Cannot open $outfile: $!\n";
			}
		elsif ($summary) {
			open OUT, ">$outdir/hydrolase.summary,txt" or die "Cannot open $outfile: $!\n";
			open OUT, ">$outdir/lipase.summary.txt" or die "Cannot open $outfile: $!\n";
			open OUT, ">$outdir/lysm.summary.txt" or die "Cannot open $outfile: $!\n";
			open OUT, ">$outdir/ntpase.transferase.summary.txt" or die "Cannot open $outfile: $!\n";
			open OUT, ">$outdir/unknown.summary.txt" or die "Cannot open $outfile: $!\n";
			}
		elsif ($append){open OUT, ">>$outdir/hydrolase.summary,txt" or die "Cannot open $outfile: $!\n";}

		else{open OUT, ">$outdir/hydrolase.summary,txt" or die "Cannot open $outfile: $!\n";}		
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
#						print ALL ">$hit\n$prots{$hit}\n"; # funtionally the same as normal append only with -type all
					}}	
				elsif ($type eq "lipase"){
					if ($columns[1] =~ /250847|254315|226200/){
					print OUT ">$hit\n$prots{$hit}\n";}
					print LI ">$hit\n$prots{$hit}\n";
				}
				elsif($type eq "hydrolase"){
					if ($columns[1] =~ /255682|226199|178040|182849|177427/){
						print OUT ">$hit\n$prots{$hit}\n";
						print HY ">$hit\n$prots{$hit}\n";
					}}
				elsif($type eq "lysm"){
					if ($columns[1] =~ /212030|197609/){
						print OUT ">$hit\n$prots{$hit}\n";
						print LY ">$hit\n$prots{$hit}\n";
					}}
				elsif($type =~ /ntpase|transferase/){
					if ($columns[1] =~ /260130|227061|224093|223187|273628|224228|184001/){
						print OUT ">$hit\n$prots{$hit}\n";
						print NT ">$hit\n$prots{$hit}\n";
					}}
				elsif($type eq "unknown"){
					print OUT ">$hit\n$prots{$hit}\n";
					print UN ">$hit\n$prots{$hit}\n";
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
sub print_usage{
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
  -o    dir		Dumps extracted seqs to specifed directory
  -prot	dir		Directory containing protein files in .faa format
			Protein and hmm file must have same basename
  -t	int		Number of threads for rps-blast stage. Default: 10
  -type	str		Type of effector to verify. Default: all; types:
			all,hydrolase,lipase,lysm,ntpase,transferase,unknown


EXAMPLES
  $prog -hmm hmmscan_vgrg -prot reference_prots -o hits_out
  $prog -h

EXIT STATUS
  0     Successful completion
  >0    An error occurred

EOF
}

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
