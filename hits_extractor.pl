#!/usr/bin/perl -w
# Aroon Chande
# Extracted HMM profiled proteins from a proteins.faa file
# Input is output of  hmmscan --tblout
# ./protein_extractor.pl -hmm hmmdir -prot proteindir -o outdir
use strict;
use Getopt::Long;
use File::Basename;
use File::Temp;
use Parallel::ForkManager;
use Pod::Usage;
#use network.pm;
my $prog = basename($0);
my ($h,$append,$threads,$scan,$summary,$rpsdb,$hmmDir)=('0','0','10','0','0','/home/achande3/bin/db/cdd/Cdd','/home/achande3/hmm/hmm_profiles');
my ($inDir,$protDir,$outdir,$hmmFile,$protFile,$outfile,$prots,@prots,%prots,@hits,$check,$type);
my @type = ( "hcp", "hydrolase", "lipase", "lysm", "ntpase", "vgrg" );
$prots = '';
if (@ARGV < 1){print_usage();exit 1;}
GetOptions ('hmm=s' => \$hmmDir, 'in=s' => \$inDir, 'prot=s' => \$protDir, 'o=s' => \$outdir, 'h' => \$h, 
				't=i' => \$threads,'type=s{1,}' => \$type, 'a' => \$append, 'scan' => \$scan, 'summary' => \$summary, 'db=s' => \$rpsdb,
				'no-check' => \$check);
if (defined $type){@type = split(/,/,$type);}
#unless (defined $check)( die print_usage() unless ((defined $hmmDir) && ($scan)) && defined $outdir || (defined $inDir && !($scan)));
if (eval $h){ print_usage();exit 1;}
my $manager = Parallel::ForkManager -> new ( $threads );
if($scan){
	my @hmms = glob ( "$hmmDir/*.hmm" );
	my @prots = glob ( "$protDir/*.faa" );
	foreach my $hmm (@type){
		print STDERR "Running hmmscans:\n$hmm\n[";
		my $tblout = join('.',$hmm,"out");
		my $hmm = join(".","$hmmDir/$hmm","hmm");
		system(`mkdir -p $outdir/$tblout`) or die "Cannot create directory: $outdir/$tblout\n";
		#$inDir="$outdir/$tblout/";
		foreach my $pr (@prots){
			$manager->start and next;
			my $protOut = join('.',basename($pr),"txt");
			my $out = temp_filename();
			system(`hmmscan --cpu $threads -o $out --tblout $outdir/$tblout/$protOut $hmm $pr`) or die "Cannot parse an input option: $!\n";
			print STDERR ".";
			$manager->finish;		
		}
		$manager->wait_all_children;
		print STDERR "]\n";
	}
}
my @file= glob( "$outdir/*.out/*");
print "Extracting and verifying hits, filtering by: @type\n";
foreach (@file){
	@hits = ( );
	$prots = '';
	%prots = (	);
	if ($_ =~ /.txt$/){
		$hmmFile = $_;
		$_ = basename($_);
		$_ =~ s/\.txt//g;
		$protFile = join('/',$protDir,$_);
		$outfile = join('/',$outdir,$_);	
		my ($strain) = $_ =~ m/([a-zA-Z]*[0-9]*[^.])/;
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
		system(`mkdir -p $outdir`);
		 if ($summary && $append){
			open OUT, ">>$outfile" or die "Cannot open $outfile: $!\n" ;
			open HY, ">>$outdir/hydrolase.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/hydrolase/i, @type));
			open LI, ">>$outdir/lipase.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/lipase/i, @type));
			open LY, ">>$outdir/lysm.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/lysm/i, @type));
			open NT, ">>$outdir/ntpase.transferase.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/ntpase|transferase/i, @type));
			#open UN, ">>$outdir/unknown.summary.txt" or die "Cannot open $outfile: $!\n";
			open VG, ">>$outdir/vgrg.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/vgrg/i, @type));
			open HCP, ">>$outdir/hcp.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/hcp/i, @type));
			}
		elsif ($summary) {
			open OUT, ">>$outfile" or die "Cannot open $outfile: $!\n" ;
			open HY, ">$outdir/hydrolase.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/hydrolase/i, @type));
			open LI, ">$outdir/lipase.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/lipase/i, @type));
			open LY, ">$outdir/lysm.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/lysm/i, @type));
			open NT, ">$outdir/ntpase.transferase.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/ntpase|transferase/i, @type));
			#open UN, ">$outdir/unknown.summary.txt" or die "Cannot open $outfile: $!\n";
			open VG, ">$outdir/vgrg.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/vgrg/i, @type));
			open HCP, ">$outdir/hcp.summary.txt" or die "Cannot open $outfile: $!\n" if(grep(/hcp/i, @type));
			}
		 else{open OUT, ">>$outfile" or die "Cannot open $outfile: $!\n";}	
		#open OUT, ">>", "$outfile" or die "Cannot open $outfile: $!\n";
		foreach my $hit (@hits){
			$manager->start and next;
			my $rpstemp=temp_filename();
			my $rpsresulttemp=temp_filename();
			open RPS, ">$rpstemp" or die "Cannot open $rpstemp: $!";
			print RPS ">$hit\n$prots{$hit}\n";
			close RPS;
			`rpsblast -db $rpsdb -max_target_seqs 1 -outfmt "6" -out $rpsresulttemp -query $rpstemp`;
			open RPSR, "$rpsresulttemp" or die "Cannot open $rpsresulttemp: $!";
			foreach (<RPSR>){	
				my @columns = split(/\s+/,$_);
				if ($columns[1] =~ /250847|254315|226200/){
					print OUT ">$hit|$strain|lipase\n$prots{$hit}" if(grep(/lipase/i, @type));
					print LI ">$hit|$strain|lipase\n$prots{$hit}" if((grep(/lipase/i, @type)) && ($summary));
				}
				elsif($columns[1] =~ /255682|226199|178040|182849|177427/){
						print OUT ">$hit|$strain|Hydrolase\n$prots{$hit}"if(grep(/hydrolase/i, @type));
						print HY ">$hit|$strain|Hydrolase\n$prots{$hit}" if((grep(/hydrolase/i, @type)) && ($summary));
					}
				elsif($columns[1] =~ /212030|197609/){
						print OUT ">$hit|$strain|LysM\n$prots{$hit}" if(grep(/lysm/i, @type));
						print LY ">$hit|$strain|LysM\n$prots{$hit}" if((grep(/lysm/i, @type)) && ($summary));
					}
				elsif($columns[1] =~ /260130|227061|224093|223187|273628|224228|184001/){
						print OUT ">$hit|$strain|ntpase\n$prots{$hit}" if(grep(/ntpase|transferase/i, @type));
						print NT ">$hit|$strain|ntpase\n$prots{$hit}" if ((grep(/ntpase|transferase/i, @type)) && ($summary));
					}
				elsif($columns[1] =~ /274542|224482/){
					print OUT ">$hit|$strain|vgrg\n$prots{$hit}" if(grep(/vgrg/i, @type));
					print VG ">$hit|$strain|vgrg\n$prots{$hit}" if((grep(/vgrg/i, @type)) && ($summary));
					}
				elsif($columns[1] =~ /213799/){
					print OUT ">$hit|$strain|hcp\n$prots{$hit}" if(grep(/hcp/i, @type));
					print HCP ">$hit|$strain|hcp\n$prots{$hit}" if((grep(/hcp/i, @type)) && ($summary));
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
  -a			Append data to same protein.faa file
  -h    		Print this help message
  -hmm	dir		Directory containing hmm profiles (profile.hmm)
  -in	dir		Directory with tab delimited HMMscan output
  -o    dir		Dumps extracted seqs to specifed directory
  -prot	dir		Directory containing protein files in .faa format
				Protein and hmm file must have same basename
  -scan 		Run hmmscan, requires -hmm. If set, -in is ignored
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
