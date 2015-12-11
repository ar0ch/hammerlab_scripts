#!/usr/bin/perl -w
# Aroon Chande
# Prepare files for network analysis
use strict;
use Math::Combinatorics;
use Getopt::Long;
use File::Basename;
use File::Temp;
my $prog = basename($0);
my($max,$protDir,$outDir,$row,$col,%matrix);
if (@ARGV < 1){print_usage();exit 1;}
my ($h,$outmat,$tsv)=('0','matrix.tsv','nodes.tsv');
GetOptions ( 'prot=s' => \$protDir, 'o=s' => \$outDir, 'h' => \$h, 'outmat=s' => \$outmat, 'tsv=s' => \$tsv);
if (eval $h){ print_usage();exit 1;}
my @file = glob ( "$protDir/*" );
for($row=0; $row<@file;$row++){
	for($col=$row+1; $col<@file;$col++){
		my $out = temp_filename();
		system(`water -asequence $file[$row] -bsequence $file[$col] -datafile EBLOSUM60 -gapopen 10 -gapextend 0.5 -outfile  $out -nobrief > /dev/null 2>&1`);
		my $id = `head -26 $out | tail -1`;
		chomp $id;
		$id =~ s/.*\(//;
		$id =~ s/%.*//;
		my $rowName = `head -1 $file[$row] | sed 's/>//'`;
		my $colName = `head -1 $file[$col] | sed 's/>//'`;
		chomp $rowName;
		chomp $colName;
		$id = 100 - $id;		
		$id = 1e-20 if($id == 0);
		$matrix{$rowName}{$colName} = $id;
		$matrix{$colName}{$rowName}= $id;
		$matrix{$rowName}{$rowName} = 0;
		$matrix{$colName}{$colName} = 0;
		}
}

my @nodes = sort keys %matrix;
open MAT, ">$outmat" or die "Cannot open $outmat: $!\n";
for (my $col = 0; $col < @nodes; $col++){
	print MAT "\t$nodes[$col]";
}
print MAT "\n";
for (my $row = 0; $row < @nodes; $row++){
	print MAT "\t$nodes[$row]";
	for (my $col = 0; $col < @nodes; $col++){
		print MAT "\t$matrix{$nodes[$row]}{$nodes[$col]}";
	}
	print MAT "\n";
}

close MAT;

open NODE, ">$tsv" or die "Cannot open $tsv: $!\n";
print NODE "node1\tc1\tt1\tdist\tnode2\tc2\tt2\n";
for (my $row = 0; $row < @nodes; $row++){
	my (undef, $c1, $t1) = split(/_/, $nodes[$row]);
	for (my $col = 0; $col < @nodes; $col++){
		my (undef, $c2, $t2) = split(/_/, $nodes[$col]);
		print NODE "$nodes[$row]\t$c1\t$t1\t$matrix{$nodes[$row]}{$nodes[$col]}\t$nodes[$col]\t$c2\t$t2\n";
	}
}
close NODE;
#####################################################################
#####################################################################

sub print_usage
{
    warn <<"EOF";
USAGE
  $prog -prot <protdir> [options]
DESCRIPTION
  This program builds a matrix of sequence identity and exports
  a tsv in cytoscape importable format. Input files must be not
  be multifasta. 
OPTIONS
  -h    		Print this help message
  -prot	dir		Directory containing protein files in fasta format
  			Description line should be formatted as shown:
  			>Strain_locus_function
  -outmat str		Output file for matrix
  -tsv	  str		Output file for nodes TSV for cytoscape import
EXAMPLES
  $prog -prot aux1/tap/ -outmat aux1.tap.matrix.tsv -tsv aux1.tap.tsv
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

	
