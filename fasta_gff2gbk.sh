#!/bin/sh
for i in $(ls $1 | sed  s/.reordered.fasta//g)
do
seqret -sequence $1/$i.reordered.fasta -feature -fformat gff -fopenfile ~/hammer_lab/GFF/$i.gff  -osformat genbank -auto -osname2 $i
echo $i converted to genbank format
done

