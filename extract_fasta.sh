#!/bin/sh
## Extract fasta sequences using unalign coordinates from
## MUMmer, extracted using extract_unaligned.coords.sh
## Aroon Chande, Hammer Lab
sed ':a;N;$!ba;s/\n//g' $2 > /tmp/nospaces.fasta
while IFS=' ' read col1 col2; do
        let delim=$col1-$col2
        echo '>'sequence from coordinate $col1 - $col2 >> fasta.txt
        cut -c $delim /tmp/nospaces.fasta >> fasta.txt
done < $1
rm /tmp/nospaces.fasta

