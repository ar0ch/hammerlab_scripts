for i in $(find . -type f  -name '*.fasta')
do
echo $i run through pilercr
/root/ReorderedContigs/pilercr -in $i -out $i.out -seq $i.seq
done
