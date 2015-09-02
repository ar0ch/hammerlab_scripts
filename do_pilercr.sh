for i in $(find . -type f  -name '*.fasta')
do
echo $i run through pilercr
$PILERPATH/pilercr -in $i -out $i.out -seq $i.seq
done
