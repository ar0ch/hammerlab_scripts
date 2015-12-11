#!/bin/bash
#Water for matrix
verbose="0"
pairs="pairs"
while getopts ":pv" opts; do
	case $opts in
		p) #Set pairs file
		pairs=$OPTARG;
		;;
		v) #Set verbose flag
		verbose="1"
		;;
	esac
done
i=1 
j=1 
length=$(wc -l pairs | awk '{print $1}')
timestamp=$(date -u +%Y%m%d_%H%M)
mkdir -p $1/temp $1/water
rm -rf $1/water_nodes_$timestamp.txt > /dev/null 2>&1
while [[ $i -le $length ]]; do
	# emoss suite is not ammenable to direct var passing, construct entire options arg in `seq`
	seq=$(awk -v folder=$1 -v i=$i 'FNR == i {printf("-asequence "folder"/"$1".fasta -bsequence "folder"/"$2".fasta")}' pairs)
	name=$(awk  -v j=$j 'FNR == j {printf($1"."$2)}' pairs)
	water $seq -datafile EBLOSUM60 -gapopen 10 -gapextend 0.5 -outfile $1/water/$name -nobrief > /dev/null 2>&1 
	echo -n "."
	sed  s/\#//g $1/water/$name > $1/temp/$name
	awk 'FNR == 19{printf $2" "}; FNR==26{print $2" "}; FNR==20{printf $2" "}'  $1/temp/$name  | tr -d '()' |  sed s'|/| |'g  > $1/temp/$name.awk.temp
	awk '{print($1,1-($3/$4),$2)}' $1/temp/$name.awk.temp >> $1/water_nodes_$timestamp.txt
        let i=$[$i +1]
        let j=$[$j +1]
done
rm -rf $1/temp $1/water

