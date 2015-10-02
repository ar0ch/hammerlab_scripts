#!/bin/bash
#Needle for matrix
i=1 
j=1 
length=$(wc -l $1)
while [ $i -le 577 ]; do
	seq=$(awk -v pairs=$1 -v i=$i 'FNR == i {printf("-asequence effector/"$1".txt -bsequence effector/"$2".txt")}' pairs)
	echo $seq
	name=$(awk -v pairs=$1 -v j=$j 'FNR == j {printf($1"."$2)}' pairs)
	needle $seq -datafile ~/bin/data/EBLOSUM60 --gapopen 10 -gapextend 0.5 -outfile $name -adirectory3 needle/ -auto -nobrief
	sed  s/\#//g needle/$name > temp/$name
	awk 'FNR == 21{printf $2" "}; FNR==28{print $2" "}; FNR==22{printf $2" "}' temp/$name | tr -d '()' |  sed s'|/| |'g > temp/$name.awk.temp
	awk '{print($1,1/($3/$4),$2)}' temp/$name.awk.temp >> effectornodes.txt
        let i=$[$i +1]
        let j=$[$j +1]
done
