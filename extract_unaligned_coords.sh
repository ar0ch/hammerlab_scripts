#!/bin/sh
## Extract unalign coordinates from promer or nucmer mapping
## using GFF. NB: If using bacterial GFF with `mapview`,
## replace `CDS` with `single-exon`
length=$(wc -l $1 | awk '{print $1}')
echo $length
i=1
j=2
while [ $j -lt $length ]; do
        awk -v j=$j 'FNR == j {printf $1" "}' $1 >> $1.unaligned && awk -v i=$i 'FNR == i {print $2}' $1 >>$1.unaligned
        let i=$[$i +i]
        let j=$[$j+1]
done
