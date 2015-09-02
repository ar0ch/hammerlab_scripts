#!/bin/sh
## Extract unalign coordinates from promer or nucmer mapping
## using GFF. NB: If using bacterial GFF with `mapview`,
## replace CDS with single-exon sed -i s/CDS/single-exon/g
## Aroon Chande, Hammer Lab

length=$(wc -l $1 | awk '{print $1}')
echo $length
i=1
j=2
while [ $j -lt $length ]; do
        awk -v i=$i 'FNR == i {printf $10" "}' $1 >>$1.unaligned #print chromosome
        awk -v j=$j 'FNR == j {printf $1" "}' $1 >> $1.unaligned #print s
        awk -v i=$i 'FNR == i {print $2}' $1 >>$1.unaligned
        let i=$[$i +1]
        let j=$[$j+1]
done
