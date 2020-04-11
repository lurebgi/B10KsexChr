genome=$1
marker=$2
genome1=$(echo ${genome##*/})

blat $genomo $marker  $genome1.$marker.blatout -noHead
cat $genome1.$marker.blatout | awk '$1/$11>0.2{ t+=($1+$12)/$11}END{print "'$genome1'\t'$marker'"t}' >> $genome1.blatout.cov
