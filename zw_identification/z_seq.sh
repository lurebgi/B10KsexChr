ref_ori=$1 ## reference genome
spe_ori=$2 ## the genome
ref_z=$3 ## list of chr/scaffold ID of the Z chromosome
cpu=4


ref=$(echo ${ref_ori##*/})
spe=$(echo ${spe_ori##*/})

cat $ref_ori | sed 's/NN/XX/g' > $ref.x # mark the gaps
ln -s $spe_ori
samtools faidx $spe

nucmer -b 400  -t $cpu $ref.x  $spe -p ${ref}-$spe

delta-filter -1 -l 400 ${ref}-$spe.delta > ${ref}-$spe.delta.filt
show-coords -H -c -l -o -r -T ${ref}-$spe.delta.filt > ${ref}-$spe.delta.filt.coords

cat  ${ref}-$spe.delta.filt.coords  | grep -w -f $ref_z - | awk '$5>500' | awk '$3<$4{print $13"\t"$3"\t"$4}$3>$4{print $13"\t"$4"\t"$3}' |  bedtools sort -i - | bedtools merge -i - -d 10000 |  awk 'BEGIN{while(getline < "'$spe'.fai"){b[$1]=$2}}{a[$1]+=$3-$2+1}END{for(i in a){print i"\t"a[i]/b[i]*100"\t"b[i]}}' > ${ref}-$spe.delta.filt.coords.Z-perc

cutoff=60 # the cutoff is to be decided according to the distribution of ${ref}-$spe.delta.filt.coords.Z-perc
cat ${ref}-$spe.delta.filt.coords.Z-perc | awk -v cutoff=$cutoff '$3>cutoff' > $spe.z-list