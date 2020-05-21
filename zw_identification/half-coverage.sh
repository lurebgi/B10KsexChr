
genome_ori=$1
bam=$2
genome=$(echo ${genome0##*/})
ln -s $genome_ori

samtools faidx $genome
cut -f 1,2 $genome.fai > $genome.fai.g
bedtools makewindows -g $genome.fai.g -w 50000 > $genome.fai.g.50k

samtools depth -m 100 -Q 60  $bam | awk '{print $1"\t"$2-1"\t"$2"\t"$3}'  | bedtools map -a $genome.fai.g.50k -g $genome.fai.g  -b - -c 4 -o median,mean,count  >  $bam.cov-50k

medianCov=$(cat $bam.cov-50k | awk '$6/($3-$2)>0.8' | cut -f 5 | sort -n | awk  { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; })

cat $bam.cov-50k | awk '$6/($3-$2)>0.6 && $3-$2>2000{a[$1]+=$3-$2+1}END{for(i in a){print i"\t"a[i]}}' > $genome.2k-length

cat $bam.cov-50k | awk '$6/($3-$2)>0.6 && $3-$2>2000' | awk -v medianCov=$medianCov '$5>medianCov/3 && $5<medianCov/3*2' | awk 'BEGIN{while(getline < "'$genome'.2k-length"){len[$1]=$2}}{a[$1]+=$3-$2+1}END{for(i in a){print i"\t"a[i]/len[i]*100}}' > $bam.cov-50k.half-perc

cat $bam.cov-50k.half-perc | awk '$2>0.8' > $bam.cov-50k.half-perc.list
