genome_ori=$1
read1=$2
read2=$3
sample=female
cpu=8

genome=$(echo ${genome_ori##*/})

[ -d index ] || mkdir -p index
if [ ! -f index/${genome}.bwt ] ; then
bwa index $genome_ori -p index/$genome
fi

bwa mem -t $cpu -R $(echo "@RG\tID:$sample\tSM:$sample\tLB:$sample"_"$sample\tPL:ILLUMINA")   index/$genome $read1 $read2  |  samtools sort -@ $cpu -O BAM -o $genome.$sample.sorted.bam  -
samtools index -@ $cpu $genome.$sample.sorted.bam

# mark duplication, please modify the path of picard.jar
java -Xmx108g -jar /apps/picard/2.21.4/picard.jar  MarkDuplicates I=$genome.$sample.sorted.bam O=$genome.$sample.dedup.bam M=$sample.m
samtools index -@ $cpu  $genome.$sample.dedup.bam

[ -f $genome.$sample.dedup.bam.bai ] && rm $genome.$sample.sorted.bam $genome.$sample.sorted.bam.bai $sample.m
