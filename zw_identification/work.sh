# create a list containing the paths of genomes and alignments
cat b10k.list | while read id; do echo -e "$id\t/hwfssz5/ST_DIVERSITY/B10K/PMO/2.Family_level/01.Assembly/$id/$id.genomic.fa\t/zfssz3/ST_DIVERSITY/F16ZQSB1SY2933/PMO/B10K_psmc/01.famliy_level/01.recal.bam/$id.clean.filter.sort.dedup.AddGroup.realign.recal.bam";done > b10k.list.genome-bam

b10k.list.genome-bam | while read id genome bam; do sh half-coverage $genome $bam; done
