ref=$1
spe=$2

export PYTHONPATH=$PYTHONPATH:/scratch/luohao/software/RaGOO-1.11/lib/python3.7/site-packages/
~/miniconda3/bin/python3.7 /scratch/luohao/software/RaGOO-1.11/bin/ragoo.py -C  $spe.z-scf.fa $ref.chrZ.fa -m /apps/minimap2/2.17/minimap2

mv ragoo_output $spe.ragoo_output

zID=`head -1 emu.chrZ.fa | awk '{print $1}' | sed 's/>//' `
cat $spe.ragoo_output/orderings/${zID}_orderings.txt  | cut -f 1,2 | awk 'BEGIN{st=0;while(getline < "'$spe'.z-scf.fa.fai"){a[$1]=$2}}{print $0"\t"st"\t"st+a[$1];st+=a[$1]+100}' | awk '{st[$1]=$3;end[$1]=$4;str[$1]=$2}END{while(getline < "'$spe'.cov-50k"){if(end[$1]>0 && str[$1]=="+"){print "chrZ\t"st[$1]+$2"\t"st[$1]+$3"\t"$5"\t"$6"\t"$1}if(end[$1]>0 && str[$1]=="-"){print "chrZ\t"end[$1]-$3"\t"end[$1]-$2"\t"$5"\t"$6"\t"$1} }}'   | sort -k2n > $spe.z.chr-cov
