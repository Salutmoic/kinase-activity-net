#!/bin/bash

LOGDIR=log/bart-lo-kin
MEM=$(echo "256*1024" | bc)
CORES=24

MERGED_FULL_DATA=out/kinase-merged-feats-clean.tsv

DIRECT_VALSET=data/direct-validation-set-omnipath.tsv

NUM=$(cat <(cut -f1 $DIRECT_VALSET) <(cut -f2 $DIRECT_VALSET) | sort | uniq | wc -l)

mkdir -p $LOGDIR

# bsub -J "bart-lo-kin-$(basename $MERGED_FULL_DATA .tsv)[1-${NUM}]" \
#      -o ${LOGDIR}/bart-lo-kin-$(basename $MERGED_FULL_DATA .tsv)-%I.out \
#      -e ${LOGDIR}/bart-lo-kin-$(basename $MERGED_FULL_DATA .tsv)-%I.err \
#      -M $MEM -P bigmem -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-lo-kin.r $MERGED_FULL_DATA $DIRECT_VALSET '\$LSB_JOBINDEX'

for n in 178; do
bsub -J "bart-lo-kin-$(basename $MERGED_FULL_DATA .tsv)-${n}" \
     -o ${LOGDIR}/bart-lo-kin-$(basename $MERGED_FULL_DATA .tsv)-${n}.out \
     -e ${LOGDIR}/bart-lo-kin-$(basename $MERGED_FULL_DATA .tsv)-${n}.err \
     -M $MEM -P bigmem -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-lo-kin.r $MERGED_FULL_DATA $DIRECT_VALSET $n
done

