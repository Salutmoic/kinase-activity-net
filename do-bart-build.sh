#!/bin/bash

LOGDIR=log/bart-build
MEM=$(echo "32*1024" | bc)
CORES=24
NUM=100

MERGED_ASSOC_DATA=out/kinase-association-feats.tsv
MERGED_DIRECT_DATA=out/kinase-direction-feats.tsv

ASSOC_VALSET=data/assoc-validation-set-omnipath.tsv
DIRECT_VALSET=data/direct-validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

bsub -J "bart-build-$(basename $MERGED_ASSOC_DATA .tsv)[1-${NUM}]" \
     -o ${LOGDIR}/bart-build-$(basename $MERGED_ASSOC_DATA .tsv)-%I.out \
     -e ${LOGDIR}/bart-build-$(basename $MERGED_ASSOC_DATA .tsv)-%I.err \
     -M $MEM  -R "rusage[mem=16000]" -n $CORES \
     Rscript src/bart-build.r $MERGED_ASSOC_DATA $ASSOC_VALSET $NEG_VALSET TRUE FALSE '\$LSB_JOBINDEX'

bsub -J "bart-build-$(basename $MERGED_DIRECT_DATA .tsv)[1-${NUM}]" \
     -o ${LOGDIR}/bart-build-$(basename $MERGED_DIRECT_DATA .tsv)-%I.out \
     -e ${LOGDIR}/bart-build-$(basename $MERGED_DIRECT_DATA .tsv)-%I.err \
     -M $MEM  -R "rusage[mem=16000]" -n $CORES \
     Rscript src/bart-build.r $MERGED_DIRECT_DATA $DIRECT_VALSET $NEG_VALSET FALSE TRUE '\$LSB_JOBINDEX'

# for n in 51 52 10 6 25 56 63 14 78; do
# bsub -J "bart-build-$(basename $MERGED_DIRECT_DATA .tsv)-${n}" \
#      -o ${LOGDIR}/bart-build-$(basename $MERGED_DIRECT_DATA .tsv)-${n}.out \
#      -e ${LOGDIR}/bart-build-$(basename $MERGED_DIRECT_DATA .tsv)-${n}.err \
#      -M $MEM  -R "rusage[mem=16000]" -n $CORES \
#      Rscript src/bart-build.r $MERGED_DIRECT_DATA $DIRECT_VALSET $NEG_VALSET FALSE TRUE $n
# done
