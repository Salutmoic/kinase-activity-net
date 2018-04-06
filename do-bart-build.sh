#!/bin/bash

LOGDIR=log/bart-build
MEM=$(echo "48*1024" | bc)
CORES=24
NUM=100

MERGED_DATA=out/kinase-merged-pred-cptac-nogo.tsv

VALSET=data/validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

# bsub -J "bart-build-$(basename $MERGED_DATA .tsv)[1-${NUM}]" \
#      -o ${LOGDIR}/bart-build-$(basename $MERGED_DATA .tsv)-%I.out \
#      -e ${LOGDIR}/bart-build-$(basename $MERGED_DATA .tsv)-%I.err \
#      -M $MEM  -R "rusage[mem=16000]" -n $CORES \
#      Rscript src/bart-build.r $MERGED_DATA $VALSET $NEG_VALSET FALSE FALSE '\$LSB_JOBINDEX'
bsub -J "bart-build-rand-negs-$(basename $MERGED_DATA .tsv)[1-${NUM}]" \
     -o ${LOGDIR}/bart-build-rand-negs-$(basename $MERGED_DATA .tsv)-%I.out \
     -e ${LOGDIR}/bart-build-rand-negs-$(basename $MERGED_DATA .tsv)-%I.err \
     -M $MEM  -R "rusage[mem=16000]" -n $CORES \
     Rscript src/bart-build.r $MERGED_DATA $VALSET $NEG_VALSET TRUE FALSE '\$LSB_JOBINDEX'
# bsub -J "bart-build-direct-$(basename $MERGED_DATA .tsv)[1-${NUM}]" \
#      -o ${LOGDIR}/bart-build-direct-$(basename $MERGED_DATA .tsv)-%I.out \
#      -e ${LOGDIR}/bart-build-direct-$(basename $MERGED_DATA .tsv)-%I.err \
#      -M $MEM  -R "rusage[mem=16000]" -n $CORES \
#      Rscript src/bart-build.r $MERGED_DATA $VALSET $NEG_VALSET FALSE TRUE '\$LSB_JOBINDEX'

