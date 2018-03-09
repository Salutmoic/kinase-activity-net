#!/bin/bash

LOGDIR=log
MEM=24000
CORES=24

MERGED_DATA=out/kinase-merged-pred.tsv

VALSET=data/validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

bsub -J bart-build-$(basename $MERGED_DATA .tsv) \
     -o ${LOGDIR}/bart-build-$(basename $MERGED_DATA .tsv).out \
     -e ${LOGDIR}/bart-build-$(basename $MERGED_DATA .tsv).err \
     -M $MEM  -R "rusage[mem=16000]" -n $CORES \
     Rscript src/bart-build.r $MERGED_DATA $VALSET $NEG_VALSET FALSE FALSE
bsub -J bart-build-rand-negs-$(basename $MERGED_DATA .tsv) \
     -o ${LOGDIR}/bart-build-rand-negs$(basename $MERGED_DATA .tsv).out \
     -e ${LOGDIR}/bart-build-rand-negs$(basename $MERGED_DATA .tsv).err \
     -M $MEM  -R "rusage[mem=16000]" -n $CORES \
     Rscript src/bart-build.r $MERGED_DATA $VALSET $NEG_VALSET TRUE FALSE
bsub -J bart-build-direct-$(basename $MERGED_DATA .tsv) \
     -o ${LOGDIR}/bart-build-direct-$(basename $MERGED_DATA .tsv).out \
     -e ${LOGDIR}/bart-build-direct-$(basename $MERGED_DATA .tsv).err \
     -M $MEM  -R "rusage[mem=16000]" -n $CORES \
     Rscript src/bart-build.r $MERGED_DATA $VALSET $NEG_VALSET FALSE TRUE
