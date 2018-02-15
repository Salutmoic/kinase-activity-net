#!/bin/bash

LOGDIR=log
MEM=24000
CORES=24

MERGED_DATA=out/kinase-merged-pred.tsv

VALSET=data/validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

bsub -J bart-build-full \
     -o ${LOGDIR}/bart-build-full.out -e ${LOGDIR}/bart-build-full.err \
     -M $MEM  -R "rusage[mem=16000]" -n $CORES \
     Rscript src/bart-build.r $MERGED_DATA $VALSET $NEG_VALSET FALSE FALSE
bsub -J bart-build-full-rand-negs \
     -o ${LOGDIR}/bart-build-full-rand-negs.out -e ${LOGDIR}/bart-build-full-rand-negs.err \
     -M $MEM  -R "rusage[mem=16000]" -n $CORES \
     Rscript src/bart-build.r $MERGED_DATA $VALSET $NEG_VALSET TRUE FALSE
bsub -J bart-build-full-direct \
     -o ${LOGDIR}/bart-build-full-direct.out -e ${LOGDIR}/bart-build-full-direct.err \
     -M $MEM  -R "rusage[mem=16000]" -n $CORES \
     Rscript src/bart-build.r $MERGED_DATA $VALSET $NEG_VALSET FALSE TRUE
