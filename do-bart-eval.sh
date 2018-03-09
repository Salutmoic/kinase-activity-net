#!/bin/bash

LOGDIR=log
MEM=40000
CORES=24

MERGED_DATA=out/kinase-merged-pred.tsv

VALSET=data/validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

bsub -J bart-eval-$(basename $MERGED_DATA .tsv) \
     -o ${LOGDIR}/bart-eval-$(basename $MERGED_DATA .tsv).out \
     -e ${LOGDIR}/bart-eval-$(basename $MERGED_DATA .tsv).err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-eval.r $MERGED_DATA $VALSET $NEG_VALSET FALSE FALSE
bsub -J bart-eval-rand-negs-$(basename $MERGED_DATA .tsv) \
     -o ${LOGDIR}/bart-eval-rand-negs-$(basename $MERGED_DATA .tsv).out \
     -e ${LOGDIR}/bart-eval-rand-negs-$(basename $MERGED_DATA .tsv).err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-eval.r $MERGED_DATA $VALSET $NEG_VALSET TRUE FALSE
bsub -J bart-eval-direct-$(basename $MERGED_DATA .tsv) \
     -o ${LOGDIR}/bart-eval-direct-$(basename $MERGED_DATA .tsv).out \
     -e ${LOGDIR}/bart-eval-direct-$(basename $MERGED_DATA .tsv).err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-eval.r $MERGED_DATA $VALSET $NEG_VALSET FALSE TRUE
