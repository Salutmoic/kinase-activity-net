#!/bin/bash

LOGDIR=log
MEM=32000
CORES=24

MERGED_DATA=out/kinase-merged-pred.tsv

VALSET=data/validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

bsub -J bart-eval-full \
     -o ${LOGDIR}/bart-eval-full.out -e ${LOGDIR}/bart-eval-full.err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-eval.r $MERGED_DATA $VALSET $NEG_VALSET FALSE FALSE
bsub -J bart-eval-full-rand-negs \
     -o ${LOGDIR}/bart-eval-full-rand-negs.out -e ${LOGDIR}/bart-eval-full-rand-negs.err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-eval.r $MERGED_DATA $VALSET $NEG_VALSET TRUE FALSE
bsub -J bart-eval-full-direct \
     -o ${LOGDIR}/bart-eval-full-direct.out -e ${LOGDIR}/bart-eval-full-direct.err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-eval.r $MERGED_DATA $VALSET $NEG_VALSET FALSE TRUE
