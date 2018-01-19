#!/bin/bash

LOGDIR=log
MEM=32000
CORES=24

FULL_DATA=out/kinact-merged_full.tsv
ASSOC_DATA=out/kinact-max-rows-merged.tsv
DIRECT_DATA=out/kinact-max-rows-direct.tsv

VALSET=data/validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

# bsub -J bart-eval \
#      -o ${LOGDIR}/bart-eval.out -e ${LOGDIR}/bart-eval.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-eval.r $ASSOC_DATA $VALSET $NEG_VALSET FALSE FALSE FALSE
# bsub -J bart-eval-rand-negs -o ${LOGDIR}/bart-eval.out -e ${LOGDIR}/bart-eval.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-eval.r $ASSOC_DATA $VALSET $NEG_VALSET TRUE FALSE FALSE

bsub -J bart-eval-full \
     -o ${LOGDIR}/bart-eval-full.out -e ${LOGDIR}/bart-eval-full.err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-eval.r $FULL_DATA $VALSET $NEG_VALSET FALSE FALSE FALSE
bsub -J bart-eval-full-rand-negs \
     -o ${LOGDIR}/bart-eval-full-rand-negs.out -e ${LOGDIR}/bart-eval-full-rand-negs.err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-eval.r $FULL_DATA $VALSET $NEG_VALSET TRUE FALSE FALSE

# bsub -J bart-eval-direct \
#      -o ${LOGDIR}/bart-eval-direct.out -e ${LOGDIR}/bart-eval-direct.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-eval.r $DIRECT_DATA $VALSET $NEG_VALSET FALSE FALSE TRUE
# bsub -J bart-eval-direct-rand-negs \
#      -o ${LOGDIR}/bart-eval-direct-rand-negs.out -e ${LOGDIR}/bart-eval-direct-rand-negs.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-eval.r $DIRECT_DATA $VALSET $NEG_VALSET TRUE FALSE TRUE
