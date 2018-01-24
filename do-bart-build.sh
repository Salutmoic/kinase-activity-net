#!/bin/bash

LOGDIR=log
MEM=16000
CORES=24

FULL_DATA=out/kinact-merged_full.tsv
ASSOC_DATA=out/kinact-max-rows-merged.tsv
DIRECT_DATA=out/kinact-max-rows-direct.tsv

VALSET=data/validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

# bsub -J bart-build \
#      -o ${LOGDIR}/bart-build.out -e ${LOGDIR}/bart-build.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-build.r $ASSOC_DATA $VALSET $NEG_VALSET FALSE FALSE
# bsub -J bart-build-rand-negs -o ${LOGDIR}/bart-build.out -e ${LOGDIR}/bart-build.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-build.r $ASSOC_DATA $VALSET $NEG_VALSET TRUE FALSE

bsub -J bart-build-full \
     -o ${LOGDIR}/bart-build-full.out -e ${LOGDIR}/bart-build-full.err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-build.r $FULL_DATA $VALSET $NEG_VALSET FALSE FALSE
bsub -J bart-build-full-rand-negs \
     -o ${LOGDIR}/bart-build-full-rand-negs.out -e ${LOGDIR}/bart-build-full-rand-negs.err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-build.r $FULL_DATA $VALSET $NEG_VALSET TRUE FALSE
bsub -J bart-build-full-direct \
     -o ${LOGDIR}/bart-build-full-direct.out -e ${LOGDIR}/bart-build-full-direct.err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-build.r $FULL_DATA $VALSET $NEG_VALSET FALSE TRUE

# bsub -J bart-build-direct \
#      -o ${LOGDIR}/bart-build-direct.out -e ${LOGDIR}/bart-build-direct.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-build.r $DIRECT_DATA $VALSET $NEG_VALSET FALSE TRUE
# bsub -J bart-build-direct-rand-negs \
#      -o ${LOGDIR}/bart-build-direct-rand-negs.out -e ${LOGDIR}/bart-build-direct-rand-negs.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-build.r $DIRECT_DATA $VALSET $NEG_VALSET TRUE TRUE
