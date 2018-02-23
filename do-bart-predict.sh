#!/bin/bash

LOGDIR=log
MEM=16000
CORES=24

MERGED_DATA=out/kinase-merged-pred.tsv

VALSET=data/validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

BART_MODELS="bart rand-negs-bart"

for m in $BART_MODELS; do
         bsub -J bart-predict-full-${m} \
              -o ${LOGDIR}/bart-predict-full-${m}.out -e ${LOGDIR}/bart-predict-full-${m}.err \
              -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
              Rscript src/bart-predict.r $MERGED_DATA out/kinase-merged-pred-${m}.Rdata
done
