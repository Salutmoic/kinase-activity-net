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
         bsub -J bart-predict-full-${m}-$(basename $MERGED_DATA .tsv) \
              -o ${LOGDIR}/bart-predict-full-${m}-$(basename $MERGED_DATA .tsv).out \
              -e ${LOGDIR}/bart-predict-full-${m}-$(basename $MERGED_DATA .tsv).err \
              -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
              Rscript src/bart-predict.r $MERGED_DATA out/$(basename $MERGED_DATA .tsv)-${m}.Rdata
done
