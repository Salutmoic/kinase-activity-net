#!/bin/bash

LOGDIR=log/bart-predict
MEM=24000
CORES=24
NUM=100

MERGED_DATA=out/kinase-merged-pred-cptac.tsv

VALSET=data/validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

# BART_MODELS="bart rand-negs-bart"
BART_MODELS="rand-negs"

for m in $BART_MODELS; do
         bsub -J "bart-predict-full-${m}-$(basename $MERGED_DATA .tsv)[1-${NUM}]" \
              -o ${LOGDIR}/bart-predict-full-${m}-$(basename $MERGED_DATA .tsv)-%I.out \
              -e ${LOGDIR}/bart-predict-full-${m}-$(basename $MERGED_DATA .tsv)-%I.err \
              -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
              Rscript src/bart-predict.r $MERGED_DATA out/$(basename $MERGED_DATA .tsv)-${m} '\$LSB_JOBINDEX'
done
