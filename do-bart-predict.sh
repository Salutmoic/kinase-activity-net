#!/bin/bash

LOGDIR=log/bart-predict
MEM=24000
CORES=24
NUM=100

MERGED_ASSOC_DATA=out/kinase-association-feats.tsv
MERGED_DIRECT_DATA=out/kinase-direction-feats.tsv
MERGED_FULL_DATA=out/kinase-merged-feats-clean-hinegs.tsv

ASSOC_VALSET=data/assoc-validation-set-omnipath.tsv
DIRECT_VALSET=data/direct-validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

bsub -J "bart-predict-$(basename $MERGED_FULL_DATA .tsv)[1-${NUM}]" \
     -o ${LOGDIR}/bart-predict-$(basename $MERGED_FULL_DATA .tsv)-%I.out \
     -e ${LOGDIR}/bart-predict-$(basename $MERGED_FULL_DATA .tsv)-%I.err \
     -M $MEM  -R "rusage[mem=12000]" -n $CORES \
     Rscript src/bart-predict.r $MERGED_FULL_DATA out/$(basename $MERGED_FULL_DATA .tsv) FALSE '\$LSB_JOBINDEX'

# bsub -J "bart-predict-$(basename $MERGED_ASSOC_DATA .tsv)[1-${NUM}]" \
#      -o ${LOGDIR}/bart-predict-$(basename $MERGED_ASSOC_DATA .tsv)-%I.out \
#      -e ${LOGDIR}/bart-predict-$(basename $MERGED_ASSOC_DATA .tsv)-%I.err \
#      -M $MEM  -R "rusage[mem=12000]" -n $CORES \
#      Rscript src/bart-predict.r $MERGED_ASSOC_DATA out/$(basename $MERGED_ASSOC_DATA .tsv) FALSE '\$LSB_JOBINDEX'

# bsub -J "bart-predict-$(basename $MERGED_DIRECT_DATA .tsv)[1-${NUM}]" \
#      -o ${LOGDIR}/bart-predict-$(basename $MERGED_DIRECT_DATA .tsv)-%I.out \
#      -e ${LOGDIR}/bart-predict-$(basename $MERGED_DIRECT_DATA .tsv)-%I.err \
#      -M $MEM  -R "rusage[mem=12000]" -n $CORES \
#      Rscript src/bart-predict.r $MERGED_DIRECT_DATA out/$(basename $MERGED_DIRECT_DATA .tsv) TRUE '\$LSB_JOBINDEX'
