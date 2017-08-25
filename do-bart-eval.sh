#!/bin/bash

LOGDIR=log
MEM=$(echo "512*1024" | bc)
CORES=24

MERGED_ASSOC_DATA=out/kinase-association-feats.tsv
MERGED_DIRECT_DATA=out/kinase-direction-feats.tsv
MERGED_FULL_DATA=out/kinase-merged-feats-clean-hinegs.tsv

ASSOC_VALSET=data/assoc-validation-set-omnipath.tsv
DIRECT_VALSET=data/direct-validation-set-omnipath.tsv
NEG_VALSET=data/validation-set-negative.tsv

mkdir -p $LOGDIR

# bsub -J bart-eval-$(basename $MERGED_ASSOC_DATA .tsv) \
#      -o ${LOGDIR}/bart-eval-$(basename $MERGED_ASSOC_DATA .tsv).out \
#      -e ${LOGDIR}/bart-eval-$(basename $MERGED_ASSOC_DATA .tsv).err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-eval.r $MERGED_ASSOC_DATA $ASSOC_VALSET $NEG_VALSET TRUE FALSE

# bsub -J bart-eval-$(basename $MERGED_DIRECT_DATA .tsv) \
#      -o ${LOGDIR}/bart-eval-$(basename $MERGED_DIRECT_DATA .tsv).out \
#      -e ${LOGDIR}/bart-eval-$(basename $MERGED_DIRECT_DATA .tsv).err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-eval.r $MERGED_DIRECT_DATA $DIRECT_VALSET $NEG_VALSET FALSE TRUE

bsub -J bart-eval-$(basename $MERGED_FULL_DATA .tsv) \
     -o ${LOGDIR}/bart-eval-$(basename $MERGED_FULL_DATA .tsv).out \
     -e ${LOGDIR}/bart-eval-$(basename $MERGED_FULL_DATA .tsv).err \
     -M $MEM -P bigmem  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/bart-eval.r $MERGED_FULL_DATA $DIRECT_VALSET $NEG_VALSET TRUE FALSE
