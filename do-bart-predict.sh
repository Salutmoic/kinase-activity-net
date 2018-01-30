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

# bsub -J bart-predict \
#      -o ${LOGDIR}/bart-predict.out -e ${LOGDIR}/bart-predict.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-predict.r $ASSOC_DATA FALSE
# bsub -J bart-predict-rand-negs -o ${LOGDIR}/bart-predict.out -e ${LOGDIR}/bart-predict.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-predict.r $ASSOC_DATA TRUE

BART_MODELS="bart directed-bart rand-negs-bart"

for m in $BART_MODELS; do
         bsub -J bart-predict-full-${m} \
              -o ${LOGDIR}/bart-predict-full-${m}.out -e ${LOGDIR}/bart-predict-full-${m}.err \
              -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
              Rscript src/bart-predict.r $FULL_DATA out/kinact-merged_full-${m}.Rdata
done

# bsub -J bart-predict-direct \
#      -o ${LOGDIR}/bart-predict-direct.out -e ${LOGDIR}/bart-predict-direct.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-predict.r $DIRECT_DATA FALSE
# bsub -J bart-predict-direct-rand-negs \
#      -o ${LOGDIR}/bart-predict-direct-rand-negs.out -e ${LOGDIR}/bart-predict-direct-rand-negs.err \
#      -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-predict.r $DIRECT_DATA TRUE
