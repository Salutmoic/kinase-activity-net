#!/bin/bash

TABLE_STRATEGIES="max-rows max-cols max-rows-max-cols"
ASSOC_METHODS="pcor pcor-filter scor scor-filter nfchisq mut_info fnn_mut_info partcor"
LOGDIR=log

mkdir -p $LOGDIR

for s in $TABLE_STRATEGIES; do
    for m in $ASSOC_METHODS;  do
        bsub -o ${LOGDIR}/validation-${s}-${m}.out -e ${LOGDIR}/validation-${s}-${m}.err \
             make img/kinact-${s}-${m}-val.pdf img/kinact-${s}-${m}-final-predictor-val.pdf
    done
    bsub -o ${LOGDIR}/validation-${s}-pssm.out -e ${LOGDIR}/validation-${s}-pssm.err \
         make img/kinact-${s}-pssm-val.pdf
done
