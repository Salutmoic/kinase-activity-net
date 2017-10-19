#!/bin/bash

TABLE_STRATEGIES="max-rows max-cols balanced"
ASSOC_METHODS="pcor pcor-filter scor scor-filter nfchisq mut_info mut_info-filter fnn_mut_info fnn_mut_info-filter partcor all paircor fvalue"
LOGDIR=log

mkdir -p $LOGDIR

for s in $TABLE_STRATEGIES; do
    for m in $ASSOC_METHODS;  do
        make clean-validation TABLE_STRATEGY=$s ASSOC_METHOD=$m
        bsub -o ${LOGDIR}/validation-${s}-${m}.out -e ${LOGDIR}/validation-${s}-${m}.err \
             make img/kinact-${s}-${m}-val.pdf TABLE_STRATEGY=$s ASSOC_METHOD=$m
    done
    bsub -o ${LOGDIR}/validation-${s}-pssm.out -e ${LOGDIR}/validation-${s}-pssm.err \
         make img/kinact-${s}-pssm-val.pdf TABLE_STRATEGY=$s
    bsub -o ${LOGDIR}/validation-${s}-string-coexp.out -e ${LOGDIR}/validation-${s}-string-coexp.err \
         make img/kinact-${s}-string-coexp-val.pdf TABLE_STRATEGY=$s
    bsub -o ${LOGDIR}/validation-${s}-string-exper.out -e ${LOGDIR}/validation-${s}-string-exper.err \
         make img/kinact-${s}-string-exper-val.pdf TABLE_STRATEGY=$s
done
