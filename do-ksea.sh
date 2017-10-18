#!/bin/bash

NUM_CONDS=668
NUM_TRIALS=10000
LOGDIR=log/ksea
USE_AUTOPHOS=TRUE

mkdir -p $LOGDIR

bsub -M 32000 -n 10 -o ${LOGDIR}/ksea-%I.out -e ${LOGDIR}/ksea-%I.err \
     -J "ksea[1-${NUM_CONDS}]" Rscript src/ksea.r '\$LSB_JOBINDEX' $NUM_TRIALS $USE_AUTOPHOS
