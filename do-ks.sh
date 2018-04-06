#!/bin/bash

NUM_CONDS=588
LOGDIR=log/ks
USE_AUTOPHOS=TRUE
TMPDIR=tmp/ks-autophos

mkdir -p $LOGDIR
mkdir -p $TMPDIR

bsub -M 8192 -R "rusage[mem=1024]" -o ${LOGDIR}/ks-%I.out -e ${LOGDIR}/ks-%I.err \
     -J "ks[1-${NUM_CONDS}]" Rscript src/kinact-ks.r '\$LSB_JOBINDEX' $USE_AUTOPHOS
