#!/bin/bash

TABLE_STRATEGIES="max-rows max-cols balanced"
LOGDIR=log

mkdir -p $LOGDIR

for s in $TABLE_STRATEGIES; do
    make clean-data TABLE_STRATEGY=$s
    bsub -o ${LOGDIR}/gen-table-${s}.out -e ${LOGDIR}/gen-table-${s}.err make data \
         TABLE_STRATEGY=$s KSEA_MIN_SITES=5 ENTROPY_FILTER=0.5 NA_THRESHOLD=0.3333
done
