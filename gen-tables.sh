#!/bin/bash

TABLE_STRATEGIES="max-rows max-cols balanced"
LOGDIR=log

mkdir -p $LOGDIR

for s in $TABLE_STRATEGIES; do
    bsub -o ${LOGDIR}/gen-table-${s}.out -e ${LOGDIR}/gen-table-${s}.err make data TABLE_STRATEGY=$s
done
