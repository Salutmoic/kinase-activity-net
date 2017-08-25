#!/bin/bash

KINACT=data/ks-kinact.Rdata

LOGDIR=log
MEM=8000

mkdir -p $LOGDIR
mkdir -p out/kinact-full-cor-by-exp

for n in $(seq 1 18); do
         bsub -J kinact-full-cor-by-exp-${n} \
              -o ${LOGDIR}/kinact-full-cor-by-exp-${n}.out \
              -e ${LOGDIR}/kinact-full-cor-by-exp-${n}.err \
              -M $MEM  -R "rusage[mem=${MEM}]" \
              Rscript src/kinact-full-cor-by-exp.r $KINACT $n
done
