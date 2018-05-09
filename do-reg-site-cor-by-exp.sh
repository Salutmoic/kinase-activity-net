#!/bin/bash

LOGDIR=log
MEM=8000

mkdir -p $LOGDIR
mkdir -p out/reg-site-cor-by-exp

for n in $(seq 1 20); do
         bsub -J reg-site-cor-by-exp-${n} \
              -o ${LOGDIR}/reg-site-cor-by-exp-${n}.out \
              -e ${LOGDIR}/reg-site-cor-by-exp-${n}.err \
              -M $MEM  -R "rusage[mem=${MEM}]" \
              Rscript src/reg-site-cor-by-exp.r $n
done
