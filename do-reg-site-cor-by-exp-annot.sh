#!/bin/bash

LOGDIR=log
MEM=4000

mkdir -p $LOGDIR
mkdir -p out/reg-site-cor-by-exp

for n in $(seq 1 2); do
         bsub -J reg-site-cor-by-exp-annot-${n} \
              -o ${LOGDIR}/reg-site-cor-by-exp-annot-${n}.out \
              -e ${LOGDIR}/reg-site-cor-by-exp-annot-${n}.err \
              -M $MEM  -R "rusage[mem=${MEM}]" \
              Rscript src/reg-site-cor-by-exp-annot.r $n
done
