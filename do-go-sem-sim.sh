#!/bin/bash

LOGDIR=log/go-sem-sim2
MEM=$(echo "2*1024" | bc)
CORES=1
NUM=$(wc -l data/human-kinome.txt | cut -f1 -d' ')

mkdir -p $LOGDIR
mkdir -p out/go-sem-sim2/

bsub -J "go-sem-sim[1-${NUM}]" \
     -o ${LOGDIR}/go-sem-sim-%I.out \
     -e ${LOGDIR}/go-sem-sim-%I.err \
     -M $MEM  -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/go-sem-sim.r '\$LSB_JOBINDEX'
