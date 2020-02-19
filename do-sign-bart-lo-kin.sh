#!/bin/bash

LOGDIR=log/sign-bart-lo-kin
MEM=$(echo "256*1024" | bc)
CORES=24

FULL_DATA=out/sign-feats.tsv

VALSET=data/sign-validation-set-omnipath.tsv

mkdir -p $LOGDIR

# bsub -J "bart-lo-kin-$(basename $MERGED_FULL_DATA .tsv)[1-${NUM}]" \
#      -o ${LOGDIR}/bart-lo-kin-$(basename $MERGED_FULL_DATA .tsv)-%I.out \
#      -e ${LOGDIR}/bart-lo-kin-$(basename $MERGED_FULL_DATA .tsv)-%I.err \
#      -M $MEM -P bigmem -R "rusage[mem=${MEM}]" -n $CORES \
#      Rscript src/bart-lo-kin.r $MERGED_FULL_DATA $DIRECT_VALSET '\$LSB_JOBINDEX'

# for kin in EGFR MTOR PDPK1 AKT1 GSK3B RPS6KB1 BRAF MAP2K1 MAPK1 AKT2 AKT3 RPS6KB2 RPS6KA3 RPS6KA4 GSK3A ARAF RAF1 MAP2K2 MAP2K3 MAPK3; do
for kin in AKT1 RPS6KB1 MAP2K1; do
bsub -J "bart-lo-kin-$(basename $FULL_DATA .tsv)-${n}" \
     -o ${LOGDIR}/sign-bart-lo-kin-$(basename $FULL_DATA .tsv)-${kin}.out \
     -e ${LOGDIR}/sign-bart-lo-kin-$(basename $FULL_DATA .tsv)-${kin}.err \
     -M $MEM -P bigmem -R "rusage[mem=${MEM}]" -n $CORES \
     Rscript src/sign-bart-lo-kin.r $FULL_DATA $VALSET $kin
done

