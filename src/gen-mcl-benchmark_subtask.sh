#!/bin/bash

COR_TSV=${1}/${2}.tsv
BENCH_DIR=${1}/benchmark-${2}-clust
MSI_DIR=${BENCH_DIR}/msi
OUT_DIR=${BENCH_DIR}/out
INFLATIONS=$(seq 1.2 0.2 10.0)
SCALE=$(echo "1.0/(1.0-$3)" | bc -l)

mkdir -p ${MSI_DIR}
mkdir -p ${OUT_DIR}

MSI_FILE=${MSI_DIR}/$(basename $COR_TSV .tsv)_co${3}.msi
TAB_FILE=${MSI_DIR}/$(basename $COR_TSV .tsv)_co${3}.tab
[[ -r ${MSI_FILE} ]] || mcxload -abc $COR_TSV \
    --stream-mirror \
    --write-binary \
    -write-tab $TAB_FILE \
    -o $MSI_FILE \
    -tf "gt(${3}),add(-${3}),mul(${SCALE})"
CO_OUT_DIR=${OUT_DIR}/co${3}
mkdir -p ${CO_OUT_DIR}
for i in $INFLATIONS; do
    OUT_FILE=${CO_OUT_DIR}/out.$(basename ${MSI_FILE}).I${i//[.]/}
    [[ -r ${OUT_FILE} ]] || mcl $MSI_FILE -odir $CO_OUT_DIR -te 6 -I $i
done
best=$(clm info $MSI_FILE ${CO_OUT_DIR}/out.$(basename ${MSI_FILE}).I* | \
    cut -f1,5 -d' ' | \
    sed '/^===/d;s/efficiency=//g;s/source=.*\(I[0-9][0-9]*\)/\1/' | \
    sort -g -k1 | \
    sed -n '${s/0\.[0-9][0-9]* //;p}')
BEST_INF=out.$(basename ${MSI_FILE}).${best}
cp ${CO_OUT_DIR}/${BEST_INF} ${BENCH_DIR}
CLUST_FILE=${BENCH_DIR}/${BEST_INF}-clusts.out
[[ -r $CLUST_FILE ]] || mcxdump -icl ${BENCH_DIR}/${BEST_INF} \
    -tabr $TAB_FILE \
    -o $CLUST_FILE
CLUST_TAB=${BENCH_DIR}/${BEST_INF}-clusts.tsv
[[ -r $CLUST_TAB ]] || awk -F'	' '{for(i=1;i<=NF;i++){printf("%s\t%d\n", $i, NR)}}' \
    $CLUST_FILE >$CLUST_TAB
