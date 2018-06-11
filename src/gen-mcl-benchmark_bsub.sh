#!/bin/bash

COR_DIR=out
COR_TBL=kinase-associations-pruned

LOG_DIR=log/${COR_TBL}
COS=$(seq 0.0 0.02 0.98)

mkdir -p ${LOG_DIR}

for co in $COS; do
    bsub -q research-rh7 -o ${LOG_DIR}/mcl-bench-${COR_TBL}-${co}.out \
        -e ${LOG_DIR}/mcl-bench-${COR_TBL}-${co}.err \
        -n 6 \
        src/gen-mcl-benchmark_subtask.sh "${COR_DIR}" "${COR_TBL}" "${co}"
done

