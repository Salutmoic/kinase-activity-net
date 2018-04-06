#!/bin/bash

LOGDIR=log/go-sem-sim2
MEM=$(echo "16*1024" | bc)
CORES=10

bsub -J go-bp-sem-sim -o log/go-bp-sem-sim.out -e log/go-cc-sem-sim.err \
     -M ${MEM} -R "rusage[mem=${MEM}]" -n $CORES \
     java -jar src/sml-toolkit-0.9.4c.jar -t sm -profile GO \
     -annots data/external/goa_human-2017-09-26.gaf \
     -go data/external/go-basic-2017-10-04.obo \
     -queries tmp/kinome-all-vs-all.tsv \
     -output out/go-bp-sem-sim-uniprot.tsv \
     -mtype g \
     -filter noEC=ISS,ISO,ISA,ISM,IGC,IBA,IBD,IKR,KRD,RCA,TAS,NAS,HTP,HDA,HMP,HGI,HEP \
     -notfound set=-1 \
     -pm schlicker -ic resnik -gm bma -aspect BP -quiet \
     -threads $CORES

bsub -J go-cc-sem-sim -o log/go-cc-sem-sim.out -e log/go-cc-sem-sim.err \
     -M ${MEM} -R "rusage[mem=${MEM}]" -n $CORES \
     java -jar src/sml-toolkit-0.9.4c.jar -t sm -profile GO \
     -annots data/external/goa_human-2017-09-26.gaf \
     -go data/external/go-basic-2017-10-04.obo \
     -queries tmp/kinome-all-vs-all.tsv \
     -output out/go-cc-sem-sim-uniprot.tsv \
     -mtype g \
     -filter noEC=ISS,ISO,ISA,ISM,IGC,IBA,IBD,IKR,KRD,RCA,TAS,NAS,HTP,HDA,HMP,HGI,HEP \
     -notfound set=-1 \
     -pm schlicker -ic resnik -gm bma -aspect CC -quiet \
     -threads $CORES

bsub -J go-mf-sem-sim -o log/go-mf-sem-sim.out -e log/go-mf-sem-sim.err \
     -M ${MEM} -R "rusage[mem=${MEM}]" -n $CORES \
     java -jar src/sml-toolkit-0.9.4c.jar -t sm -profile GO \
     -annots data/external/goa_human-2017-09-26.gaf \
     -go data/external/go-basic-2017-10-04.obo \
     -queries tmp/kinome-all-vs-all.tsv \
     -output out/go-mf-sem-sim-uniprot.tsv \
     -mtype g \
     -filter noEC=ISS,ISO,ISA,ISM,IGC,IBA,IBD,IKR,KRD,RCA,TAS,NAS,HTP,HDA,HMP,HGI,HEP \
     -notfound set=-1 \
     -pm schlicker -ic resnik -gm bma -aspect MF -quiet \
     -threads $CORES
