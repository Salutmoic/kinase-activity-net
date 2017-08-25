declare -a meth=("Louvain" "greedy") 

ts=$(seq 0.6 0.05 0.95)


for m in "${meth[@]}"
do
echo $m
	
	for t in $ts
		
	do
	echo $t
	Rscript clust.r full-annotated-predictor-pruned.tsv $t $m F
	Rscript clust.r full-annotated-predictor-pruned.tsv $t $m T
	
	Rscript clust.r full-annotated-predictor.tsv $t $m F
	Rscript clust.r full-annotated-predictor.tsv $t $m T

	Rscript clust.r full-annotated-predictor-to-cluster.tsv $t $m F
	Rscript clust.r full-annotated-predictor-to-cluster.tsv $t $m T
	done
	
done
		
