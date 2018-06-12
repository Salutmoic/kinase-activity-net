declare -a meth=("Louvain" "greedy") 
declare -a types=("adj" "other")
ts=$(seq 0.6 0.05 0.98)


for m in "${meth[@]}"
do
echo $m
	for ty in "$types"
	do
		echo $ty
		for t in $ts
		
		do
		echo $t
		Rscript netw-clust.r full-annotated-predictor.tsv $t $m $ty  
		done
	done
done
		
