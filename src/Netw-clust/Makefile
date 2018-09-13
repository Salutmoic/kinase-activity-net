
#This pipeline should generate the cluster files, enrichments and relevant plots

TABLE = merged-predictor-hinegs-annot.tsv 
IDS = Entrez-ids.txt

plot: make-undir-netw.r clust.r enrich-clusts.r plot-clust.r $(TABLE)
	
	Rscript make-undir-netw.r $(TABLE)

	Rscript clust.r undir-netw/symmetric-netw.tsv 0.5 Louvain T F
	
	Rscript clust.r undir-netw/symmetric-netw.tsv 0.5 Louvain T T Separated_clusters/scaled-Louvain-0.5-1.tsv
	Rscript clust.r undir-netw/symmetric-netw.tsv 0.5 Louvain T T Separated_clusters/scaled-Louvain-0.5-2.tsv
	Rscript clust.r undir-netw/symmetric-netw.tsv 0.5 Louvain T T Separated_clusters/scaled-Louvain-0.5-3.tsv
	Rscript clust.r undir-netw/symmetric-netw.tsv 0.5 Louvain T T Separated_clusters/scaled-Louvain-0.5-4.tsv
	Rscript clust.r undir-netw/symmetric-netw.tsv 0.5 Louvain T T Separated_clusters/scaled-Louvain-0.5-5.tsv	
	
	Rscript enrich-clusts.r louvain $(TABLE) $(IDS)
	 
	Rscript plot-clust.r clusters/scaled-Louvain-0.5.tsv $(IDS)
	
