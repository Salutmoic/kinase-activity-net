
#This pipeline should generate the cluster files, enrichments and relevant plots

TABLE = merged-predictor-hinegs-annot.tsv 
IDS = Entrez-ids.txt

enrich: make-undir-netw.r clust.r enrich-clusts.r plot-clust.r $(TABLE)
	
	Rscript make-undir-netw.r $(TABLE)

	Rscript clust.r undir-netw/symmetric-netw.tsv 0.5 Louvain T F
		
	Rscript enrich-clusts.r louvain $(TABLE) $(IDS)
	 

	
