library(dplyr)
library(ggplot2)

filenames = c("mclpathway-enrichments.tsv","louvainpathway-enrichments.tsv","greedypathway-enrichments.tsv")

load("number-of-genes-per-cutoff.RData")

for(filename in filenames){ 
	m = read.delim(filename)

	files = unique(m[,1])
	#this table collect the various metrics for the different cutoffs and networks (median coverage/cluster ........)
	tbl.stats = c()

	if(length(grep("louvain",filename)) > 0|length(grep("greedy",filename)) ){

		for(f in files){
			sub.tbl = m[m[,1] ==f ,]

			substrs = strsplit(f,"-")
			cutoff = substrs[[1]][length(substrs[[1]])]
			cutoff = strsplit(cutoff,"\\.tsv")[[1]][1]

			ind =grep("Louvain",substrs[[1]])
			if(length(ind)>0){
				prefix = paste(substrs[[1]][1:(ind-1)],collapse = "-")
			}else{
				ind =grep("greedy",substrs[[1]])
				if(ind>1){
					prefix = paste(substrs[[1]][1:(ind-1)],collapse = "-")
				}else{
					prefix = "no-prefix"
				}
				print(prefix)
			}

			grp = group_by(sub.tbl, cluster)
			#median coverage score accross clusters with enriched patways.
			median = median(unlist(summarise(grp, median=no_genes_in_top_enriched[1]/no_genes_in_cluster[1])[,2]))
			mean = mean(unlist(summarise(grp, median=no_genes_in_top_enriched[1]/no_genes_in_cluster[1])[,2]))
			#sum of coverage score across different clusters			
			sum = sum(unlist(summarise(grp, sum=no_genes_in_top_enriched[1]/no_genes_in_cluster[1])[,2]))
			# weighted mean coverage score sum(genes in the top patway across clusters)/(number of genes found within these clusters) 
			w.mean= sum(unlist(summarise(grp, sum=no_genes_in_top_enriched[1])[,2]))/sum(unlist(summarise(grp, sum=no_genes_in_cluster[1])[,2]))
			# sum of genes found within clusters with enriched pathways			
			sum.genes = sum(unlist(summarise(grp, sum=sum(no_genes_in_cluster[1]))[,2]))
			# median number of genes/(no. cluster with enriched pathways)			
			median.genes = median(unlist(summarise(grp, median=median(no_genes_in_top_enriched[1]))[,2]))
	
			#number of genes in a network at a given cutoff	
			ngen = gene.nos[[prefix]][[as.character(as.numeric(cutoff))]]

			tbl.stats = rbind(tbl.stats,c(f,prefix,cutoff,nrow(sub.tbl),sub.tbl[1,4],length(unique(sub.tbl[,3])),length(unique(sub.tbl[,2])),
			mean,sum/sub.tbl$no_clusters[1],w.mean,sum.genes, (w.mean*sum.genes)/ngen,sqrt(mean*median.genes)))

			print(gene.nos[[prefix]][[as.character(as.numeric(cutoff))]])
		}
	}else{

		for(f in files){
			sub.tbl = m[m[,1] ==f ,]

			substrs = strsplit(f,"_co")
			cutoff = strsplit(substrs[[1]][2],".\\msi")[[1]][1]
			cutoff = strsplit(cutoff,"\\.tsv")[[1]][1]


			prefix = strsplit(substrs[[1]][1],"-net-")[[1]][2]
			if(is.na(prefix)){
				prefix = "no-prefix"
			}
			#same variables as above
			grp = group_by(sub.tbl, cluster)
			median = median(unlist(summarise(grp, median=median(no_genes_in_top_enriched[1]/no_genes_in_cluster[1]))[,2]))
			mean = 	mean(unlist(summarise(grp, median=median(no_genes_in_top_enriched[1]/no_genes_in_cluster[1]))[,2]))		
			sum = sum(unlist(summarise(grp, sum=no_genes_in_top_enriched[1]/no_genes_in_cluster[1])[,2]))

			w.mean= sum(unlist(summarise(grp, sum=sum(no_genes_in_top_enriched[1]))))/sum(unlist(summarise(grp, sum=sum(no_genes_in_cluster[1]))))
			sum.genes = sum(unlist(summarise(grp, sum=sum(no_genes_in_cluster[1]))))
			median.genes = median(unlist(summarise(grp, median=median(no_genes_in_top_enriched[1]))[,2]))

			ngen = gene.nos[[prefix]][[as.character(as.numeric(cutoff))]]
			
			tbl.stats = rbind(tbl.stats,c(f,prefix,cutoff,nrow(sub.tbl),sub.tbl[1,4],length(unique(sub.tbl[,3])),
			length(unique(sub.tbl[,2])),mean,sum/sub.tbl$no_clusters[1],
			w.mean,sum.genes,(w.mean*sum.genes)/ngen,sqrt(mean*median.genes)))
			print(prefix)
			#print(gene.nos[[prefix]][[as.character(as.numeric(cutoff))]])
		}


	}

	colnames(tbl.stats) = c("file","prefix","cutoff","no.enrichments","no.clusts"
	,"no.paths","no.enrich.clusts","genes.enrich.clust","genes.all.clusts","w.mean.cov","sum.genes","comb.score","geom_mean")



	tbl.stats = tbl.stats[order(-as.numeric(tbl.stats[,4]),-as.numeric(tbl.stats[,5]),-as.numeric(tbl.stats[,3])),]

	tbl.stats = data.frame(tbl.stats)

	for(i in 3:ncol(tbl.stats)){
		tbl.stats[,i] = as.numeric(as.character(tbl.stats[,i]))

	}	




	if(length(grep("louvain",filename)) > 0){

		#tbl.stats = tbl.stats[tbl.stats$prefix == "pruned-to-cluster-scaled",]
		#tbl.stats = tbl.stats[tbl.stats$prefix == "scaled",]
		tbl.lov = tbl.stats[grep("Louvain",tbl.stats[,1]),]


		tbl.g = tbl.stats[grep("greedy",tbl.stats[,1]),]

		#pdf("louvain-unique-pathways.pdf")
		#g=ggplot(tbl.lov, aes(x = cutoff, y = no.paths,color = prefix)) + geom_point()
		#print(g)
		#dev.off()

		#pdf("louvain-enrich-clusters.pdf")
		#g= ggplot(tbl.lov, aes(x = cutoff, y = no.enrich.clusts,color = prefix)) + geom_point()
		#print(g)
		#dev.off()

		#pdf("louvain-enrich-clusters-vs-depleted-clusters.pdf")
		#g= ggplot(tbl.lov, aes(x = cutoff, y = no.enrich.clusts/no.clusts,color = prefix)) + geom_point()
		#print(g)
		#dev.off()


		pdf("prunde-to-clust-louvain-geom-mean.pdf")
		g= ggplot(tbl.lov, aes(x = cutoff, y = geom_mean)) + geom_point()
		print(g)
		dev.off()

		pdf("prunde-to-clust-greedy-geom-mean.pdf")
		g= ggplot(tbl.g, aes(x = cutoff, y = geom_mean)) + geom_point()
		print(g)
		dev.off()

		pdf("prunde-to-clust-louvain-median-coverage-per-enrich-clusters.pdf")
		g=ggplot(tbl.lov, aes(x = cutoff, y = genes.enrich.clust,color = prefix)) + geom_point()
		print(g)
		dev.off()


		pdf("prunde-to-clust-louvain-mean-coverage-per-all-clusters.pdf")
		g=ggplot(tbl.lov, aes(x = cutoff, y = genes.all.clusts,color = prefix)) + geom_point()
		print(g)
		dev.off()



		pdf("prunde-to-clust-greedy-median-coverage-per-enrich-clusters.pdf")
		g=ggplot(tbl.g, aes(x = cutoff, y = genes.enrich.clust,color = prefix)) + geom_point()
		print(g)
		dev.off()


		pdf("prunde-to-clust-greedy-mean-coverage-per-all-clusters.pdf")
		g= ggplot(tbl.g, aes(x = cutoff, y = genes.all.clusts,color = prefix)) + geom_point()
		print(g)
		dev.off()

		pdf("prunde-to-clust-greedy-weighted-coverage-all-clusters.pdf")
		g=ggplot(tbl.g, aes(x = cutoff, y = w.mean.cov,color = prefix)) + geom_point()
		print(g)
		dev.off()

		#pdf("greedy-sum-genes-all-clusters.pdf")
		#g=ggplot(tbl.g, aes(x = cutoff, y = sum.genes,color = prefix)) + geom_point()
		#print(g)
		#dev.off()


		pdf("prunde-to-clust-lov-weighted-coverage-all-clusters.pdf")
		g=ggplot(tbl.stats, aes(x = cutoff, y = w.mean.cov,color = prefix)) + geom_point()
		print(g)
		dev.off()

		#pdf("lov-sum-genes-all-clusters.pdf")
		#g=ggplot(tbl.stats, aes(x = cutoff, y = sum.genes,color = prefix)) + geom_point()
		#print(g)
		#dev.off()

		pdf("prunde-to-clust-lov-com-score.pdf")
		g=ggplot(tbl.lov, aes(x = cutoff, y = comb.score,color = prefix)) + geom_point()
		print(g)
		dev.off()


		#pdf("greedy-sum-genes-all-clusters.pdf")
		#g=ggplot(tbl.g, aes(x = cutoff, y = sum.genes,color = prefix)) + geom_point()
		#print(g)
		#dev.off()

		pdf("prunde-to-clust-greedy-com-score.pdf")
		g=ggplot(tbl.g, aes(x = cutoff, y = comb.score,color = prefix)) + geom_point()
		print(g)
		dev.off()


	}else{
		print(tbl.stats)
		tbl.stats$prefix = gsub("-","",tbl.stats$prefix)
		#tbl.stats = tbl.stats[tbl.stats$prefix == "noprefix",]
		tbl.mcl = tbl.stats
		#tbl.stats = tbl.stats[tbl.stats$prefix == "to-cluster-pruned",]
		

		pdf("prunde-to-clust-mcl-medina-coverage-per-enrich-clusters.pdf")
		g= ggplot(tbl.stats, aes(x = cutoff, y = genes.enrich.clust,color = prefix)) + geom_point()
		print(g)
		dev.off()


		pdf("prunde-to-clust-mcl-mean-coverage-per-all-clusters.pdf")
		g=ggplot(tbl.stats, aes(x = cutoff, y = genes.all.clusts,color = prefix)) + geom_point()
		print(g)
		dev.off()

		pdf("prunde-to-clust-mcl-weighted-coverage-all-clusters.pdf")
		g=ggplot(tbl.stats, aes(x = cutoff, y = w.mean.cov,color = prefix)) + geom_point()
		print(g)
		dev.off()

		#pdf("mcl-sum-genes-all-clusters.pdf")
		#g=ggplot(tbl.stats, aes(x = cutoff, y = sum.genes,color = prefix)) + geom_point()
		#print(g)
		#dev.off()

		pdf("prunde-to-clust-mcl-comb-score.pdf")
		g=ggplot(tbl.stats, aes(x = cutoff, y = comb.score,color = prefix)) + geom_point()
		print(g)
		dev.off()

		pdf("prunde-to-clust-mcl-geom-mean.pdf")
		g=ggplot(tbl.stats, aes(x = cutoff, y = geom_mean,color = prefix)) + geom_point()
		print(g)
		dev.off()

	}
}
