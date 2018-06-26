cutoffs = c(seq(0,0.96, by = 0.02),seq(0.6,0.95,0.05))
cutoffs = cutoffs[order(cutoffs)]

network = read.delim("full-annotated-predictor.tsv")
network = network[,1:3]
network.prune = read.delim("full-annotated-predictor-pruned.tsv")
network.prune = network.prune[,1:3]
network.clust = read.delim("full-annotated-predictor-to-cluster.tsv")
network.clust = network.clust[,1:3]

network.prun.clust  = read.delim("assoc-net-to-cluster-pruned.tsv")
gene.nos = list()

prefixes = c("pruned","pruned-scaled","to-cluster","to-cluster-scaled","scaled","Louvain","no-prefix","pruned-to-cluster","to-cluster-pruned")

for(cutoff in cutoffs){
	
	network = network[which(network[,3] >= cutoff),]
	no.genes = length(unique(c(network[,1],network[,2])))

	gene.nos[["no-prefix"]][[as.character(cutoff)]] = no.genes
	gene.nos[["Louvain"]][[as.character(cutoff)]] = no.genes
	gene.nos[["scaled"]][[as.character(cutoff)]] = no.genes

	network.prune = network.prune[which(network.prune[,3] >= cutoff),]
	no.genes = length(unique(c(network.prune[,1],network.prune[,2])))

	gene.nos[["pruned"]][[as.character(cutoff)]] =no.genes
	gene.nos[["pruned-scaled"]][[as.character(cutoff)]] =no.genes

	network.clust = network.clust[which(network.clust[,3] >= cutoff),]
	no.genes = length(unique(c(network.clust[,1],network.clust[,2])))

	gene.nos[["to-cluster"]][[as.character(cutoff)]] = no.genes
	gene.nos[["to-cluster-scaled"]][[as.character(cutoff)]] = no.genes

	network.prun.clust = network.prun.clust[which(network.prun.clust[,3] >= cutoff),]
	no.genes = length(unique(c(network.prun.clust[,1],network.prun.clust[,2])))
	
	gene.nos[["to-cluster-pruned"]][[as.character(cutoff)]]=no.genes
	gene.nos[["pruned-to-cluster"]][[as.character(cutoff)]]=no.genes
	gene.nos[["pruned-to-cluster-scaled"]][[as.character(cutoff)]]=no.genes

}


save(gene.nos,file = "number-of-genes-per-cutoff.RData")



