library(igraph)
library(reshape2)
library("ReactomePA")
library("clusterProfiler")
library("UpSetR")

convert = function(gene.list,id.map){
	entrez.list = id.map[which(id.map[,1] %in% gene.list),2]
	return(entrez.list)
}

argv <- commandArgs(TRUE)
clusters = read.delim(argv[1])

id.map = read.delim(argv[2])
id.map = id.map[which(id.map$Species == "Homo sapiens"),]

clusts = names(table(clusters[,2]))
clust.kins = list()

for(clust in clusts){

kins = clusters[clusters[,2] == clust,1]
entrez = convert(kins,id.map)

clust.kins[[clust]] = entrez

}

kinome = read.delim("merged-predictor-hinegs-annot.tsv")
kinome = unique(c(as.character(kinome[,1]),as.character(kinome[,2])))

universe = convert(kinome,id.map)

png(paste("plots/",gsub(".tsv","",gsub("clusters/","",argv[1])),"-Louvain-enrich-results-lovain-0.55.png",sep = ""))
res <- compareCluster(clust.kins, fun="enrichPathway",pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(universe))
	g = dotplot(res,font.size = 8)
	print(g)
	dev.off()


paths.tbl = read.delim("enrichments/louvain-pathway-enrichments.tsv")
paths.tbl = paths.tbl[paths.tbl[,1] == gsub("clusters/","",argv[1]),]

clustr = paths.tbl[,2]
clustr = paste("Cluster_",clustr,sep = "")
paths = paths.tbl[,3]

set.matrix = matrix(0,ncol = length(unique(clustr)),nrow = length(paths))

rownames(set.matrix) = paths
colnames(set.matrix) = unique(clustr)

for(clust in unique(clustr)){

set.matrix[which(paths %in% paths[which(clustr==clust)]),clust] = 1
}

print(head(set.matrix))

png(paste("plots/",gsub("clusters/","",argv[1]),"-Louvain-upset-rev.png",sep = ""))

method.matrix = data.frame(set.matrix)
colsums = apply(method.matrix,2,sum)
	
upset(method.matrix,sets.bar.color = "#56B4E9",order.by = "freq")

dev.off()


