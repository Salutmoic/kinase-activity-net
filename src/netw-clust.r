library(igraph)
library(reshape2)
library("ReactomePA")
library("clusterProfiler")
library("UpSetR")

#same as pathway.r in general but uses different clustering methods

convert = function(gene.list,id.map){

	entrez.list = id.map[which(id.map[,1] %in% gene.list),2]

}

make_graph = function(network,type){

	if(type == "adj"){
		adj = acast(network,network[,1]~network[,2])	
		graph = graph_from_adjacency_matrix(adj,mode = "undirected",weighted = T)
		weights = E(graph)$weight
	}else{
		netw = network[which(as.numeric(network[,3]) != 0),]
		scale.w = scale(as.numeric(netw[,3]))
		graph  = graph_from_edgelist(as.matrix(netw[,1:2]), directed = F)
		weights = as.numeric(scale.w)
	}
	return(list(graph,weights))
}

cluster = function(graph, method,w){
	clusters = list()
	if(method == "greedy"){
		community = cluster_fast_greedy(graph,weights = w)

	}
	if(method == "Louvain"){
		community = cluster_louvain(graph,weights = w)	
	}
	
	cnames =  as.numeric(which(table(membership(community))>= 10))

	for(i in 1:length(cnames)){
		clusters[[as.character(cnames[i])]] = names(which(membership(community)==cnames[i]))
	}
	
	return(clusters) 
}

enrich = function(clusters,id.map,entrez,method){
	
	entrez.list = list()
	pathways = list()
	cluster.enrichment = list()
	for(name in names(clusters)){
		genes = clusters[[name]]
		entrez.genes = convert(genes,id.map)
		print(name)
		x = enrichPathway(as.character(entrez.genes),pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(entrez))
		
		if(is.null(x)){
		
			cluster.enrichment[[name]] = NULL
			pathways[[name]] = NULL
			next
		}
		print(x)
		dir.create(paste("img3/",method,"/",sep = ""))
		dir.create(paste("img3/",method,"/",threshold,"/",sep = ""))
		pdf(paste("img3/",method,"/",threshold,"/plot-final-",name,"-",threshold,"-",length(entrez.genes),"-2.pdf",sep =""),width = 15,height = 8)
		g  =barplot(x,showCategory=10)
		print(g)
		dev.off()
		
		pathways[[name]] = x$Description
		entrez.list[[name]] = entrez.genes
		
	}
	return(list(pathways, entrez.list))
}

make.set.mat = function(clusts,paths){

	set.matrix = matrix(0,ncol = length(clusts),nrow = length(paths))
	rownames(set.matrix) = paths
	colnames(set.matrix) = names(clusts)

	for(clust in names(clusts)){
		set.matrix[which(paths %in% clusts[[clust]]),clust] = 1
	}
	return(data.frame(set.matrix))
}

scaling = function(values){

	s =(values-min(values))/(max(values) - min(values))
	return(s)

}

argv <- commandArgs(TRUE)

netfile = argv[1]
network = read.delim(netfile)
network = network[,1:3]
threshold = as.numeric(argv[2])
method = argv[3]
graph.type = argv[4]

network[which(is.na(network[,3])),3] =0
network[which(network[,3] <threshold),3] = 0
network = network[,1:3]

id.map = read.delim("Entrez-ids.txt")
id.map = id.map[which(id.map$Species == "Homo sapiens"),]
net.sub = network[which(as.numeric(network[,3]) > threshold),]
genes = unique(c(as.character(net.sub[,1]),as.character(net.sub[,2])))
entrez = convert(genes,id.map)

graph = make_graph(network,graph.type)

weights = graph[[2]]
graph = graph[[1]]
clusters = cluster(graph, method,weights)

save(clusters,file =paste("img3/clusters/",threshold,"-",method,"-",graph.type,".Rdata",sep = ""))

pathways = enrich(clusters,id.map,entrez,method)
entrez.list = pathways[[2]]
pathways = pathways[[1]]
set.matrix = make.set.mat(clusters,unique(unlist(pathways)))

dir.create(paste("img3/upset/",sep = ""))
pdf(paste("img3/upset/set-test--final-plot-final-",method,"-",threshold,"-2.pfd",sep = ""),width = 15,height = 8)
upset(set.matrix,sets.bar.color = "#56B4E9",order.by = "freq")
dev.off()

dir.create(paste("img3/clust-com/",sep = ""))
pdf(paste("img3/clust-com/cluster-test-adj-comparison-",method,"-",threshold,".pdf",sep = ""))
cluster = entrez.list
names(cluster) = as.character(seq(1,length(cluster)))
res <- compareCluster(cluster, fun="enrichPathway",pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(entrez))
g = dotplot(res,font.size = 6)
print(g)
dev.off()








