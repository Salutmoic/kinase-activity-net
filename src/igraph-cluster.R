library(igraph)
library(reshape2)
library("ReactomePA")
library('MCL')
library('VennDiagram')
library("clusterProfiler")

#same as pathway.r in general but uses different clustering methods

convert = function(gene.list,id.map){

entrez.list = id.map[which(id.map[,1] %in% gene.list),2]

}

Venn= function(a,b){
	c = length(intersect(a,b))
	a = length(a)
	b = length(b)

	draw.pairwise.venn(a,b,c,fill= c("green","orange"))
}

id.map = read.delim("Entrez-ids.txt")
id.map = id.map[which(id.map$Species == "Homo sapiens"),]

argv <- commandArgs(TRUE)

netfile = argv[1]
network = read.delim(netfile)

directed = as.logical(argv[2])
random = as.logical(argv[3])
threshold = as.numeric(argv[4])

print(threshold)

network[which(is.na(network[,3])),3] =0
network[which(network[,3] <threshold),3] = 0
network = network[,1:3]
adj = acast(network,network[,1]~network[,2])
network = network[which(network[,3] >=threshold),]



dir = ""
if(random == T){
	#do random decided to do that elsewhere
}


if(directed == T){
	network =network[,c("prot1","prot2","direct.pred.mean")] 
	directed = "directed"
}else{
	
	directed = "undirected"
}


genes = unique(c(as.character(network[,1]),as.character(network[,2])))
entrez = convert(genes,id.map)

graph  = graph_from_edgelist(as.matrix(network[,1:2]), directed = F)
graph2 = graph_from_adjacency_matrix(adj,mode = "undirected",weighted = T)

greedy.community = cluster_fast_greedy(graph2,weights = E(graph2)$weight )
louvain.community = cluster_louvain(graph,weights = network[,3] )


write.table(t(c("louvain",as.character(random),as.character(modularity(louvain.community)))),"img2/cluster-final-modularity.tsv",quote =F,sep = "\t",col.names = F,append = T)
write.table(t(c("greedy",as.character(random),as.character(modularity(greedy.community)))),"img2/cluster-final-modularity.tsv",quote =F,sep = "\t",col.names = F,append = T)

lov.names =  as.numeric(which(table(membership(louvain.community))>= 10))
greedy.names =  as.numeric(which(table(membership(greedy.community))>= 10))

names.vect = c()
coms = list()

for(i in 1:(length(lov.names)+length(greedy.names))){

	index = i
	if(i >length(lov.names)){
		index = i- length(lov.names)
		coms[[i]] = names(which(membership(greedy.community)==greedy.names[index]))
		names.vect = c(names.vect,paste("greedy",as.character(greedy.names[index]),sep = "-"))
	}else{
		coms[[i]] = names(which(membership(louvain.community)==lov.names[index]))
		names.vect = c(names.vect,paste("lovain",as.character(lov.names[index]),sep = "-"))		
	}
	

}
names(coms) = names.vect

#use MCL

adj[which(adj != 0)]= 1
print(table(adj ==1))

k = mcl(adj,addLoops = T)
clusters = unique(k$Cluster)
print(clusters)
entrez.list = list()

for(i in 1:k$K){
	if(length(which(k$Cluster ==clusters[i])) >= 10){
		coms[[paste("mcl",as.character(i),sep = "-")]] = rownames(adj)[which(k$Cluster ==clusters[i])]
	}
}

cluster.enrichment = list()
pathways = list()
entrez.list = list()

for(name in names(coms)){
	genes = coms[[name]]
	
	entrez.genes = convert(genes,id.map)
	print(name)
	x = enrichPathway(as.character(entrez.genes),pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(entrez))
	method = strsplit(name,"-")[[1]][1]
	
	if(is.null(x)){
		
		cluster.enrichment[[name]] = NULL
		pathways[[name]] = NULL
		next
	}
	print(x)
	pdf(paste("img2/",dir,method,"/",threshold,"/plot-final-",name,"-",threshold,"-",directed,"-",length(entrez.genes),"-2.pdf",sep =""),width = 15,height = 8)
	g  =barplot(x,showCategory=10)
	print(g)
	dev.off()
	cluster.enrichment[[name]] = x
	pathways[[name]] = x$Description
	entrez.list[[name]] = entrez.genes
	
}


save(cluster.enrichment,file = paste("cluster-enrichment-",as.character(random),"-",threshold,"-",directed,".Rdata",sep = ""))


#mcode
#markov
#imap
#pdf("img2/lovain-venn.pdf")

#g =Venn(pathways[[1]],pathways[[2]])

#print(g)
#dev.off()

#pdf("img2/greedy-venn.pdf")

#g =Venn(pathways[[3]],pathways[[4]])

#print(g)
#dev.off()

#Visualize mclust which has 14 clusters harder to visualize
library("UpSetR")
paths= unique(unlist(pathways[1:length(pathways)]))
clusts = pathways

set.matrix = matrix(0,ncol = length(clusts),nrow = length(paths))

rownames(set.matrix) = paths
colnames(set.matrix) = names(clusts) 

for(clust in names(clusts)){

set.matrix[which(paths %in% clusts[[clust]]),clust] = 1
}

methods = c("lovain","mcl","greedy")

for(method in methods){
	pdf(paste("img2/",dir,"upset/set-final-plot-final-",method,"-",threshold,"-",directed,"-2.pfd",sep = ""),width = 15,height = 8)

	method.matrix = data.frame(set.matrix[,grep(method,colnames(set.matrix))])
	colsums = apply(method.matrix,2,sum)
	if(sum(colsums) == 0 | length(which(colsums > 0)) < 2){
		next
	}
	upset(method.matrix,sets.bar.color = "#56B4E9",order.by = "freq")

	dev.off()

	pdf(paste("img2/",dir,"clust-com/cluster-comparison-",method,"-",threshold,".pdf",sep = ""))
	cluster = entrez.list[grep(method,names(entrez.list))]
	names(cluster) = as.character(seq(1,length(cluster)))
	res <- compareCluster(cluster, fun="enrichPathway",pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(entrez))
	g = dotplot(res,font.size = 6)
	print(g)
	dev.off()

}















