library("igraph")
library("ggplot2")
library("reshape2")
library("ReactomePA")
library("UpSetR")

gen.random = function(kin.pairs,values){

	values = c(values,rep(0,nrow(kin.pairs)-length(values)))
	rand.net = data.frame(kin.pairs[,1],kin.pairs[,2],sample(values))

	rand.net = rand.net[rand.net[,3] > 0,]
	return(rand.net)
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

convert = function(gene.list,id.map){

entrez.list = id.map[which(id.map[,1] %in% gene.list),2]

}

id.map = read.delim("Entrez-ids.txt")
id.map = id.map[which(id.map$Species == "Homo sapiens"),]

network = read.delim("/home/borgthor/full-annotated-predictor.tsv")
network[which(is.na(network[,3])),3] =0
network[which(network[,3] <= 0.8),3]= 0
network = network[,1:3]

modul = c()

for(i in 1:100){

	create random network and sub network from initial network
	net.sub = network[which(network[,3] > 0),]
	samp =sample(nrow(net.sub),ceiling(0.2*nrow(net.sub)))
	net.sub[samp,3] = 0
	net.sub = net.sub[net.sub[,3] > 0,]
	all.kins = unique(c(as.character(net.sub[,1]),as.character(net.sub[,2])))
	kin.pairs <- expand.grid(all.kins, all.kins, stringsAsFactors=FALSE)
	kin.pairs= kin.pairs[-which(kin.pairs[,1] ==kin.pairs[,2]),]
	net.rand = gen.random(kin.pairs,net.sub[,3])
	colnames(net.rand) = colnames(network)
	print(head(net.rand))
	#greate gene universe for each set
	genes = unique(c(as.character(net.sub[,1]),as.character(net.sub[,2])))
	entrez = convert(genes,id.map)

	genes.rand = unique(c(as.character(net.rand[,1]),as.character(net.rand[,2])))
	entre.rand = convert(genes.rand,id.map)
	
	#make graph
	graph  = graph_from_edgelist(as.matrix(net.sub[,1:2]), directed = F)
	
	graph.random  = graph_from_edgelist(as.matrix(net.rand[,1:2]), directed = F)

	#cluster
	louvain.community = cluster_louvain(graph,weights = net.sub[,3] )
	lov.names =  as.numeric(which(table(membership(louvain.community))>= 10))

	louvain.random = cluster_louvain(graph.random,weights = net.rand[,3] )
	rand.names =  as.numeric(which(table(membership(louvain.random))>= 10))
	modul = rbind(modul,c(i,"non-random",modularity(louvain.community)))
	modul = rbind(modul,c(i,"random",modularity(louvain.random)))
	#write modularity number
	write.table(t(c(i,modularity(louvain.community),modularity(louvain.random))),"Network_clusters/modularity-values.tsv",quote = F,sep = "\t",row.names = F,col.names = F,append = T)

	clusts = list()
	paths = list()
	for(name in lov.names){

		clust.genes = names(which(membership(louvain.community)==name))
		clusts[[as.character(name)]] = clust.genes
		write.table(clust.genes,paste("Network_clusters/clusts/nonr-",name,"clust-genes.tsv",sep = ""),quote = F ,sep = "\t",row.names = F,col.names = F)
		#enrich pathways our network		 
		entrez.genes = convert(clust.genes,id.map)
		x = enrichPathway(as.character(entrez.genes),pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(entrez))
		paths[[as.character(name)]] = x$Description
		
	}
	
	path.rand = list()
	for(name in rand.names){

		clust.genes = names(which(membership(louvain.random)==name))
		write.table(clust.genes,paste("Network_clusters/rand-clusts/",as.character(i),"-random-",name,"clust-genes-no subset.tsv",sep = ""),quote = F ,sep = "\t",row.names = F,col.names = F)
		#enrich pathways random network
		entrez.rand = convert(clust.genes,id.map)
		x = enrichPathway(as.character(entrez.rand),pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(entrez.rand))
		path.rand[[as.character(name)]] = x$Description
		
	}
	
	print(i)
	print(head(modul))

	#make set plot for our network
	path.list = unique(unlist(paths))
	set.matrix = make.set.mat(paths,path.list)
	colsums = apply(set.matrix,2,sum)

	if(sum(colsums) == 0 | length(which(colsums > 0)) < 2){
		
	}else{
	pdf(paste("Network_clusters/upset/clusts/whole-plot-",as.character(i),"-nonrandom-",name,".pfd",sep = ""),width = 15,height = 8)

	upset(set.matrix,sets.bar.color = "#56B4E9",order.by = "freq")

	dev.off()
	}
	#make set plot from random network

	rand.list = unique(unlist(path.rand))
	print(rand.list)
	set.rand = make.set.mat(path.rand,rand.list)
	table(set.rand)
	colsums = apply(set.rand,2,sum)

	if(sum(colsums) == 0 | length(which(colsums > 0)) < 2|length(rand.list) == 0){
		
	}else{

	pdf(paste("Network_clusters/upset/rand-clusts/set-plot-",as.character(i),"-random-",name,".pfd",sep = ""),width = 15,height = 8)

	upset(set.rand,sets.bar.color = "#56B4E9",order.by = "freq")

	dev.off()
	}
	
}

modul.tbl = data.frame(as.character(modul[,1]),as.character(modul[,2]),as.numeric(modul[,3]))

pdf("img2/modularity-boxplot.pdf")

colnames(modul.tbl) = c("run-no","status","modularity")
p <- ggplot(modul.tbl, aes(x=status, y=modularity,color=status)) + geom_boxplot() + labs(title= paste("modularity of random and non-random networks n =",as.character(nrow(modul.tbl)/2),sep = ""))

print(p)
dev.off()


wilcox.test(modul.tbl[which(modul.tbl$status == "non-random"),3],modul.tbl[which(modul.tbl$status == "random"),3])


