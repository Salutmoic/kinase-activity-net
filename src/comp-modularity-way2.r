library("igraph")
library("ggplot2")
library("reshape2")
library("ReactomePA")
library("UpSetR")

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

network = read.delim("undir-netw/symmetric-netw.tsv")
#network = network[-which(is.na(network[,3])),3]
network = network[,1:3]

modul = c()

network.sub = network[network[,3] > 0.5,]
network.sub[,3] = (network.sub[,3]-min(network.sub[,3]))/(max(network.sub[,3]) -min(network.sub[,3]))
nodes = unique(c(as.character(network.sub[,1]),as.character(network.sub[,2])))
graph  = graph_from_edgelist(as.matrix(network.sub[,1:2]), directed = F)
louvain.community = cluster_louvain(graph,weights = network.sub[,3] )
E(graph)$weight = network.sub[,3]
modul = rbind(modul,c("1","non-random",modularity(louvain.community)))

for(i in 1:1000){

	print(i)
	#net.rand = network.sub
	#net.rand[,3] = sample(net.rand[,3])
	
	#graph.random  = graph_from_edgelist(as.matrix(net.rand[,1:2]), directed = F)

	#cluster
	#graph.random = permute(graph, sample(length(nodes))) 
	graph.random  = sample_degseq(as.numeric(degree(graph)),method = "vl")
	louvain.random = cluster_louvain(graph.random,weights = E(graph)$weight)
	rand.names =  as.numeric(which(table(membership(louvain.random))>= 10))
	modul = rbind(modul,c(i,"random",modularity(louvain.random)))
	print("modularity")	
	print(modul)
#	path.rand = list()
#	for(name in rand.names){

#		clust.genes = names(which(membership(louvain.random)==name))
		#write.table(clust.genes,paste("Network_clusters/rand-clusts/",as.character(i),"-random-",name,"clust-genes-no-subset.tsv",sep = ""),quote = F ,sep = "\t",row.names = F,col.names = F)
		#enrich pathways random network
#		entrez.rand = convert(clust.genes,id.map)
#		print(length(entrez.rand))
#		x = enrichPathway(as.character(entrez.rand),pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(entrez.rand))
#		path.rand[[as.character(name)]] = x$Description
		
#	}

#	rand.list = unique(unlist(path.rand))
#	print(rand.list)
#	set.rand = make.set.mat(path.rand,rand.list)
#	table(set.rand)
#	colsums = apply(set.rand,2,sum)

#	if(sum(colsums) == 0 | length(which(colsums > 0)) < 2|length(rand.list) == 0){
		
#	}else{

#	pdf(paste("Network_clusters/upset/rand-clusts/set-plot-",as.character(i),"-random-",name,".pfd",sep = ""),width = 15,height = 8)

#	upset(set.rand,sets.bar.color = "#56B4E9",order.by = "freq")

#	dev.off()
#	}
	

}

save(path.rand,filename= "random-enrichments.Rdata")

write.table(modul,"modularity-stats.tsv",quote=F,sep = "\t",row.names = F)


modul = data.frame(modul[,2],as.numeric(as.character(modul[,3])))
vline = data.frame(value= modul[1,2], label ="Our Network")
#modul[,1] = c("Our Network",rep("Randomized Networks",nrow(modul)-1))
colnames(modul) = c("Network_Type","Modularity")
modul = modul[2:nrow(modul),]

#ggplot(modul,aes(modul[,2],color = Network_Type)) + stat_ecdf(geom = "point")+stat_ecdf(geom = "step") + labs(x = "modularity",y ="F(modularity)") + theme_classic()


pdf("modularity_plot.pdf")
ggplot(modul,aes(modul[,2])) + stat_ecdf(geom = "step",color = "#56B4E9")+labs(x = "Network modularity",y ="F(modularity)") + theme_classic() +geom_vline(xintercept=vline$value, linetype="dashed", color = "red",show.legend=T)+geom_hline(yintercept=0.95, linetype="dashed", color = "black",show.legend=T)  +
  geom_text(aes(0,0.95,label = "y=0.95", vjust = -1))
dev.off()
