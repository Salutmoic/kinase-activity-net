suppressMessages(library(igraph))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))


make_graph = function(netw,scaled){

	# Remove B-A if A-B is in the network
	net = t(apply(netw[,1:2],1,sort))
	inds = which(duplicated(net))
	if(length(inds)>0){
		netw = netw[-inds,]
	}
	#make graph from the edgelist
	graph  = graph_from_edgelist(as.matrix(netw[,1:2]), directed = F)
	weights = as.numeric(netw[,3])
	if(scaled == T){
		weights = scale(weights)
	}	
	return(list(graph,weights))
}



cluster = function(graph, method,w){
	
	if(method == "greedy"){
		community = cluster_fast_greedy(graph,weights = w)

	}
	if(method == "Louvain"){
		community = cluster_louvain(graph,weights = w)	
	}
	
	cnames =  names(table(membership(community)))
	clust.mat = c()

	for(i in 1:length(cnames)){
		print(cnames)
		genes = names(which(membership(community)==cnames[i]))
		gene.mat = cbind(genes,rep(i,length(genes)))
		clust.mat = rbind(clust.mat,gene.mat)
	}
	
	return(clust.mat) 
}


scale = function(values){
	print(head(values))
	scale = 1.0/(1.0-min(values))
	values = values-min(values)

	values = scale*values

	return(values)

}

argv <- commandArgs(TRUE)

netfile = argv[1]
network = read.delim(netfile)

network = network[,1:3]
print(head(network))
threshold = as.numeric(argv[2])
print(threshold)
method = argv[3]
print(method)

scaled = as.logical(argv[4])
separated = as.logical(argv[5])


network = network[which(as.numeric(as.character(network[,3])) >=threshold),]

if(separated == T){
	cl = read.delim(argv[6])
	network = network[which(network[,1] %in% cl[,1] & network[,2] %in% cl[,1]),]
}

graphv = make_graph(network,scaled)
w = graphv[[2]]
print(summary(w))
graph = graphv[[1]]

clusts = cluster(graph, method,w)

prefix = "" 
if(length(grep("pruned",netfile))>0){
	prefix = "pruned-"
}

if(length(grep("to-cluster",netfile))>0){
	prefix = paste(prefix,"to-cluster-",sep = "")
}

if(scaled == T){
	prefix = paste(prefix,"scaled-",sep = "")
}

if(separated ==F){
write.table(clusts,paste("clusters/",prefix,method,"-",threshold,".tsv",sep =""),sep = "\t",quote = F,row.names = F)

clusters = names(table(clusts[,2]))

for(cluster in clusters){
	write.table(clusts[which(clusts[,2] %in% cluster),],paste("Separated_clusters/",prefix,method,"-",threshold,"-",cluster,".tsv",sep =""),sep = "\t",quote = F,row.names = F)	
}

}else{
write.table(clusts,paste("Sperated_scluster_res/",prefix,method,"-",threshold,"-",length(list.files("Sperated_scluster_res"))+1,".tsv",sep =""),sep = "\t",quote = F,row.names = F)
}








