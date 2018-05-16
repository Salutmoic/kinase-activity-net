library(igraph)
library(reshape2)
library("ReactomePA")
#same as pathway.r in general but uses different clustering methods

convert = function(gene.list,id.map){

entrez.list = id.map[which(id.map[,1] %in% gene.list),2]

}

id.map = read.delim("Entrez-ids.txt")
id.map = id.map[which(id.map$Species == "Homo sapiens"),]

argv <- commandArgs(TRUE)

netfile = argv[1]
network = read.delim(netfile)
directed = as.logical(argv[2])

network = network[,1:3]

adj.m <- acast(network, network[,1] ~ network[,2])
adj.m[which(is.na(adj.m))] =0

graph  = graph_from_adjacency_matrix(adj.m,mode = "undirected",weighted = T)

greedy.community = cluster_fast_greedy(graph,weights = E(graph)$weight )
louvain.community = cluster_louvain(graph,weights = E(graph)$weight )

lc2 = names(which(membership(louvain.community)==2))
lc1 = names(which(membership(louvain.community)==1))

gc1 = names(which(membership(greedy.community)==1))
gc2 = names(which(membership(greedy.community)==2))

coms = list(lc1,lc2,gc1,gc2)
names(coms) = c("lovain-1","lovain-2","greedy-1","greedy-2")
cluster.enrichment = list()

for(name in names(coms)){
genes = coms[[name]]
entrez = convert(genes,id.map)

x = enrichPathway(as.character(entrez),pvalueCutoff=0.05, readable=T)
method = strsplit(name,"-")[[1]][1]
pdf(paste("img2/",method,"/plot-",name,".pdf",sep =""),width = 15,height = 8)
g  =barplot(x,showCategory=10)
print(g)
dev.off()

cluster.enrichment[[name]] = x

}


#mcode
#markov
#imap


