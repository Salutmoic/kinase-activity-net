library("igraph")
library("ggplot2")
library("reshape2")
library("ReactomePA")
library("UpSetR")

network = read.delim("undir-netw/symmetric-netw.tsv")
network = network[,1:3]

modul = c()

network.sub = network[which(network[,3] > 0.5),]
network.sub[,3] = (network.sub[,3]-min(network.sub[,3]))/(max(network.sub[,3]) -min(network.sub[,3]))
nodes = unique(c(as.character(network.sub[,1]),as.character(network.sub[,2])))
graph  = graph_from_edgelist(as.matrix(network.sub[,1:2]), directed = F)
louvain.community = cluster_louvain(graph,weights = network.sub[,3] )
modul = rbind(modul,c("1","non-random",modularity(louvain.community)))

for(i in 1:1000){

	graph.random  = sample_degseq(as.numeric(degree(graph)),method = "vl")
	louvain.random = cluster_louvain(graph.random,weights = sample(network.sub[,3]))
	rand.names =  as.numeric(which(table(membership(louvain.random))>= 10))
	modul = rbind(modul,c(i,"random",modularity(louvain.random)))
	print("modularity")	
	print(modul)

}



write.table(modul,"modularity-stats.tsv",quote=F,sep = "\t",row.names = F)

modul = data.frame(modul[,2],as.numeric(as.character(modul[,3])))
vline = data.frame(value= modul[1,2], label ="Our Network")
colnames(modul) = c("Network_Type","Modularity")
modul = modul[2:nrow(modul),]


pdf("modularity_plot.pdf")
ggplot(modul,aes(modul[,2],size=1)) + stat_ecdf(geom = "step",color = "#56B4E9")+labs(x = "Network modularity",y ="F(modularity)") + theme_classic() +geom_vline(xintercept=vline$value, linetype="dashed", color = "red",show.legend=T)+geom_hline(yintercept=0.95, linetype="dashed", color = "black",show.legend=T)  +
  geom_text(aes(0,0.95,label = "y=0.95", vjust = -1)) + theme( axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20),axis.title.x=element_text(size = 20))
dev.off()
