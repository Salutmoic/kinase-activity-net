
calc_average_degree = function(network){

	nodes = unique(c(network[,1],network[,2]))
	edges = c(network[,1],network[,2])
	degrees = lapply(nodes,function(x) length(which(edges == x)))

	return(median(unlist(degrees)))

}

ref_netw = read.delim("../BioPlex_interactionList_v4a.tsv")
ref_netw = ref_netw[,c("SymbolA","SymbolB")]

avg_netw = calc_average_degree(ref_netw)

netw = read.delim("merged-predictor-hinegs-annot.tsv")

ts = seq(0.7,0.5,-0.025)
averages = c()

for(t in ts){
net = netw[,1:3]
net = net[which(net[,3] > t),1:2]

averages = c(averages,calc_average_degree(net))

}



