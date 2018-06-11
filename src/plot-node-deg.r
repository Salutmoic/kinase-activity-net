library("ggplot2")

gen.random = function(kin.pairs,values){

	values = c(values,rep(0,nrow(kin.pairs)-length(values)))
	rand.net = data.frame(kin.pairs[,1],kin.pairs[,2],sample(values))

	rand.net = rand.net[rand.net[,3] > 0,]
	return(rand.net)
}

make.deg.mat = function(mat){
	kins = c(as.character(unique(mat[,1])))
	counts = c()
	for(kin in kins){
		counts = c(counts,length(which(mat[,1] == kin)))	
	}
	return(counts)
}

network = read.delim("full-annotated-predictor.tsv")

network[which(is.na(network[,3])),3] =0
network[which(network[,3] <= 0.8),3]= 0
network = network[,1:3]

res = c()

for(i in 1:100){

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
	
	counts.rand = make.deg.mat(net.rand)
	counts = make.deg.mat(net.sub)
	res = rbind(res,cbind(rep("random",length(counts.rand)),counts.rand))
	res = rbind(res,cbind(rep("not-rand",length(counts)),counts))	
}

df = data.frame(res[,1],as.numeric(res[,2]))
colnames(df) = c("status","counts")
head(df)
#is(df)

pdf("deg-density-plots.pdf")

g = ggplot(df,aes(x = counts, fill=status)) + geom_density()
print(g)

dev.off()

