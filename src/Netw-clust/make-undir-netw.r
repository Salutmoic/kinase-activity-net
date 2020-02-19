argv <- commandArgs(TRUE)

network = read.delim(argv[1])
if(length(which(is.na(network[,3])))){
	network = network[-which(is.na(network[,3])),]
}
prot.pairs = t(apply(network[,1:2],1,sort))
network[,1:2] = prot.pairs
network = network[order(network[,1],network[,2],as.numeric(as.character(network[,3])),decreasing = T),]

pairs = paste(network[,1],network[,2])
sym.netw = network[which(!duplicated(pairs)),1:3]

colnames(sym.netw) = c("prot1","prot2","score")
write.table(sym.netw,"undir-netw/symmetric-netw.tsv",row.names = F,sep = "\t",quote =F)
