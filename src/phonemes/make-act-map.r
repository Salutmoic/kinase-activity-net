make.prior = function(netw,map){
	network = c()
	for(i in 1:nrow(netw)){
		sub = paste(as.character(netw[i,2]),"_s1",sep = "")
		kin = as.character(netw[i,1])
		nrow = c(gsub("_s1","",sub),gsub("_s1","",sub),kin,kin,as.character(map[which(map[,2] == sub),]$res),map[which(map[,2] == sub),]$pos,paste("e",nrow(network) +1,sep = ""), as.character(map[which(map[,2] == sub),]$S.cc))
		network = rbind(network,nrow)
	}	
	colnames(network) = c("S.AC","S.ID","K.AC","K.ID", "res", "pos", "SID", "S.cc")
	return(network)
}

argv <- commandArgs(TRUE)

val.set = read.delim(argv[1])
print(head(val.set))

if(length(argv) > 1){
	threshold = as.numeric(argv[2])
	val.set = val.set[which(as.numeric(val.set[,3]) > threshold),1:2]
}else{
	threshold = "omni"	
}

kins= read.delim("netw-kinases.tsv")
map = matrix(nrow=nrow(kins),ncol = 7)
colnames(map) = c("X", "dataID", "UPID", "site", "res",  "pos", "S.cc")

map[,1] = seq(1,nrow(kins))
map[,2] = paste(as.character(kins[,1]),"_s1",sep = "")
map[,3] = as.character(kins[,1])
map[,4] = rep("1s",nrow(kins))
map[,5] = rep("s",nrow(kins))
map[,6]= rep("1",nrow(kins))
map[,7] = paste(as.character(kins[,1]),"_s1",sep = "")
map = as.data.frame(map)
print(head(map))

network = make.prior(val.set,map)
rownames(network) = NULL
prior = network
print(head(prior))
write.table(map,"map-act.tsv",quote= F ,sep = "\t",row.names = F)
save("prior",file= paste("my-act-networks/prior-matrix-act-",as.character(threshold),".Rdata",sep = ""))








