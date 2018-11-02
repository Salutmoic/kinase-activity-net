library(ggplot2)

merge = function(t,files,dir){
	files.sub = files[grep(paste(t,"-",sep = ""),files)]
	pairs = c()
	merged.netw = c()	
	
	for(file in files.sub){
		m = read.delim(paste(dir,file,sep = "/"))
		m[,1] = gsub("_s1","",m[,1])
		m[,3] = gsub("_s1","",m[,3])
		pairs.f = paste(m[,1],m[,3])
		merged.netw = rbind(merged.netw, m[which(!(pairs.f %in% pairs)),c(1,3)])
		pairs = c(pairs,pairs.f)
	}
	return(merged.netw)

}

make.comp.res = function(filt.net,net,val.set,t){
	
	val.pairs = paste(val.set[,1],val.set[,2])
	prun.pairs = paste(filt.net[,1],filt.net[,2])
	net.pairs = paste(net[,1],net[,2])
	
	net.edges= nrow(net)
	prun.edges = nrow(filt.net)
	net.prec = length(which(net.pairs %in% val.pairs))/length(net.pairs)
	prun.prec = length(which(prun.pairs %in% val.pairs))/length(prun.pairs)
	net.rec = length(which(net.pairs %in% val.pairs))/length(val.pairs)
	prun.rec = length(which(prun.pairs %in% val.pairs))/length(val.pairs)	
	
	return(c(t,net.rec,net.prec,net.edges,prun.rec,prun.prec,prun.edges))
}

argv <- commandArgs(TRUE)

dir = argv[1]
files = list.files(dir)

ts =seq(0.6,0.9,0.025)
network = read.delim("merged-predictor-annot.tsv")
val.set = read.delim("validation-set-omnipath.tsv")

results = c()
 
for(t in ts){
	net = merge(t,files,dir)
	netw = network[which(as.numeric(network[,3])>t),]
	row = make.comp.res(net,netw,val.set,t)	
	results = rbind(results,row)
}

colnames(results)= c("cutoff","net.recall","net.prec","net.edges","pruned.recall","pruned.prec","pruned.edges")

results = as.data.frame(results)
write.table(results,paste(dir,"results-comparison.tsv",sep =""),sep = "\t",quote=F,row.names =F)

pdf(paste(dir,"pruned-comparison-precision.pdf",sep=""))
boxplot(results$net.prec,results$pruned.prec,names = c("network","fitted-network"))

dev.off()

pdf(paste(dir,"-pruned-comparison-recall.pdf",sep =""))
boxplot(results$net.rec,results$pruned.rec,names = c("network","fitted-network"))

dev.off()


pdf(paste(dir,"pruned-comparison-f1-score.pdf",sep = ""))

f1.fit = 2*(results$pruned.prec*results$pruned.rec)/(results$pruned.prec+results$pruned.rec)
f1.rec = 2*(results$net.prec*results$net.rec)/(results$net.prec+results$net.rec)

boxplot(f1.rec,f1.fit,names = c("network","fitted-network"))

dev.off()


