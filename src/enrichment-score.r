library(ggplot2)
source("read-netw-files.r")

kol.test = function(network,name,subset = ""){
	if(subset != ""){
		if(!(subset %in% network[,1]) ){
			return(c(subset,NA,NA))
		}
		network = network[which(network[,1] == subset),]
	}
	if(nrow( network[which(network[,4] == TRUE ),]) < 3 ){
		return(c(subset,NA,NA))
	}
	results = ks.test(network[,3], network[which(network[,4] == TRUE ),3],alternative ="greater")

	return(c(name,results[[1]],results[[2]]))

}

kins = unique(network[,1])
network.undir = network[-which(is.na(network$assoc.pred.mean)),c("prot1","prot2","assoc.pred.mean","assoc.is.true")]
network.dir =network[-which(is.na(network$direct.pred.mean)),c("prot1","prot2","direct.pred.mean","direct.is.true")] 
network.sign = network[-which(is.na(network$prob.act.mean)),c("prot1","prot2","prob.act.mean","sign.is.correct")]
network.sign= network.sign[-which(is.na(network.sign[,4])),]
pdf("img/undir-netw-enrichment.pdf")

ggplot(network.undir, aes(x=prot1, y=assoc.pred.mean, fill= assoc.is.true)) + geom_boxplot()

dev.off()

pdf("img/dir-netw-enrichment.pdf")

ggplot(network.dir, aes(x=prot1, y=direct.pred.mean, fill= direct.is.true)) + geom_boxplot()

dev.off()


pdf("img/sign-netw-enrichment.pdf")

ggplot(network.sign, aes(x=prot1, y=prob.act.mean, fill= sign.is.correct)) + geom_boxplot()

dev.off()


colnames(network.dir) = colnames(network.undir)
colnames(network.sign) = colnames(network.undir)

network.whole = rbind(network.undir,network.dir,network.sign)

network.whole = cbind(network.whole,c(rep("undir",nrow(network.undir)),rep("dir",nrow(network.dir)),rep("signed",nrow(network.sign))))

colnames(network.whole) = c("prot1","prot2","prediction","Is_TRUE","Network_type")
pdf("img/whole-netw-enrichment.pdf")

ggplot(network.whole, aes(x=Network_type, y=prediction, fill= Is_TRUE)) + geom_boxplot()

dev.off()

res.undir =  kol.test(network.undir,"whole-net")
res.dir = kol.test(network.dir,"dir-net")
res.sign= kol.test(network.sign,"sign-net")
 
for(kin in kins){

	res.undir = rbind(res.undir, kol.test(network.undir,kin,subset = kin))
	res.dir = rbind(res.dir, kol.test(network.undir,kin,subset = kin))
	res.sign = rbind(res.sign, kol.test(network.undir,kin,subset = kin))

}

undir.adjust = p.adjust(as.numeric(res.undir[,3]),method ="bonferroni",n = length(res.undir[-which(is.na(res.undir[,3])),3])
)
res.undir = cbind(res.undir,undir.adjust)

dir.adjust = p.adjust(as.numeric(res.dir[,3]),method ="bonferroni",n = length(res.dir[-which(is.na(res.dir[,3])),3])
)
res.dir = cbind(res.dir,dir.adjust)

sign.adjust = p.adjust(as.numeric(res.sign[,3]),method ="bonferroni",n = length(res.sign[-which(is.na(res.sign[,3])),3])
)
res.sign = cbind(res.sign,sign.adjust)

colnames(res.undir) = c("kin","statistic","pval","padj")
colnames(res.dir) = c("kin","statistic","pval","padj")
colnames(res.sign) = c("kin","statistic","pval","padj")

write.table(res.undir,"out/undirected-predi-enrichment.tsv",quote = F,sep = "\t",row.names = F)

write.table(res.dir,"out/directed-predi-enrichment.tsv",quote = F,sep = "\t",row.names = F)

write.table(res.sign,"out/sign-predi-enrichment.tsv",quote = F,sep = "\t",row.names = F)


