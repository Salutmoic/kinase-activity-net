library(ggplot2)
library(viridis)

enrichments = read.delim("enrichments/reactome-paths.tsv")
enrichments = enrichments[which(enrichments[,1] == "scaled-Louvain-0.5.tsv"),]

dists = list.files("reactome_2/")

dist.vals = c()
names = gsub(".txt","",dists)

i = 1

for(dist in dists){
	paths = strsplit(names[i],"_")[[1]] 
	path1 = paths[1]
	path2 = paths[2]
	
	if((path1 %in% enrichments$reactome) &  (path2 %in% enrichments$reactome)){
		m = read.delim(paste("reactome_2/",dist,sep = ""),header = F)
		index = which(enrichments$reactome == path1)
		index2 =which(enrichments$reactome == path2)
	
		value = m[1,1]
		dist.vals = rbind(dist.vals,c(path1,path2,value))

		i = i+1
		print(i)
	}
}

clusters = names(table(enrichments[,2])) 

within.cl = c()
rest = dist.vals

for(cluster in clusters){
	tbl =  enrichments[which(enrichments[,2] == cluster),]
	repaths = tbl$reactome
	inds = which(rest[,1] %in% repaths & rest[,2] %in% repaths)
	wc.sub = rest[inds,]
	if(length(inds)>0){rest = rest[-inds,]}
	print(dim(rest))
	indexes = c()
	if(nrow(wc.sub) > 1){
	for(i in 1:nrow(wc.sub)){
		en.row1 = enrichments[which(enrichments$reactome == wc.sub[i,1]),]
		en.row2 = enrichments[which(enrichments$reactome == wc.sub[i,2]),]
		
		genes1 = strsplit(as.character(en.row1$geneids),"/",fixed = T)[[1]]
		genes2 = strsplit(as.character(en.row2$geneids),"/",fixed = T)[[1]]
		
		if(length(intersect(genes1,genes2)) > 0){
			indexes = c(indexes,i)
		}else if(length(intersect(genes1,genes2)) == 0){
			print(cluster)
		}
		
	}
	wc.sub = wc.sub[-indexes,]
	}
	
	within.cl = rbind(within.cl,wc.sub)
}

print(wilcox.test(as.numeric(rest[,3]),as.numeric(within.cl[,3])))

measurment = c(rep("within clusters",nrow(within.cl)),rep("between clusters",nrow(rest)))

dist.mat = rbind(within.cl,rest)
print(dim(dist.mat))
dist.mat = cbind(dist.mat,measurment)
dist.mat = as.data.frame(dist.mat)
dist.mat[,3] = as.numeric(as.character(dist.mat[,3]))
colnames(dist.mat) = c("path1","path2","Distance","distance")


pdf("inter-intra-cluster-path-dist-2.pdf")
p= ggplot(dist.mat, aes(x=distance, y=Distance,fill=distance)) + geom_boxplot(outlier.size=0.25,show.legend = FALSE) + theme_classic() +scale_fill_manual("", values=c("#b3b3b3", "#8f2150")) +theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20),axis.title.y = element_text(size = 20),axis.title.x=element_blank())

print(p)

dev.off()





