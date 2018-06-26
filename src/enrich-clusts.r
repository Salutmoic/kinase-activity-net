suppressMessages(library("ReactomePA"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("UpSetR"))

convert = function(gene.list,id.map){
	entrez.list = id.map[which(id.map[,1] %in% gene.list),2]
	return(entrez.list)
}

find.paths= function(tbl,path.tbl,file,universe,id.map){

	clusters = which(table(tbl[,2])>= 10)
	
	tbl.sect = c()
	for(cluster in clusters){
		
		genes = tbl[which(tbl[,2] ==cluster),1]
		entrez.genes = convert(genes,id.map)
		x = enrichPathway(as.character(entrez.genes),pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(universe))
		paths = x$Description
		tbl.sect =rbind(tbl.sect,cbind(rep(cluster,length(paths)),paths,rep(length(clusters),length(paths))))
		
	}
	
	if(!is.null(tbl.sect)){
		tbl.sect = cbind(rep(file,nrow(tbl.sect)),tbl.sect)
		path.tbl = rbind(path.tbl,tbl.sect)
	}
	return(path.tbl)
}

argv <- commandArgs(TRUE)
methods = argv[1]

if(methods == "louvain-greedy"){
	dirs = list.dirs("clusters") 
	dirs = dirs[2:length(dirs)]

	print(dirs)
}

if(methods == "mcl"){
	dirs = c("mcl/benchmark-assoc-net-clust","mcl/benchmark-assoc-pruned-clust","mcl/benchmark-assoc-net-to-cluster-clust")
	
}


kinome = read.delim("human-kinome.txt")
id.map = read.delim("Entrez-ids.txt")
id.map = id.map[which(id.map[,3] == "Homo sapiens"),]
universe = convert(kinome[,1],id.map)
path.tbl = c()

for(dir in dirs){
	print(dir)

	files = list.files(dir)
	files = files[grep("\\.tsv",files)]
	for(file in files){
		print(file)
		tbl = read.delim(paste(dir,file,sep ="/"))
		path.tbl = find.paths(tbl,path.tbl,file,universe,id.map)
		print(dim(path.tbl))		
	}
}


write.table(path.tbl,paste(methods,"pathway-enrichments.tsv",sep = ""),quote = F,sep = "\t",row.names = F)
