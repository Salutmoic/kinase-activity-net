suppressMessages(library("ReactomePA"))
suppressMessages(library("clusterProfiler"))
suppressMessages(library("UpSetR"))
suppressMessages(library("UpSetR"))

convert = function(gene.list,id.map){
	entrez.list = id.map[which(id.map[,1] %in% gene.list),2]
	return(entrez.list)
}

find.paths= function(tbl,path.tbl,file,universe,id.map){

	clusters = names(table(tbl[,2]))[which(table(tbl[,2])>= 10)]
	
	tbl.sect = c()
	for(cluster in clusters){
		
		genes = tbl[which(tbl[,2] ==cluster),1]
		entrez.genes = convert(genes,id.map)
		x = enrichPathway(as.character(entrez.genes),pvalueCutoff=0.05,pAdjustMethod = "BH", readable=T,universe = as.character(universe))
		paths = x$Description
		react = x$ID

		pvals = x$p.adjust
		
		if(length(paths) > 0){
			no.genes = sapply(x$GeneRatio,function(x) strsplit(x,"/")[[1]][1])
			geneids = x$geneID
 		}else{
			no.genes = NULL
			geneids = NULL
		}
		
		tbl.sect =rbind(tbl.sect,cbind(rep(cluster,length(paths)),paths,react,pvals,rep(length(clusters),length(paths)),no.genes,rep(length(genes),length(paths)),geneids))

	}
	
	if(!is.null(tbl.sect)){
		tbl.sect = cbind(rep(file,nrow(tbl.sect)),tbl.sect)
		path.tbl = rbind(path.tbl,tbl.sect)
	}
	return(path.tbl)
}

argv <- commandArgs(TRUE)
method = argv[1]

if(method == "louvain"){
	dirs = "clusters" 
	print(dirs)
}


if(method == "greedy"){
	dirs = list.dirs("clusters/greedy") 

	print(dirs)
}

if(method == "Separated"){
	dirs = "Separated_clusters"
}

if(method == "mcl"){
	dirs = list.dirs("clusters/mcl", recursive = F) 
	dirs = dirs[2:length(dirs)]

	print(dirs)
}



kinome = read.delim(argv[2])
kinome = unique(c(as.character(kinome[,1]),as.character(kinome[,2])))

id.map = read.delim(argv[3])
id.map = id.map[which(id.map[,3] == "Homo sapiens"),]

universe = convert(kinome,id.map)
path.tbl = c()

for(dir in dirs){
	print(dir)

	files = list.files(dir)
	files = files[grep("\\.tsv",files)]
	for(file in files){	
		path.tbl = c()

		print(file)
		tbl = read.delim(paste(dir,file,sep ="/"))
		path.tbl = find.paths(tbl,path.tbl,file,universe,id.map)
		print(dim(path.tbl))

	}

}

colnames(path.tbl) = c("file","cluster","pathways","reactome","p-values","no_clusters","no_genes_in_top_enriched","no_genes_in_cluster","geneids")
print(head(path.tbl))

write.table(path.tbl,paste("enrichments/reactome-",method,"-pathway-enrichments-round9.tsv",sep = ""),quote = F,sep = "\t",row.names = F)
write.table(path.tbl,paste("Random_test/reactome-",method,"-pathway-enrichments-round9.tsv",sep = ""),quote = F,sep = "\t",row.names = F)

