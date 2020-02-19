library(igraph)

argv = commandArgs(T)

if(length(grep("unknown",argv[1])) > 0){
	set.name = "unknown"
}else{
	set.name= "known"
}

paths = read.delim(argv[1],header = F)
starts = grep("x" ,paths[,1])
starts = c(starts,nrow(paths) +1)
val.set = read.delim("../data/validation-set-omnipath.tsv",header = F)
val.set2 = read.delim("../data/direct-validation-set-omnipath-1.tsv",header = F)
val.pairs = c(paste(val.set[,1],val.set[,2]))
net = read.delim("../out/merged-predictor.tsv")

vivo = read.delim("cutillas-pairs.tsv")
vitro = read.delim("haruna-full-int-revised.tsv")


path.pairs = c()
pathways = c()
psites = c()

i = 1
for(start in starts){

	finish = starts[i+1]
	path = c()
	path2 = c()	
	for(j in (start+1):(finish-2)){	
		path2 = c(path2,as.character(paths[j,1])) 
	}
	path = paths[(start+1):(finish-2),1]
	print(path)	
	print(table(path==path2))
		
	if(length(path)==1){
		i = i+1
		next
	}

#	pathways =c(pathways, paste(path,collapse=" "))

	t = 0
	psite = as.character(paths[(finish-1),1])
	print(psite)

	for(p in 1:(length(path)-1)){
		pair = paste(path[p],path[p+1])	
#		print(pair)
		if(!(pair %in% val.pairs)){
			print(pair)
			print(i)
			path.pairs = c(path.pairs,pair)
			psites = c(psites,psite)
			pathways =c(pathways, paste(path,collapse=" "))

			t=1
		}
	}
	if(t==1){
		print(path)
	}
	i = i+1
	if(i == length(starts)){
		break
	}
}

net = net[which(paste(net[,1],net[,2]) %in% path.pairs),1:3]


vivo.pairs = paste(vivo[,1],vivo[,2])
vitro.pairs = paste(vitro$kin_gene,vitro$gene)
int.pairs = intersect(vitro.pairs,vivo.pairs)
int.paths = pathways[which(path.pairs %in% int.pairs)]
int.ps = path.pairs[which(path.pairs %in% int.pairs)] 
int.psites = psites[which(path.pairs %in% int.pairs)]

cutoff = strsplit(argv[1],"-")[[1]][1]
wdt = strsplit(argv[1],"-")[[1]][2]

vivo.ps = path.pairs[which(path.pairs %in% vivo.pairs)]
vivo.psites = psites[which(path.pairs %in% vivo.pairs)]
vivo.paths = pathways[which(path.pairs %in% vivo.pairs)]

vitro.ps = path.pairs[which(path.pairs %in% vitro.pairs)]
vitro.psites = psites[which(path.pairs %in% vitro.pairs)]
vitro.paths = pathways[which(path.pairs %in% vitro.pairs)]

if(length(int.ps) > 0){
	write.table(data.frame(file = argv[1],cutoff=cutoff,path = int.paths,pairs = int.ps,ps = int.psites,evidence = "int",weight = wdt),paste("out/pik-mek-",set.name,"-new-regulations-PIK3-limma3.tsv",sep=""),quote=F,sep="\t",row.names = F,append = T)

}
if(length(vivo.ps) > 0){
	write.table(data.frame(file = argv[1],cutoff=cutoff,path = vivo.paths,pairs = vivo.ps,ps = vivo.psites,evidence = "vivo",weight = wdt),paste("out/pik-mek-",set.name,"-new-regulations-PIK3-limma3.tsv",sep=""),,quote=F,sep="\t",row.names = F,append = T)

}
if(length(vitro.ps) > 0){
	write.table(data.frame(file= argv[1],cutoff=cutoff,path = vitro.paths,pairs = vitro.ps,ps = vitro.psites,evidence = "vitro",weight = wdt),paste("out/pik-mek-",set.name,"-new-regulations-PIK3-limma3.tsv",sep=""),,quote=F,sep="\t",row.names = F,append = T)
}



