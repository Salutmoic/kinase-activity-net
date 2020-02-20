
val.set = read.delim("../data/validation-set-omnipath.tsv",header = F)
val.set2 = read.delim("../data/direct-validation-set-omnipath-1.tsv",header = F)

val.pairs = paste(val.set[,1],val.set[,2])

val.pairs2 = paste(val.set2[,1],val.set2[,2])

m = read.delim("out3/pik-mek-unknown-new-regulations-PIK3-limma2.tsv")

m = m[m[,1] != "file",]

cyto.file = data.frame(prot1=character(),prot2 = character(),val.stat = character())

prot1 = c()
prot2 = c()
val.stat = c()

i = 1
for(path in m$path){
	prots = strsplit(path," ")[[1]]
	psite = as.character(m$ps[i])
	evi = as.character(m$evidence[i])
	new.pair = m$pairs[i]

	for(j in 1:(length(prots)-1)){
		prot1= c(prot1,prots[j])
		prot2 = c(prot2,prots[j+1])	
		pair = paste(prots[j],prots[j+1])		
		if(pair == as.character(new.pair)){
			val.stat = c(val.stat,evi)
		}else if(pair %in% val.pairs){
			val.stat = c(val.stat,"2_Src")
		}else if(pair %in% val.pairs2){
			 val.stat = c(val.stat,"1_Src")			
		}else{
			val.stat = c(val.stat,"unknown")			
		}
	}

	prot1 = c(prot1,prots[j+1])
	prot2 = c(prot2,psite)
	val.stat = c(val.stat,"pplus")
	i = i+1	
}

cyto.file = data.frame(prot1=prot1,prot2 = prot2,val.stat = val.stat)
cyto.pairs = paste(cyto.file[,1],cyto.file[,2],cyto.file[,3])
cyto.file = cyto.file[-which(duplicated(cyto.pairs)),]

write.table(cyto.file,"out3/new-net-known-cytofile.tsv",quote=F,sep="\t",row.names = F)


