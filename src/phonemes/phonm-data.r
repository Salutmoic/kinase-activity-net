source("phon-dat-functions.r")

argv <- commandArgs(TRUE)

t.a = argv[1]
t.d = argv[2]
status = argv[3]

data = read.csv("Data_AND_Background_Network/log_measurments.csv")
uniprot = read.delim("uniprot-id-map.tsv",header = F)

prots = strsplit(as.character(data[,1])," ")
m = lapply(prots,length)
prots.m = lapply(prots,function(x) paste(x[1],x[2]))

inds = which(m > 4)
ps = list()
indexes = c()

if(status == "reg"){
#Only inlcude regulatory phosphosites in data
	regs = read.delim("Regulatory_sites_2018-05-01",skip = 2)
	regs = regs[which(regs$ORGANISM =="human"),]

	mod.rsd = strsplit(as.character(regs$MOD_RSD),"-")
        mrsd = lapply(mod.rsd,function(x) paste(x[2],x[1],sep = "-"))
	rown = paste(regs$PROTEIN,mrsd)

	for(i in inds){
	prot = prots[i]
	n = 2
	v = length(ps)
	
	while(length(grep("\\(",prots[[i]][n])) == 0){
			psite = paste(prots[[i]][1],prots[[i]][n])
			if(length(ps) ==v){	
				ps[[length(ps)+1]] = psite
			}else{
			#ps[[v+1]] = c(ps[[length(ps)]],psite)
				ps[[length(ps)+1]] = psite
			}
			n = n+1
			if(psite %in% rown){
				#indexes = c(indexes,i)
			}
		}
		indexes = c(indexes,i)
	}
	
	indexes = unique(indexes)
	rows = rep(0,nrow(data))
	rows[which(prots.m %in% rown)] =1
	rows[indexes] = 1

	data = data[which(rows==1),] 

}

if(status == "func"){
#Only consider top X% phosphosites in regards to function
	psites= read.delim("psites-to-use.tsv")
	psites = psites[psites[,2] > 0.2,]
	psites.genes = paste(lapply(psites[,1],function(x) as.character(uniprot[which(as.character(uniprot[,1]) == gsub("_.*","",x)),2])),"_HUMAN",sep = "")
	sites = lapply(psites[,1],function(x) strsplit(as.character(x),"_")[[1]][2])
	for(i in inds){
		prot = prots[i]
		n = 2
		v = length(ps)
	
		while(length(grep("\\(",prots[[i]][n])) == 0){
			psite = paste(prots[[i]][1],prots[[i]][n])
			if(length(ps) ==v){	
				ps[[length(ps)+1]] = psite
			}else{
				#ps[[v+1]] = c(ps[[length(ps)]],psite)
				ps[[length(ps)+1]] = psite
			}
			n = n+1
			#if(ps %in% rown){
				#indexes = c(indexes,i)
			#}
		}
		indexes = c(indexes,i)
	}

	prots.m = c(unlist(prots.m),unlist(ps))
}


netw = read.delim("separated-predictor-annot.tsv")
kins = unique(as.character(netw[,1]))

netw = netw[which(netw[,3] > t.a & netw$direct.pred.mean > t.d),]

d.rows = strsplit(as.character(data[,1])," ")
d.rows = lapply(d.rows,function(x) x[1])

data = data[which(d.rows %in% kins),]
d.rows = d.rows[which(d.rows %in% kins)]

my.map = make.map(data,d.rows)
my.map = as.data.frame(my.map)

if(status == "func"){
	my.map = my.map[which(paste(my.map$UPID,my.map$pos) %in% paste(psites.genes, sites)),]
}

data = data[data[,1] %in% my.map$dataID,]
d.rows = strsplit(as.character(data[,1])," ")
d.rows = lapply(d.rows,function(x) x[1])

netw = netw[which(netw[,2] %in% d.rows),]

load("Data_AND_Background_Network/allD_HUMAN_Omnipath.Rdata")
prior = make.prior(netw,allD,uniprot,data,d.rows,my.map)
prior = as.data.frame(prior)
allD = allD[which(allD$S.cc %in% my.map$S.cc),]



if(status == "func"){
	write.table(my.map,"my-map-func.tsv",quote= F ,sep = "\t",row.names = F)
	write.csv(data,"log-subset-func.csv",row.names = F)
	save("allD",file = "omni-prior-func.Rdata")
	save("prior",file= "prior-matrix-func.Rdata")
}else if(status == "reg"){
	write.table(my.map,"my-map-reg.tsv",quote= F ,sep = "\t",row.names = F)
	write.csv(data,"log-subset-reg.csv",row.names = F)
	save("allD",file = "omni-prior-reg.Rdata")
	save("prior",file= "prior-matrix-reg.Rdata")
}else{
	write.table(my.map,"my-map.tsv",quote= F ,sep = "\t",row.names = F)
	write.csv(data,"log-subset.csv",row.names = F)
	save("allD",file = "omni-prior.Rdata")
	save("prior",file= "prior-matrix.Rdata")
}

