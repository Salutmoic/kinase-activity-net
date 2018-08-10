make.map = function(data,d.rows){
		
	map = matrix(ncol=7,nrow = nrow(data))
	colnames(map) =  c("X", "dataID", "UPID", "site","res", "pos", "S.cc")
	for(i in 1:nrow(data)){
		site = strsplit(as.character(data[i,1])," ")[[1]][2]
		site = gsub("p-","",site)
		upid =paste(d.rows[i],"HUMAN",sep = "_")
		row = c(i,as.character(data[i,1]),upid,site,substr(site,1,1),substr(site,2,nchar(site)),paste(upid,substr(site,1,1),substr(site,2,nchar(site)),sep = "."))	
		map[i,] = row
	}	
	return( map)	
}


make.prior = function(netw,allD,uniprot,data,d.rows,map){
	
	uni.d = allD[,c(1,2)]
	uni.p = allD[,c(3,4)]
	colnames(uni.p) = colnames(uni.d)
	uni = rbind(uni.d,uni.p)
	prior = c()
	m = 0
	for(i in 1:nrow(netw)){
		s.id = paste(as.character(netw[i,2]),"HUMAN",sep ="_")
		if(s.id %in% as.character(uni[,2])){
			s.ac = as.character(uni[which(as.character(uni[,2]) ==as.character(s.id))[1],1])
		}else{
			s.ac = as.character(uniprot[which(as.character(uniprot[,2]) == as.character(netw[i,2]))[1],1])
		}
		k.id = paste(as.character(netw[i,1]),"HUMAN",sep = "_")
		if(k.id %in% as.character(uni[,2])){
			k.ac = as.character(uni[which(as.character(uni[,2]) ==as.character(k.id))[1],1])
		}else{
			k.ac = as.character(uniprot[which(as.character(uniprot[,2]) == as.character(netw[i,1]))[1],1])
		}
		map.rows = data.frame(map[which(map$UPID == s.id),])
		rows = c()
		n = 1
		for(h in 1:nrow(map.rows)){
			r = c(s.ac,s.id,k.ac,k.id,as.character(map.rows[h,]$res),as.character(map.rows[h,]$pos),
			paste("e",as.character(m+n),sep= ""),as.character(map.rows[h,]$S.cc))			
			rows = rbind(rows,t(r))
			n = n+1		
		}		
		prior = rbind(prior,rows)		
		m= nrow(prior)
	}
	colnames(prior) = c("S.AC", "S.ID", "K.AC", "K.ID","res", "pos", "SID", "S.cc")
	return(prior)
}

argv <- commandArgs(TRUE)

data = read.csv("Data_AND_Background_Network/log_measurments.csv")
#regs = read.delim("Regulatory_sites_2018-05-01",skip = 2)
#regs = regs[which(regs$ORGANISM =="human"),]

#mod.rsd = strsplit(as.character(regs$MOD_RSD),"-")
#mrsd = lapply(mod.rsd,function(x) paste(x[2],x[1],sep = "-"))
#rown = paste(regs$PROTEIN,mrsd)

prots = strsplit(as.character(data[,1])," ")
m = lapply(prots,length)
prots.m = lapply(prots,function(x) paste(x[1],x[2]))

inds = which(m > 4)
ps = list()
indexes = c()

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
psites= read.delim("psites-to-use.tsv")
uniprot = read.delim("uniprot-id-map.tsv",header = F)

psites.genes = paste(lapply(psites[,1],function(x) as.character(uniprot[which(as.character(uniprot[,1]) == gsub("_.*","",x)),2])),"_HUMAN",sep = "")
sites = lapply(psites[,1],function(x) strsplit(as.character(x),"_")[[1]][2])

#indexes = unique(indexes)
#rows = rep(0,nrow(data))
#rows[which(prots.m %in% rown)] =1
#rows[indexes] = 1

#data = data[which(rows==1),] 

netw = read.delim("separated-predictor-annot.tsv")
t.a = argv[1]
t.d = argv[2]
kins = unique(as.character(netw[,1]))

netw = netw[which(netw[,3] > t.a & netw$direct.pred.mean > t.d),]

d.rows = strsplit(as.character(data[,1])," ")
d.rows = lapply(d.rows,function(x) x[1])

data = data[which(d.rows %in% kins),]
d.rows = d.rows[which(d.rows %in% kins)]

my.map = make.map(data,d.rows)
my.map = as.data.frame(my.map)
my.map = my.map[which(paste(my.map$UPID,my.map$pos) %in% paste(psites.genes, sites)),]

data = data[data[,1] %in% my.map$dataID,]
d.rows = strsplit(as.character(data[,1])," ")
d.rows = lapply(d.rows,function(x) x[1])

netw = netw[which(netw[,2] %in% d.rows),]

write.table(my.map,"my-map-func.tsv",quote= F ,sep = "\t",row.names = F)
load("Data_AND_Background_Network/allD_HUMAN_Omnipath.Rdata")
prior = make.prior(netw,allD,uniprot,data,d.rows,my.map)
prior = as.data.frame(prior)
save("prior",file= "prior-matrix-func.Rdata")
allD = allD[which(allD$S.cc %in% my.map$S.cc),]
save("allD",file = "omni-prior-func.Rdata")

write.csv(data,"log-subset-func.csv",row.names = F)
