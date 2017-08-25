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
