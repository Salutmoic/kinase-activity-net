library('mclust')

find_cutoff = function(pos,neg){

	cutoff = 0.99
	pos = pos[order(pos,decreasing = T)]

	for(i in 1:length(pos)){

		if(length(which(pos >= pos[i])) > length(which(neg >= pos[i]))){
		cutoff = (pos[i] + max(neg[which(neg <= pos[i])]))/2

		}

	}

return(cutoff)
}

define_clust = function(x){

	m = Mclust(x)
	params = m$parameters

	return( params)

}

calc.prob = function(x,pos,neg){

	pospar = define_clust(pos)
	negpar = define_clust(neg)
    
    pprob = 0
	for(i in 1:length(pospar$pro)){ 
       if(length(pospar$variance$sigmasq) == 1){
           variance =pospar$variance$sigmasq[1]
       }else{
           variance =pospar$variance$sigmasq[i]
       }
       pprob = pprob + dnorm(x,pospar$mean[i],sqrt(variance))*pospar$pro[i]
    }
 
    nprob = 0
	for(i in 1:length(negpar$pro)){
       if(length(negpar$variance$sigmasq) == 1){
           variance =negpar$variance$sigmasq[1]
       }else{
           variance =negpar$variance$sigmasq[i]
       }
       nprob = nprob + dnorm(x,negpar$mean[i],sqrt(variance))*negpar$pro[i]
    }

   return(pprob/(nprob+pprob))

}

argv <- commandArgs(TRUE)
pred.file = argv[1]
pred.m = read.delim(pred.file,header = F)
pred.pairs = paste(pred.m[,1],pred.m[,2])

pos.val = read.delim("data/validation-set-omnipath.tsv")

kinases = union(pred.m[,1],pred.m[,2]) 
kinpairs = expand.grid(kinases,kinases)
pos.pairs = paste(pos.val[,1],pos.val[,2])
neg.val = kinpairs[-which(paste(kinpairs[,1],kinpairs[,2]) %in% pos.pairs),]
neg.pairs = paste(neg.val[,1],neg.val[,2])

pred.filt = c()
clust.list = list()
pred.prob = c()

for(kinase in kinases){
	for(j in 3:ncol(pred.m)){

		kin.m = pred.m[pred.m[,1] == kinase,]
		kin.m = kin.m[!is.na(kin.m[,3]),]
		kin.pairs = paste(kin.m[,1],kin.m[,2])

		if(length(which(kin.pairs %in% pos.pairs))>3){

			print(kinase)
			pos.set = kin.m[kin.pairs %in% pos.pairs ,j]
			neg.set = kin.m[kin.pairs %in% neg.pairs,j]

			clust.prob = lapply(kin.m[,3],calc.prob,pos.set,neg.set )
			clust.list[[kinase]] = unlist(clust.prob)
            kin.m[,3] = unlist(clust.prob)
            pred.prob = rbind(pred.prob,kin.m) 
           
			cutoff = find_cutoff(pos.set,neg.set)
			print(cutoff)
			pf = pred.m[pred.m == kinase,]
			pf = pf[which(as.numeric(pf[,3]) > cutoff),]
			pred.filt = rbind(pred.filt,pf)

			print(pos.set)

			x=density(as.numeric(pos.set)) 
			y = density(as.numeric(neg.set))

			pdf(paste("img/comp-curves/",kinase,"-",colnames(pred.file)[j],"-prediction-curves.pdf"))

			hy =max(c(range(y[[2]])[2],range(x[[2]])[2])) 
			hx = max(c(range(y[[1]])[2],range(x[[1]])[2])) 
			lx = min(c(range(y[[1]])[1],range(x[[1]])[1])) 

			plot(x,main = paste(kinase,"Comparison between known and unknown edges n=",length(pos.set)),xlab = paste("score"),xlim = c(lx,hx),ylim = c(0,hy + 0.2*hy),col = "blue") 
			lines(y,col = "red") 
			points(pos.set,rep(0,length(pos.set)),col = "blue")
			points(neg.set,rep(0,length(neg.set)),col = "red")
			legend(hx*0.7,hy*0.7,c("known edges","random edges"), text.col = c("blue","red"), bty = "n") 
			dev.off()


		}	
	}
}

filt.pairs = paste(pred.filt[,1],pred.filt[,2])
print(table(filt.pairs %in% pos.pairs))
print(dim(pred.filt))

prob = pred.prob[pred.prob[,3] > 0.6,]
print(dim(prob))
ppairs = paste(pred.prob[,1],pred.prob[,2])
print(table( %in% pos.pairs))

write.table(pred.filt,"out/kinase-level-cutoff-network.tsv",quote =F,sep = "\t",row.names = F)
write.table(pred.prob,"out/kinase-level-prob-network.tsv",quote =F,sep = "\t",row.names = F)


