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

for(kinase in kinases){
for(j in 3:ncol(pred.m)){

kin.m = pred.m[pred.m[,1] == kinase,]
kin.pairs = paste(kin.m[,1],kin.m[,2])

if(length(which(kin.pairs %in% pos.pairs))>3){

print("hello")
print(kinase)
pos.set = kin.m[kin.pairs %in% pos.pairs ,j]
neg.set = kin.m[kin.pairs %in% neg.pairs,j]

pos.set = pos.set[!is.na(pos.set)]
neg.set = neg.set[!is.na(neg.set)]

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


print(dim(pred.filt))

write.table(pred.filt,"out/kinase-levle-filtered-network.tsv",quote =F,sep = "\t",row.names = F)





