#known_interactions provides list of known kinase substrate information, protien-protein interactions
source("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Scripts/known_interactions.R")

files = list.files("/nfs/research/petsalaki/users/borgthor//Kinase_activities/Kinase-activity/")

dir = "/nfs/research/petsalaki/users/borgthor//Kinase_activities/Kinase-activity/"
clust.methods = c("whole","by row","manual") 

for(file in files){

#Read kinase activities
kin.table = read.table(paste(dir,file,sep = ""))

library(reshape2)
kin.num =acast(kin.table, kin.table[,1] ~ kin.table[,2])

write.table(kin.num,paste("/nfs/research/petsalaki/users/borgthor/Kinase_activities/",substr(file,1,nchar(file)-4),"reshaped.txt",sep = ""),quote =F,sep ="\t")

# Values discretisized by various methods
for(u in clust.methods){
source("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Scripts/cluster.R")
if(u == "whole"){
discrete=whole(kin.num)
}


if(u == "by row"){
discrete=by_row(kin.num)
}

if(u == "manual"){
discrete= manual(kin.num)



}



library("FunChisq")

library("MASS")

#Functional Chi-sq values tested between each kinase-pair

funchisq = list()
pval = list()

gena = list()
genb = list()
ind = 1

for(i in 1:nrow(discrete)){
for(j in 1:nrow(discrete)){
if(i !=j){
tbl = table(as.numeric(discrete[i,]),as.numeric(discrete[j,])) 

r = fun.chisq.test(tbl,method = "nfchisq")

gena[[ind]] = rownames(kin.num)[i]
genb[[ind]]=rownames(kin.num)[j]
funchisq[[ind]] = r[[1]]
pval[[ind]] = r[[3]]

ind = ind +1
}
}
}


scores = cbind(gena,genb,funchisq,pval)

write.table(scores,paste("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Results/FunchiSq/funchisq_",file,"_", u,"_clustering matrix.txt",sep = ""),quote = F, sep = "\t")

#subsequent steps test if associations are known and plot the distribution of the socres of known-edges compared to random sample 
source("check_if_known")

vals =extract_known(scores,ksub)
known = vals[[1]]
all = vals[[2]]


source("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Scripts/plot_distributions")

plot_dist(known,all,"Funchisq",u,file)

library('infotheo')
#Find mutual information between each kinase-pair

mut.info = c()
gena = c()
genb = c()


for(i in 1:(nrow(discrete)-1)){
for(j in(i+1):nrow(discrete)){

if(i != j){

mut.info = c(mut.info,mutinformation(as.numeric(discrete[i,]),as.numeric(discrete[j,])))

gena = c(gena,rownames(kin.num)[i])
genb = c(genb,rownames(kin.num)[j])
}
}
print(i)
}


scores = cbind(gena,genb,mut.info)

write.table(scores,paste("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Results/Mut_info/mutual information",file, u,"clustering matrix.txt"),quote = F, sep = "\t")



#subsequent steps test if associations are known and plot the distribution of the socres of known-edges compared to random sample 


vals =extract_known(scores,ksub)
known = vals[[1]]
all = vals[[2]]


plot_dist(known,all,"mutual-information",u,file)

}
}
