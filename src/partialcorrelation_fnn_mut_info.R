library(ppcor)
library(FNN)
#known_interactions provides list of known kinase substrate information, protien-protein interactions
source("known_interactions")

files = list.files("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Kinase-activity/")

dir = "/nfs/research/petsalaki/users/borgthor/Kinase_activities/Kinase-activity/"

for(file in files){

#Read kinase activities
kin.table = read.table(paste(dir,file,sep = ""))

library(reshape2)
kin.num =acast(kin.table, kin.table[,1] ~ kin.table[,2])

#Find partial correlation between each kinase-pair
m = pcor(t(kin.num),method = "s")

scores = m$estimate

pvals = m$p.val


write.table(scores,paste("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Results/pcor/pcor coefficients",file,".txt"),quote = F, sep = "\t")

write.table(pvals,paste("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Results/pcor/pcor pvalues",file,".txt"),quote = F, sep = "\t")

#Here correlations of known kinase-substrate pairs are extracted
all = c()
known = c()

for(i in 1:(nrow(scores)-1)){

all = c(all,scores[i,(i+1):ncol(scores)])

#knownsubs stores the substrate of kinase i
knownsubs = ksub[which(ksub[,1] == rownames(scores)[i]),]
#knownkin has the kinases that react with kinase i
knownkin = ksub[which(ksub[,2] == rownames(scores)[i]),]

#combine lists of A->B and B-> A

subs = which(colnames(scores)[(i+1):ncol(scores)] %in% knownsubs[,2])+i
kins = which(colnames(scores)[(i+1):ncol(scores)] %in% knownkin[,1])+i
inter = c(subs,kins)
inter = inter[!duplicated(inter)]

known = c(known,scores[i,inter])


}





source("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Scripts/plot_distributions")

plot_dist(known,all,"pcor",u,file)


library(FNN)
#Find mutual information between each kinase-pair
mut.info = c()
gena = c()
genb = c()


for(i in 1:(nrow(kin.num)-1)){
for(j in(i+1):nrow(kin.num)){



mut.info = c(mut.info,mutinfo(as.numeric(kin.num[i,]),as.numeric(kin.num[j,])))

gena = c(gena,rownames(kin.num)[i])
genb = c(genb,rownames(kin.num)[j])

}
print(i)
}




scores = cbind(gena,genb,mut.info)

write.table(scores,paste("/nfs/research/petsalaki/users/borgthor/Kinase_activities/Results/fnn_mutinfo/FNNmutual information",file,".txt"),quote = F, sep = "\t")

source("check_if_known")





plot_dist(known,all,"FNN-mut_info",u,file)


}


