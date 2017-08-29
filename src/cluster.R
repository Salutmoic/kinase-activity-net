
library("mclust")

function <- whole(kin.num){
#kin.num stores kinase activities
#All values fed to Mclust as a one dimensional vector
p = as.vector(kin.num)
clust = Mclust( p )


discrete = matrix(clust$classification,nrow = nrow(kin.num),ncol = ncol(kin.num))
return(discrete)

}


function <- by_row(kin.num){

######
# kin.num stores the kinase activites acorss conditions
#discrete contains the discretized values
######

library(mclust,quietly=TRUE)

disc = c()

for(i in 1:nrow(kin.num)){

res = Mclust(kin.num[i,],G=3)$classification 


disc = rbind(disc,res)


}

discrete = disc
return(discrete)
}


function <- manual(kin.num){
#kin.num stores kinase activites acorss conditions

discrete = c()
for(i in 1:nrow(kin.num)){

row = rep(0,ncol(kin.num))
p = kin.num[i,]

row[which(p <= 1 & p >= -1)] = 2
row[which(p < -1)] = 1
row[which(p > 1)] = 3

if(0 %in% row){
break
}

discrete = rbind(discrete,row)
}

return(discrete)
}
