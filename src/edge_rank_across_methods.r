#Thoughts on how to intergrate the correlation methods together, The arguenemts for all the functions is a list of the correlation matrices
#first method is to average the ranks of all the edges 
#Second is to scale the scores of each method to a score on the scale [0.1]
#Third methods takes the z scores of correlation coefficients within each matrix and 


averagerank <-function(matrices){
#matries is a list of correlation matrices which are ranked by corr coefficients or something equivalent (for instance, chisq score )
# returns matrix containing the average rank of each edge: kinase <-> kinase interaction

for(i in 1:length(matrices)){
rank = seq(1,nrow(matrices[[i]]))
matrices[[i]] = rbind(matrices[[i]],rank)

}
average.rank = c()


rankmat = cbind(as.character(matrices[[1]][,1]),as.character(matrices[[1]][,2]))

rows = unlist(lapply(matrices,nrow))
min = min(rows)
max = max(rows)

genepairs = paste(rankmat[,1],rankmat[,2])

ranks = c()

for(i in 1:length(genepairs)){

rank = 0
for(j in 1:4){

if(nrow(matrices[[j]]) == min){
rank = rank + which(paste(matrices[[j]][,1],matrices[[j]][,2]) == genepairs[i])
}

if(nrow(matrices[[j]]) == max){
#This is the funchisq matrix here I add the scores of A-B and B-A together, perhaps max() of those would be better?
rank = rank + which(paste(matrices[[j]][,1],matrices[[j]][,2]) == genepairs[i])
rank = rank + which(paste(matrices[[j]][,2],matrices[[j]][,1]) == genepairs[i])
}


}

rank = rank/5

ranks = c(ranks,rank)

}

rankmat = cbind(rankmat,ranks)
return(rankmat)

}


normalintergration <- function(matrices){
#matries is a list of correlation matrices which are ranked by corr coefficients or something equivalent (for instance, chisq score )
# returns the normalized sum of normalized values.

genepairs = paste(matrices[[1]][,1],matrices[[1]][,2])
normalvals = list()

rows = unlist(lapply(matrices,nrow))
min = min(rows)
max = max(rows)


for(i in matrices){

normalvals[[length(normalvals)+1]] = as.numeric(i[,3])/max(as.numeric(i[,3]))

}


ig.normvals =c()

for(i in 1:length(genepairs)){

nval = 0
for(j in 1:length(matrices)){

if(nrow(matrices[[j]]) == min){
nval = nval + normalvals[[j]][which(paste(matrices[[j]][,1],matrices[[j]][,2]) == genepairs[i])]
}

if(nrow(matrices[[j]]) == max){
#This is the funchisq matrix here I add the scores of A-B and B-A together
nval = nval + normalvals[[j]][which(paste(matrices[[j]][,1],matrices[[j]][,2]) == genepairs[i])]
nval = nval + normalvals[[j]][which(paste(matrices[[j]][,2],matrices[[j]][,1]) == genepairs[i])]
}

}
ig.normvals = c(ig.normvals,nval/5)
}




normat = cbind(as.character(matrices[[1]][,1]),as.character(matrices[[1]][,2]),ig.normvals)
return(normat)

}




zintergration <- function(matrices){

genepairs = paste(matrices[[1]][,1],matrices[[1]][,2])
zvals = list()

rows = unlist(lapply(matrices,nrow))
min = min(rows)
max = max(rows)


for(i in matrices){
#calc z vals
zvals[[length(zvals)+1]] = (as.numeric(i[,3])  - mean(as.numeric(i[,3])))/sd(as.numeric(i[,3]))

}

#combine z vals
ig.zvals =c()

for(i in 1:length(genepairs)){

zval = 0
for(j in 1:length(matrices)){

if(nrow(matrices[[j]]) == min){
zval = zval + zvals[[j]][which(paste(matrices[[j]][,1],matrices[[j]][,2]) == genepairs[i])]
}

if(nrow(matrices[[j]]) == max){
#This is the funchisq matrix here I add the scores of A-B and B-A together
zval = zval + zvals[[j]][which(paste(matrices[[j]][,1],matrices[[j]][,2]) == genepairs[i])]
zval = zval + zvals[[j]][which(paste(matrices[[j]][,2],matrices[[j]][,1]) == genepairs[i])]
}

}
ig.zvals = c(ig.zvals,zval/5)


}

zmat = cbind(as.character(matrices[[1]][,1]),as.character(matrices[[1]][,2]),pnorm(ig.zvals))
return(zmat)

}

