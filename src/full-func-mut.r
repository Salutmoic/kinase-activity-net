library(FNN) 
library("FunChisq")
library('infotheo')
library(reshape2) 
library("OneR")

####
#funchisq Performance similar to spearman,will test with other predictors

###
suppressMessages(library(Biobase))
args = commandArgs(TRUE)


perm = function(x,y,n,method,num.clusts){

dist = c()
for(i in 1:n){
x.perm = sample(x)
y.perm = sample(y)

if(method == "Fnn"){
m= mutinfo(x.perm,y.perm,k= num.clusts)
dist = c(dist,m)
print("here")
}
if(method == "mutinfo"){
x.perm = bin(x.perm,num.clusts,labels = 1:num.clusts)
y.perm = bin(y.perm,num.clusts,labels = 1:num.clusts)

m= mutinformation(x.perm,y.perm)
dist = c(dist,m)
}
}

return(dist)
}

stat.test = function(x,y,method,num.clusts){
#cluster into bins for funchi

if(method == "funchisq"){
x=bin(x,num.clusts,labels = 1:num.clusts)
y =bin(y,num.clusts,labels = 1:num.clusts)
cont.tbl <- table(x, y)     
r <- fun.chisq.test(cont.tbl, method="nfchisq")     
return(r$p.value)

}

if(method == "Fnn"){

m= mutinfo(x,y,k = num.clusts)

dist = perm(x,y,100,method,num.clusts)
func  = ecdf(dist)
return(1-func(m))

}

if(method == "mutinfo"){
x=bin(x,num.clusts,labels = 1:num.clusts)
y =bin(y,num.clusts,labels = 1:num.clusts)

m= mutinformation(x,y)
dist = perm(x,y,100,method,num.clusts)
func  = ecdf(dist)

return(1-func(m))

}

}

kinact.file = args[1]
method = args[2]
num.clusts = as.numeric(args[3])

kinact.tbl = read.delim(kinact.file,sep = " ")
kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                as.is=TRUE)

kins1 = c()
kins2 = c()
scores = c()

kins = row.names(kinact.tbl)

for(i in 1:nrow(kinact.tbl)){

for(j in 1:nrow(kinact.tbl)){

if(i != j){

kin1cond = kinase.conditions[which(kinase.conditions$Kinase == kins[i]),"Condition"]
kin2cond = kinase.conditions[which(kinase.conditions$Kinase == kins[j]),"Condition"]

kin1val= kinact.tbl[i,] 
kin1val[kin1cond ] = "NA"
kin2val=kinact.tbl[j,]
kin2val[kin2cond] = "NA"

cond1 = which(!is.na(as.numeric(kin1val)))
cond2 = which(!is.na(as.numeric(kin2val)))

c.mut = intersect(cond1,cond2)

if(length(c.mut) >= max(5,num.clusts)){

p = stat.test(as.numeric(kin1val[c.mut]),as.numeric(kin2val[c.mut]),method,num.clusts)
score = -log(p+2.2e-16)
kins1 = c(kins1,kins[i])
kins2 = c(kins2,kins[j])
scores = c(scores,score)
}

}
print(j)
}
print(i)
}

res.tbl = cbind(kins1,kins2,scores)
colnames(res.tbl) = c("kin1","kin2","scores") 
write.table(res.tbl,paste("out/kinact-",method,"score.tsv",sep = ""),quote = F,sep ="\t",row.names =F)
