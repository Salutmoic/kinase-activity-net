#Idea based on synthetic lethality
#Correlation of expression with shRNA ko survival
#correlation might be indicative of A-> B interaction

expr = read.delim("data/external/data_expression_median.txt")
kinases = read.delim("data/human-kinome.txt")
 
expr = expr[which(expr[,1] %in% kinases[,1]),]
rownames(expr) = expr[,1]
expr = expr[,3:ncol(expr)]

achilles = read.delim("data/external/Achilles_v2.20.2_GeneSolutions.gct",skip = 2)

achilles = achilles[which(achilles[,1] %in% kinases[,1]),]
rownames(achilles) = achilles[,1]
achilles = achilles[,3:ncol(achilles)]


common.cc = intersect(colnames(expr),colnames(achilles))

a.sub = achilles[,common.cc]
e.sub = expr[,common.cc]

e.sub = e.sub[,order(colnames(e.sub))]
a.sub = a.sub[,order(colnames(a.sub))]

cors = c()
kin.pairs = expand.grid(kinases[,1],kinases[,1])

for(i in 1:nrow(kin.pairs)){

if(!(kinases[i,1] %in% rownames(mut.mat))){
cors = c(cor,rep("NA",nrow(kinases)))
next
}

e.row =as.numeric( e.sub[which(rownames(e.sub)) == kinases[i,1] ),])

for(j in 1:nrow(kinases)){

if(!(kinases[j,1] %in% rownames(achilles))){
cors = c(cor,"NA")
next
}

a.row = as.numeric(a.sub[which(rownames(a.sub) == kinases[j,1]),])
cor = cor.test(a.row[which(m.row == 1)],a.row[which(m.row == 0)], method = "s")

#take the log(p) for consistency
cors = c(cors,cor$p.value)

}
print(i)

}

kin.table = cbind(kin.pairs,cors)
write.table(kin.table,"out/mutation_essential.tsv",sep = "\t",quote = F, row.names = F)



