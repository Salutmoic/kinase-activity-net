
#Possible predictor, haven't tried yet. Based on the principle of synthetic lethality
#So if gene A is mutated in some cell lines and there is difference in survival between those cell lines and the rest when gene B is ko-> possible interaction between A and B
#note to self, incorporate effect size and.... the mutated sample is usually much smaller.


mut.data = read.delim("data/external/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf")

kinases = read.delim("data/human-kinome.txt")

mut.data = mut.data[mut.data[,1] %in% kinases[,1],]

#remove mutations that are unlikely to affect protein structure or expression level, this is open for debate
mut.data = mut.data[-which(mut.data$Variant_Classification %in% c("Silent","Intron")),]

achilles = read.delim("data/external/Achilles_v2.20.2_GeneSolutions.gct",skip = 2)

achilles = achilles[which(achilles[,1] %in% kinases[,1]),] 

rownames(achilles) = achilles[,1]
#I only want numeric values in achilles
achilles = achilles[,3:ncol(achilles)]

mut.matrix = matrix(0, nrow(achilles),ncol(achilles))
rownames(mut.matrix) = rownames(achilles)
colnames(mut.matrix) = colnames(achilles)

for(i in 1:nrow(mut.matrix)){

mut.sub = mut.data[(mut.data$Hugo_Symbol == rownames(mut.matrix)[i] ),]
mut.matrix[i,which(colnames(mut.matrix) %in% mut.sub$Tumor_Sample_Barcode)] = 1


}

m = apply(mut.matrix,1,sum)
mut.mat = mut.matrix[m > 20,]

p.scores = c()
kin.pairs = expand.grid(kinases[,1],kinases[,1]) 

for(i in 1:nrow(kinases)){

if(!(kinases[i,1] %in% rownames(mut.mat))){
p.scores = c(p.scores,rep("NA",nrow(kinases)))
next
}

m.row =as.numeric( mut.mat[which(rownames(mut.mat) == kinases[i,1] ),])

for(j in 1:nrow(kinases)){

if(!(kinases[j,1] %in% rownames(achilles))){
p.scores = c(p.scores,"NA")
next
}

a.row = as.numeric(achilles[which(rownames(achilles) == kinases[j,1]),])
cor = wilcox.test(a.row[which(m.row == 1)],a.row[which(m.row == 0)])

#take the log(p) for consistency
p.scores = c(p.scores, -log10(cor$p.val))

}
print(i)
}


cbind(kin.pairs,p.scores)
