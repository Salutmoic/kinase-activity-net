#phosphositeplus and reactome data sets

reactome = read.delim("/home/borgthor/Kinase_activities/data/reactome.homo_sapiens.interactions.tab-delimited.txt")

reactome = reactome[intersect(grep("uniprot",reactome[,1]), grep("uniprot",reactome[,4])),]

reactome[,1] = gsub("uniprotkb:","",reactome[,1])
reactome[,4] = gsub("uniprotkb:","",reactome[,4])

#translate from uniprot symbols to gene symbols
library('org.Hs.eg.db')

a = select(org.Hs.eg.db, keys=reactome[,1] , columns= c('UNIPROT', 'SYMBOL'), keytype="UNIPROT")

b= select(org.Hs.eg.db, keys= reactome[,4], columns= c('UNIPROT', 'SYMBOL'), keytype="UNIPROT")

ksub = c()

for(i in 1:nrow(reactome)){

gene1 =reactome[i,1] 
gene2 =reactome[i,4] 

ksub = rbind(react,c(a[which(a[,1] == gene1)[1], 2],b[which(b[,1] == gene2)[1], 2]))

}


write.table("ksub","kinase_substrates.txt",quote =F,sep = "\t")


