uniprot = read.delim("data/uniprot-id-map.tsv",header = F)

setwd("data/external/y2h/")

files = list.files()
interactome = c()

for(file in files){

mat = read.delim(file)
if(ncol(mat) == 42){
if(!is.null(interactome)){
print(table(colnames(mat) == colnames(interactome)))
}
interactome = rbind(interactome,mat)
}else{
print(file)
}
}

interactome = interactome[-which(interactome[,1] == "-" | interactome[,2] == "-"),]

mat = read.delim("LitBM-17.psi",header = F)

interactome.sub = interactome[,1:2]

names(mat) = names(interactome.sub)
interactome.sub = rbind(interactome.sub,mat[,1:2])

interactome.sub[,1] = gsub("uniprotkb:","",interactome.sub[,1])
interactome.sub[,1] = gsub("-.*","",interactome.sub[,1])

interactome.sub[,2] = gsub("uniprotkb:","",interactome.sub[,2])
interactome.sub[,2] = gsub("-.*","",interactome.sub[,2])

pairs = paste(interactome.sub[,1],interactome.sub[,2])

interactome.sub = interactome.sub[!duplicated(pairs),]

names(interactome.sub) = c("SymbolA", "SymbolB")

write.table(interactome.sub,"../../y2h-interactome.tsv",row.names = F,sep = "\t",quote = F)

