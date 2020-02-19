library(openxlsx)
library(limma)
library( preprocessCore)


argv = commandArgs(T)
### data+annotation
data=read.xlsx("data/PI3Ki_MEKi.xlsx",sheet=1)
anno = read.xlsx("data/PI3Ki_MEKi.xlsx",sheet=2)

rownames(data) = data[,1]
data = data[,-1]
#remove rows with no proteins
data = data[-grep("09\\/|03\\/",rownames(data)),]
a =grep("\\(T|\\(S|\\(Y",rownames(data))

data = data[a,]
data = log(data,base=2)
data.norm = normalize.quantiles(data.matrix(data))
rownames(data.norm) = rownames(data)
colnames(data.norm) = colnames(data)
data=data.norm

table(anno[,1] == colnames(data))

control = grep("DMSO",anno[,2])
mek = grep("MEK",anno[,2])
pki= grep("PI3K",anno[,2])

data.mek = data[,c(mek,control)]
data.pki = data[,c(pki,control)]

########limma ####
#data.pki = log(data.pki,base =2)
#data.mek = log(data.mek,base =2)

if(argv[1] == "MEK"){
	design  = cbind(MEK=c(1,1,1,1,0,0,0,0),Cont=c(0,0,0,0,1,1,1,1))

#data.mek = log(data.mek,base =2)

	fit = lmFit(data.mek,design = design)

	cont.mat = makeContrasts(MEK-Cont,levels=design)
	fit2 = contrasts.fit(fit, cont.mat)
	fit2 = eBayes(fit2)

	table = as.data.frame(topTable(fit2, adjust.method = "BH",number = nrow(data.mek)))


}else{
	design  = cbind(PIK=c(1,1,1,1,0,0,0,0),Cont=c(0,0,0,0,1,1,1,1))

#	data.pki = log(data.pki,base =2)

	fit = lmFit(data.pki,design = design)

	cont.mat = makeContrasts(PIK-Cont,levels=design)
	fit2 = contrasts.fit(fit, cont.mat)
	fit2 = eBayes(fit2)

	table = as.data.frame(topTable(fit2, adjust.method = "BH",number = nrow(data.pki)))

}

table = table[order(rownames(table)),]

prots = c()
positions = c()
logfc = c()
pvals = c()

fdrs = c()
rots.logfc = c()

i = 1
for(row in rownames(table)){
	#Here I split rows with two ore more phosphosites
	proteins = strsplit(row,";")[[1]]
	for(protein in proteins){
		ps = strsplit(protein,"\\(")[[1]]
		p = ps[1]
		pos = gsub("\\)","",ps[2]) 
		if(length(grep("no_mod",pos)) > 0){
			next
		}
		prots = c(prots,p)
		positions= c(positions,pos)
		logfc = c(logfc,table$logFC[i])
		pvals = c(pvals,table$adj.P.Val[i])	
	}
	i = i+1
}


psite.table = data.frame(protein = prots,positions = positions,logfc = logfc, adj.pval = pvals)

write.table(psite.table,paste("out/",argv[1],"-all-psites-new.tsv",sep=""),quote=F,sep="\t",row.names=F)




