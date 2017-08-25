library(org.Hs.eg.db)
library(plyr)
library(openxlsx)
library(data.table)
library(dplyr)

#Our kin-kin network+ included kinases
kinome = read.delim("../data/human-kinome.txt",header =F)
net= read.delim("../out/merged-predictor.tsv")
net = net[-which(is.na(net[,3])),]

#Downloaded from ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
uniprot= read.delim("data/HUMAN_9606_idmapping.dat",header=F)

unis= uniprot[grep("UniProtKB-ID",uniprot[,2]),]
genes=uniprot[grep("Gene_Name",uniprot[,2]),]
#functional score uniprot+site+ value
phosfun = read.xlsx("media-3.xlsx?download=true")

#these include kinase id mapping from the id they used from https://www.carnabio.com/english/product/protein4.html, ids mapped with unprot, and manually
uni.kins= read.delim("out/haruna-kins-uniprot.tsv")
manual.kins = read.delim("alias-gene-2.tsv")

#kinase citations to look at
kin.cit = read.delim("../out/kinase-citations.tsv")

#loads of kinases are part of a complex, I split the kinases in these rows into two/three
if(!("res-split-tbl-2.tsv" %in% list.files())){
	res = read.delim("data/kin-sub-relationship")
	res = res[-1,]
	res = res[-grep("mutant",res[,1]),] 
	rows = grep("/",res[,2])

	for(row.no in rows){
		row = res[row.no,]
		kins= strsplit(as.character(unlist(row[2])),"/")[[1]]
		if(kins[1] == "AMPKa2"){kins=c("AMPKa2","AMPKb1","AMPKg1")}
		if(kins[1] == "AMPKa1"){kins=c("AMPKa1","AMPKb1","AMPKg1")}
		for(i in 1:length(kins)){
			row[2] = kins[i]
			print(row)
			print(nrow(res))
			res = rbind(res,row)
			print(i)
		}
	}

	#remove rows with the ambigous status
	res = res[-rows,]
	write.table(res,"res-split-tbl-2.tsv",quote=F,sep="\t",row.names = F)
}else{
	res = read.delim("res-split-tbl-2.tsv")
}

###################
#convert uniprot to managebla id-map
#and add tot able
######
unis = unis[which(unis[,1] %in% genes[,1]),]
genes = genes[which(genes[,1] %in% unis[,1]),]

table(unis[,1] == genes[,1])

tbl = join(genes[,c(1,3)],unis[,c(1,3)],by = "V1")

colnames(tbl) = c("uniprot","gene","Uniprot.ID")
write.table(tbl,"data/uni-gene-map.tsv",quote=F,sep="\t",row.names = F)

res = join(res,tbl,by = "Uniprot.ID")
#subsequently I add functionality of phosphosites
psites.res = paste(res$uni,gsub("T|Y|S","",res$Position))
psites.func = paste(phosfun$uniprot,phosfun$position)

res = cbind(res,psites.res)
colnames(res)[ncol(res)] = "psites"
phosfun = cbind(phosfun,psites.func)
colnames(phosfun)[ncol(phosfun)]= "psites"

res = join(res,phosfun,by="psites")

#here I map kinases onto gene -names and add to the res table
manual.kins = cbind(alias=as.character(manual.kins[,1]),gene = as.character(manual.kins[,2]))

for(i in 1:nrow(manual.kins)){
	row = as.character(unlist(manual.kins[i,]))
	if(row[1] %in% kinome[,1] & !(row[2] %in% kinome[,1])){
		print("here")	
		print(row)
		print(i)
		manual.kins[i,2] = row[1]
	}
}

manual.kins = data.frame(alias=as.character(manual.kins[,1]),gene = as.character(manual.kins[,2]))

uni.kins = uni.kins[-1,]
uni.kins = uni.kins[-grep(",",uni.kins[,1])]

Gene = c()
for(i in 1:nrow(uni.kins)){
	gene.names = as.character(uni.kins$Gene.names[i])
	gene.names  = strsplit(gene.names," ")[[1]]
	Gene = c(Gene,NA)
	if(as.character(uni.kins[i,1]) %in% kinome[,1]){
		Gene[i] = as.character(uni.kins[i,1])
	}else{	
		t=0	
		for(genes in gene.names){
			if( genes %in% kinome[,1] ){
				t = t+1
				Gene[i]= genes
			}
		}	
		print(t)
		if(t==2){
			print(uni.kins[i,])
			print(gene.names)
		}
	}
}

uni.kins = cbind(uni.kins,Gene)
uni.kins = uni.kins[,c("yourlist.M20190917E5A08BB0B2D1C45B0C7BC3B55FD26556446D65U","Gene")]
colnames(manual.kins) = c("Kinase","kin_gene")
colnames(uni.kins)= colnames(manual.kins)
kin.gene.tbl = rbind(uni.kins,manual.kins)
table(res[,2] %in% kin.gene.tbl[,1])

res = join(res,kin.gene.tbl,by = "Kinase")
head(res)

#### this is to only include pho

#Here I add the bart predictions and citations behind upstream kinase to the reults
colnames(net)[1:2] = c("kin_gene" ,"gene")
res = join(res,net[,1:3],by = c("kin_gene" ,"gene"))
write.table(res,"out/haruna-full-int-tbl.tsv",quote=F,sep="\t",row.names = F)
res = res[-which(is.na(res$bart.pred.mean)),]

kin.cit = kin.cit[kin.cit[,2]=="num.subs",]
colnames(kin.cit)[1] = "kin_gene"
res = join(res,kin.cit[,c(1,5)],by="kin_gene")
write.table(res,"haruna-kins-res-3.tsv",quote=F,sep="\t", row.names = F)
res.zero = res[which(res$count  ==0),]

res.zero= res.zero[,c("kin_gene","gene","Position","position","functional_score","bart.pred.mean", "count")]
write.table(res.zero,"zero-kin-tbl-3.tsv",quote=F,sep="\t",row.names = F)

res.zero = data.table(res.zero)

res  = res.zero %>% group_by(kin_gene,gene) %>% slice(which.max(functional_score))
res[order(res$bart.pred.mean,decreasing = T),]

write.table(res,"zero-sub-kin-top-prediction.tsv",quote=F,sep="\t",row.names = F)
#####

