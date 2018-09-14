library("preprocessCore")
library("ksea")

make.regulon = function(kinases,val.set,psites,row.no,data.rows){

    subs = paste(val.set$SUB_GENE, val.set$SUB_MOD_RSD)
    regulon = list()
    for(kinase in kinases[,1]){
        regs = data.rows[row.no[which(psites %in%  subs[as.character(val.set[,1]) ==kinase])]]
        if(length(regs) > 0){
            regulon[[kinase]] = unique(as.character(regs))
        }
    }

    return(regulon)
}

val.set = read.delim("Kinase_Substrate_Dataset_2018-05-01",skip =2)
val.set = val.set[val.set$KIN_ORGANISM == "human" & val.set$SUB_ORGANISM == "human",]

data = read.csv("Data_AND_Background_Network/log_measurments.csv")
kinases = read.delim("netw-kinases.tsv")

prots = strsplit(as.character(data[,1])," ")
m = lapply(prots,length)
prots.m = lapply(prots,function(x) paste(x[1],gsub("p-","",x[2])))

row.no = seq(1,nrow(data))
inds = which(m > 4)

for(i in inds){
	prot = prots[i]
	n = 2
	while(length(grep("\\(",prots[[i]][n])) == 0){
			psite = paste(prots[[i]][1],gsub("p-","",prots[[i]][n]))
			prots.m[[length(prots.m) + 1]] = psite
			row.no = c(row.no,i) 
			n = n+1
			
	}
}


psites = unlist(prots.m)
regulon = make.regulon(kinases,val.set,psites,row.no,data[,1])
m = lapply(regulon,length)

regulon = regulon[[unlist(m)>5]]

acts = c()

for(j in 2:ncol(data)){
names = data[-which(is.na(data[,j])),1]
vals = data[-which(is.na(data[,j])),j]
names = names[order(vals,decreasing=TRUE)]
vals = vals[order(vals,decreasing=TRUE)]
col = ksea_result = ksea_batchKinases(names,vals ,regulon , trial=1000)

acts = cbind(acts,unlist(col))

print(j)

}
#for(kinase in kinases[,1]){
#	row = c()
#	for(j in 2:ncol(data)){
#		if(length(which(!is.na(data[which(data[,1] %in% regulon[[kinase]]),j]))) ==0 ){
#			row = c(row,NA)
#			next
#		}
#		names = data[-which(is.na(data[,j])),1]
#		vals = data[-which(is.na(data[,j])),j]
#		names = names[order(vals,decreasing=TRUE)]
#		vals = vals[order(vals,decreasing=TRUE)]
#		signature = regulon[[kinase]][which(regulon[[kinase]] %in% names)]
#		if(length(signature) == 0){
#			row = c(row,NA)
#			next
#		}
#		ksea_result = ksea(names,vals ,signature , trial=1000, display=F,significance = TRUE)
#		row = c(row,as.numeric(ksea_result[1]))
			
#	}
#	print(kinase)
#	acts = rbind(acts,row)
#}

m  = rowSums(is.na(acts))
acts = acts[m<144,]


kinases = kinases[which(m<144),1]

acts = cbind(paste(as.character(kinases),"_s1",sep = ""),acts)
colnames(acts) = colnames(data)

write.table(acts,"act-tbl.tsv",row.names =F,sep = "\t",quote = F)

acts = read.delim("act-tbl.tsv") 
acts.norm = normalize.quantiles(as.matrix(acts[,2:ncol(acts)]))
colnames(acts.norm) = colnames(data)
acts.norm = cbind(paste(as.character(kinases),"_s1",sep = ""),acts.norm)
write.table(acts.norm,"act-tbl-norm.tsv",row.names =F,sep = "\t",quote = F)









