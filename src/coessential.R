
library(reshape2)
library(WGCNA)

bin <- function(d){

d.bin = rep(1,length(d))
d.bin[which(d > 1)] =2
d.bin[which(d < -1)] =0 

return(d.bin)
}


do.chisq <- function(kin1, kin2, tbl){
    if (kin1 == kin2){
        return(NA)
    }
    tbl = table(tbl[kin1,], tbl[kin2,]) 
    r<- chisq.test(tbl)
    
    return((-1)*log(r[[3]]))
}


argv <- commandArgs(TRUE)

data.file <- argv[1]
kinases <- arbv[2]
out.file <- argv[3]
method <- argv[4]

if(length(grep(".rds",data.file)) > 0 & paste(gsub("../data/external/","",data.file),"_temp-table.tsv",sep ="") %in% list.files()){

system(paste("Rscript format-drive-data.r",data.file,"../data/kinases-to-use.txt","temp-table.tsv",sep = " "))

tbl = read.delim(paste(gsub("../data/external/","",data.file),"_temp-table.tsv",sep =""))
}
else{
tbl <- read.table(data.file, as.is=TRUE)
kinases = read.delim(kinases)
}


tbl =tbl[tbl[,1] %in% kinases[,1],]

tbl = tbl[which(tbl[,1] %in% kinases[,1]),]

if(ncol(tbl) ==3){
names(tbl) <- c("kinase", "condition",  "SV")
tbl<- acast(tbl, kinase ~ condition)
}


if(method == "chisq"){

	bin.tbl <- apply(tbl,1,bin)
	v.do.chisq <- Vectorize(do.chisq, vectorize.args=c("kin1", "kin2"))
	res.tbl<- outer(rownames(bin.tbl), rownames(bin.tbl),v.do.chisq, bin.tbl)

	colnames(res.tbl) = rownames(tbl)
	rownames(res.tbl) = rownames(tbl)

	write.table(res.tbl,paste(out.file,"_",method,".txt"))
}


if(method = "WGCNA"){
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(tbl, powerVector = powers, verbose = 5)

softPower = sft[[2]][which(sft[[2]][,4] == max(sft[[2]][,4])),1]
adjacency = adjacency(tbl, power = softPower)
TOM = TOMsimilarity(adjacency,TOMType = "signed")

colnames(TOM) = rownames(tbl)
rownames(TOM) = rownames(tbl)

write.table(TOM,paste(out.file,"_",method,".txt"))
}


