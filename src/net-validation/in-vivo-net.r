library(openxlsx)
library(plyr)
library(ggplot2)
library(RColorBrewer )
library(forcats)


#Supp 
data.table = read.xlsx("data/41587_2019_391_MOESM5_ESM.xlsm")
#Predicted network
net = read.delim("../out/merged-predictor.tsv")
valset = rbind(read.delim("../data/validation-set-omnipath.tsv",header=F))

net = net[,1:3]
net = net[-which(is.na(net[,3])),]

valpairs = paste(valset[,1],valset[,2])

pairs = paste(data.table[,1],gsub("\\(.*","",data.table[,2]))
pairs = unique(pairs)
down.prot = lapply(data.table[,2],function(x) strsplit(x,"\\(")[[1]][1])
down.site= lapply(data.table[,2],function(x) strsplit(x,"\\(")[[1]][2])
down.site = gsub("\\)|S|T|Y","", down.site)

tbl = data.frame(kin_gene=as.character(data.table[,1]),gene=as.character(down.prot),site= as.character(down.site))

write.table(tbl,"cutillas-pairs.tsv",quote=F,sep="\t",row.names = F)

net.pairs = paste(net[,1],net[,2])
pairs.in.net = pairs[which(pairs %in% net.pairs)]
net.not.predicted = net[which(!(net.pairs %in% pairs)),]
net.not.predicted.pairs = paste(net.not.predicted[,1],net.not.predicted[,2])

##statistical testing of 
tbl.vect =rep(0,length(unique(pairs.in.net)))
tbl.vect[which(unique(pairs.in.net) %in% valpairs )]=1

net.vect= rep(0,nrow(net.not.predicted))
net.vect[which(net.not.predicted.pairs %in% valpairs)]=1

cont.tbl = matrix(c(c(table(tbl.vect)[2],table(tbl.vect)[1],c(table(net.vect)[2],table(net.vect)[1]))),nrow=2)

stat = fisher.test(cont.tbl,alternative="greater")
print(stat)


##############################
pairs = pairs[-which(pairs %in% valpairs)]

net.t=net[which(net.pairs %in% pairs),]
net.f= net[which(!(net.pairs %in% pairs)),]
net.f = net.f[which(!(net.pairs %in% valpairs)),]

pdf("cutillas=set.pdf")
boxplot(net.t[,3],net.f[,3],ylab = "probability",xlab = "set", names = c("cutillas-set","Unknown"))
dev.off()

phosfun = read.xlsx("data/media-3.xlsx?download=true")
uni.map = read.delim("data/uni-gene-map.tsv",header = F)
colnames(uni.map )=c("uni","gene","acc")
print(head(uni.map))
haruna.data = read.delim("haruna-full-kin-revised.tsv")

h.pairs = unique(paste(haruna.data$kin_gene,haruna.data$gene))
net.not.pred = net[-which(net.pairs %in% h.pairs),]
net.not.pred.pairs = paste(net.not.pred[,1],net.not.pred[,2])

cont.tbl = matrix(c(table(h.pairs %in% valpairs)[2],table(h.pairs %in% valpairs)[1],c(table(net.not.pred.pairs %in% valpairs)[2] ,table(net.not.pred.pairs %in% valpairs)[1])),nrow=2)

stat = fisher.test(cont.tbl,alternative="greater")

print(stat)
##############################

haruna.data = haruna.data[-which(paste(haruna.data$kin_gene, haruna.data$gene)  %in% valpairs),]

haruna.func = haruna.data[which(haruna.data$functional_score > 0.5) ,]

tbl = tbl[-which(paste(tbl$kin_gene, tbl$gene)  %in% valpairs),]

tbl = join(tbl,uni.map,by= "gene")
tbl = cbind(tbl,paste(tbl$uni,gsub("S|T|Y","",tbl$site)))
colnames(tbl)[6] = "psites"
phosfun = cbind(phosfun,paste(phosfun$uniprot, phosfun$position))
colnames(phosfun)[4] = "psites"
tbl = join(tbl,phosfun,by = "psites")
colnames(net)[1:2] = c("kin_gene","gene")
tbl = join(tbl,net,by= c("kin_gene","gene"))
tbl = tbl[-which(is.na(tbl$bart.pred.mean)),]

kin.cit = read.delim("../out/kinase-citations.tsv")
kin.cit = kin.cit[kin.cit[,2]=="num.subs",]
colnames(kin.cit)[1] = c("kin_gene")

tbl = join(tbl,kin.cit[,1:5],by="kin_gene")

tbl = tbl[-which(duplicated(paste(tbl[,1],tbl[,2],tbl[,3]))),]

tbl.func = tbl[which(tbl$functional_score >0.5),]

tbl.psites = cbind(tbl,paste(tbl[,2],tbl[,3]))
haruna.data.int = haruna.data[which(paste(haruna.data$kin_gene,haruna.data$gene) %in% paste(tbl$kin_gene,tbl$gene)),]
haruna.data.int.func = haruna.func[which(paste(haruna.func$kin_gene,haruna.func$gene) %in% paste(tbl.func$kin_gene,tbl.func$gene)),]

#basically escracting kinase -kinase pairs found in the different sets
h.pairs = paste(haruna.data$kin_gene,haruna.data$gene)
hf.pairs = paste(haruna.func$kin_gene,haruna.func$gene)
tbl.pairs =paste(tbl$kin_gene,tbl$gene)
tf.pairs =paste(tbl.func$kin_gene,tbl.func$gene)
hi.pairs = paste(haruna.data.int$kin_gene,haruna.data.int$gene)
hif.pairs = paste(haruna.data.int.func$kin_gene,haruna.data.int.func$gene)

net.f = net[-which(net.pairs %in% c(h.pairs,hf.pairs,tbl.pairs,tf.pairs,hi.pairs)),]
net.val = net[net.pairs %in% valpairs,]
net.h=net[which(net.pairs %in% h.pairs),]
net.hf= net[which(net.pairs %in% hf.pairs),]
net.t= net[which(net.pairs %in% tbl.pairs),]
net.tf= net[which(net.pairs %in% tf.pairs),]

#ggplot
ggplot.data = rbind(net.f,net.h,net.hf,net.t,net.tf,net.val)
sets = c(rep("bkg",nrow(net.f)),rep("vitro",nrow(net.h)),rep("vitro func",nrow(net.hf)),rep("vivo",nrow(net.t)),rep("vivo func",nrow(net.tf)),rep("val.set",nrow(net.val)))


ggplot.data = cbind(ggplot.data,sets)

names(ggplot.data)[3] = "probability"
print(head(ggplot.data))

pdf("gplot-net-val-valset-plot.pdf")
ggplot(ggplot.data,aes(x = reorder(sets,probability,FUN=mean),y = probability,fill = sets))+xlab("sets") + geom_boxplot() + theme_classic()+ scale_fill_brewer(palette="PRGn")
dev.off()



