
tissue_spec = function(data,t){

if( median(data) >= t){
print(1)
return(1)

}else{
print(0)
return(0)
}

}

argv <- commandArgs(TRUE)

e.th = argv[1]
a.th = argv[2]

expr = read.delim("TCGA/tcge_expr_brca_ov.txt")
samps = read.delim("TCGA/GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt",header = F)
rownames(expr) = expr[,1]
expr = data.matrix(expr[,2:ncol(expr)])
colnames(expr) = gsub("\\.","-",colnames(expr))

tfs = read.delim("TFs/TFs-human.txt")[,1]
kins = read.delim("TFs/TF-kinases.tsv")[,1]

e.tf = expr[rownames(expr) %in% tfs,]
e.kin = expr[rownames(expr) %in% kins,]

s.ov = samps[samps[,2] == "OV" ,]
s.b = samps[samps[,2] == "BRCA",]

tf.act = read.delim("TFs/tcga-TF-activities.tsv")

e.tf.o = e.tf[,colnames(e.tf) %in% s.ov[,1]]
e.tf.b = e.tf[,colnames(e.tf) %in% s.b[,1]]

ov = apply(e.tf.o,1,tissue_spec,t= e.th)
brca = apply(e.tf.b,1,tissue_spec,t= e.th)

tf.ov = rownames(e.tf.o)[which(ov == 1)]
tf.brca = rownames(e.tf.b)[which(brca == 1)]

save(tf.ov,file ="TFs/tfs-expressed-in-ov.RData")
save(tf.brca,file ="TFs/tfs-expressed-in-brca.RData")

e.kin.o = e.kin[,colnames(e.kin) %in% s.ov[,1]]
e.kin.b = e.kin[,colnames(e.kin) %in% s.b[,1]]

ov = apply(e.kin.o,1,tissue_spec,t= e.th)
brca = apply(e.kin.b,1,tissue_spec,t= e.th)

kin.ov = rownames(e.kin.o)[which(ov == 1)]
kin.brca = rownames(e.kin.b)[which(brca == 1)]

save(kin.ov,file ="TFs/kin-expressed-in-ov.RData")
save(kin.brca,file="TFs/kin-expressed-in-brca.RData")

act.ov = tf.act[,colnames(tf.act) %in% s.ov]
act.brca = tf.act[,colnames(tf.act) %in% s.b]

ov = apply(act.ov,1,tissue_spec,t= a.th)
brca = apply(act.brca,1,tissue_spec,t= a.th)

ao =  rownames(act.ov)[which(ov == 1)]
ab = rownames(act.brca)[which(brca == 1)]

save(tf.ov,file="TFs/tfs-active-in-ov.RData")
save(tf.brca,file="TFs/tfs-active-in-brca.RData")



