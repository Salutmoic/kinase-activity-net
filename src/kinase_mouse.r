
kinome = read.delim("../../data/external/Kinase_Substrate_Dataset_2017-09-08",skip = 2)

kin.mouse = kinome[kinome$KIN_ORGANISM == "mouse",]

write.table(kin.mouse,"../data/mouse-kinome.tsv",quote = F,sep = "\t")

kinases = unique(kin.mouse[,1])
write.table(kinases,"../data/mouse-kinases.tsv",quote = F,sep = "\t")


system("Rscript ../../src/kinase-overlap.r ../data/mouse-kinome.tsv ../data/mouse-kinases.tsv")

ktable = read.delim("../../data/kinase_substrate_mouse_overlap.tsv")

red.kin.t = ktable[ktable[,6] > 0.5,]

red.kins = c()

for(i in 1:nrow(red.kin.t)){

if(red.kin.t[i,1] %in% red.kins | red.kin.t[i,3] %in% red.kins){
next
}

if(red.kin.t[i,2] < red.kin.t[i,4]){
red.kins = c(red.kins,as.character(red.kin.t[i,1]))

}else{

red.kins = c(red.kins,as.character(red.kin.t[i,3]))

}

}

kinases = kinases[-which(kinases %in% red.kins)]

write.table(kinases,"../data/mouse-kinases.tsv",sep = "\t",quote = F)

