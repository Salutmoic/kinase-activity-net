#phos stores kinase activiteis acorss conditions (Ochoa et al. 2016)
phos = read.delim("/nfs/research/petsalaki/users/borgthor/Kinase_activities/data/Supplementary_Table_2.csv")
ksub = read.delim("/nfs/research/petsalaki/users/borgthor/Kinase_activities/data/Kinase_Substrate_Dataset.gz",skip = 3)

ksub = ksub[which(ksub[,4] == "human" & ksub[,9] == "human"),]

kinases = phos[,1]
kinases = kinases[-1]


overlap = c()

v= 0

#see the overlap between each and every kinase pair
for(i in 1:214){
m =nrow(overlap)

for(j in (i+1):215){

k1 = ksub[which(toupper(ksub[,1]) == toupper(as.character(kinases[i]))),]
k2 = ksub[which(toupper(ksub[,1]) ==toupper(as.character(kinases[j]))),]

#paste modified residue to substrate name
subs1 = paste(k1[,5],k1$SUB_MOD_RSD)
subs2 = paste(k2[,5],k2$SUB_MOD_RSD)

subs1 =subs1[!duplicated(subs1)]
subs2 =subs2[!duplicated(subs2)]

inter = length(intersect(toupper(subs1),toupper(subs2)))

if(length(which(toupper(subs1) %in% toupper(subs2))) != inter){
break
}

row = c(as.character(kinases[i]),as.character(length(subs1)),as.character(kinases[j]),as.character(length(subs2)),as.character(inter), as.character(inter/(min(length(subs2),length(subs1)))))
overlap = rbind(overlap,row)

v = v+1
print(c(i,j))
}

if(is.null(m)){
m = 0
}

if(nrow(overlap) != m + (215-i)){
break
}

}

colnames(overlap) = c("Kinase","no substrates","Kinase","no. substrates","no. shared substrates","proportion of substrates shared")
overlap.sub = overlap[which(as.numeric(overlap[,6]) !=0), ]


overlap.sub = overlap.sub[which(as.numeric(overlap.sub[,2])  > 10 & as.numeric(overlap.sub[,4]) > 10), ]





overlap.sub = overlap.sub[order(as.numeric(overlap.sub[,6]),decreasing = TRUE),]
write.table(overlap.sub,"Overlap of kinase_substrate relationships.txt", sep = "\t",quote = F,row.names = FALSE)
