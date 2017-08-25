signor = read.delim("data/external/human_phosphorylations_05_03_18.tsv")
spairs = cbind(as.character(signor$ENTITYA),as.character(signor$ENTITYB))
sp = paste(spairs[,1],spairs[,2])
spairs = spairs[-which(duplicated(sp)),]
sp = paste(spairs[,1],spairs[,2])

argv <- commandArgs(TRUE)
print(argv[1])

merge = read.delim(argv[1])
kinases = union(merge[,1],merge[,2])
kinpairs = expand.grid(kinases,kinases) 

kpairs = paste(kinpairs[,1],kinpairs[,2])

pos.set = kinpairs[kpairs %in% sp,]
neg.set = kinpairs[-which(kpairs %in% sp),]

write.table(pos.set,"validation-set-signor.tsv",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(neg.set,"validation-random-signor.tsv",quote = F,sep = "\t",row.names = F,col.names = F)

