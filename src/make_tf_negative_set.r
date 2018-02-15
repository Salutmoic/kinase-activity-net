#Generation of negative validation set for Kinase-> TF interaction
#so far, I only remove known positives, might add other strategies
pos= read.delim("data/validation-TFs-omnipath.tsv")
kinases = read.delim("TFs/TF-kinases.tsv.tsv")
tfs = read.delim("TFs/TFs-human.txt")


pot.pairs = expand.grid(kinases[,1],tfs[,1])

kpairs = paste(pos[,1],pos[,2])

ppairs = paste(pot.pairs[,1],pot.pairs[,2])

neg.pairs = pot.pairs[-which(ppairs %in% kpairs),]

write.table(neg.pairs,"validation-TF-negative.tsv",sep = "\t",quote = F,row.names = F)
