library("ggplot2")
load("data/louvain-clusters.Rdata")

enriched.clusts = unique(unlist(clusts[c(2,3,4)]))
depleted.clusts = unique(unlist(clusts[c(1,5,6)]))

kin.citations <- read.table("data/kinome-entrez-data.tsv", as.is=TRUE) 
names(kin.citations) <- c("kinase", "citations", "citations.filt")

enrich.nos = cbind(kin.citations[which(as.character(kin.citations[,1]) %in% enriched.clusts),3],rep("Enriched"))
depleted.nos = cbind(kin.citations[which(as.character(kin.citations[,1]) %in% depleted.clusts),3],rep("Depleted"))
ctbl = rbind(enrich.nos,depleted.nos)

cit.no = data.frame(as.numeric(ctbl[,1]),as.character(ctbl[,2]))


pdf("img/cit-cluster-no-outliers.pdf")

colnames(cit.no) = c("citations","status")
p <- ggplot(cit.no, aes(x=status, y=citations,color=status)) + geom_boxplot() + labs(title= "number of citations per kinase in enriched vs depleted clusters")
ylim1 =  boxplot.stats(cit.no$citations)$stats[c(1, 5)]
p = p+  coord_cartesian(ylim = ylim1*1.5)

print(p)
dev.off()


