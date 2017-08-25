sml.results <- read.delim("out/go-sem-sim-uniprot.tsv", as.is=TRUE)
uniprot.map <- read.table("data/uniprot-id-map.tsv", as.is=TRUE)
colnames(uniprot.map) <- c("uniprot", "hgnc")

names(sml.results) <- c("prot1", "prot2", "sem.sim")
sml.results$sem.sim[which(sml.results$sem.sim==-1.0)] <- NA

sml.results$prot1 <- sapply(sml.results$prot1,
                            function(uniprot){
                                unique(subset(uniprot.map, uniprot==uniprot)$hgnc)
                            })
sml.results$prot2 <- sapply(sml.results$prot2,
                            function(uniprot){
                                unique(subset(uniprot.map, uniprot==uniprot)$hgnc)
                            })

write.table(sml.results.final[order(sml.results.final$prot1),],
            "out/go-sem-sim.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")
