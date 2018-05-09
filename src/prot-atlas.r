library(tidyr)

kinome <- read.table("data/human-kinome.txt", as.is=TRUE)$V1

tissue <- read.delim("data/external/protein-atlas-normal-tissue.tsv", as.is=TRUE)
tissue.kin <- subset(tissue, Gene.name %in% kinome & Reliability=="Approved",
                     select=c("Gene.name", "Tissue", "Cell.type"))

tissue.kin$tis.cell <- paste(tissue.kin$Tissue, tissue.kin$Cell.type, sep="-")

tissue.kin.table <- table(tissue.kin$Gene.name, tissue.kin$tis.cell)
## tissue.kin.cor <- cor(t(tissue.kin.table))
## tissue.kin.cor.df <- melt(tissue.kin.cor)
## colnames(tissue.kin.cor.df) <- c("node1", "node2", "tissue.cor")

tissue.kin.dist <- dist(tissue.kin.table, method="manhattan")
tissue.kin.dist.df <- melt(as.matrix(tissue.kin.dist))
## tissue.kin.dist.df$value <- 1.0-tissue.kin.dist.df$value
tissue.kin.dist.df$value <- 1.0-(tissue.kin.dist.df$value/ncol(tissue.kin.table))
colnames(tissue.kin.dist.df) <- c("node1", "node2", "tissue.sim")

write.table(tissue.kin.dist.df, "out/protein-atlas-tissue-sim.tsv",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


subcell <- read.delim("data/external/protein-atlas-subcellular-loc.tsv", as.is=TRUE)
subcell.kin.dirty <- subset(subcell, Gene.name %in% kinome & Reliability=="Approved",
                     select=c("Gene.name", "Approved"))

subcell.kin <- subcell.kin.dirty %>% separate_rows(Approved, sep=";")

subcell.kin.table <- table(subcell.kin$Gene.name, subcell.kin$Approved)

subcell.kin.dist <- dist(subcell.kin.table, method="manhattan")
subcell.kin.dist.df <- melt(as.matrix(subcell.kin.dist))
## subcell.kin.dist.df$value <- 1.0-subcell.kin.dist.df$value
subcell.kin.dist.df$value <- 1.0-(subcell.kin.dist.df$value/ncol(subcell.kin.table))
colnames(subcell.kin.dist.df) <- c("node1", "node2", "subcell.sim")


write.table(subcell.kin.dist.df, "out/protein-atlas-subcell-sim.tsv",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

