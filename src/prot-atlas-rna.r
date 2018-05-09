library(reshape2)
library(e1071)

rna.tissue <- read.delim("data/rna_tissue.tsv", as.is=TRUE)
rna.tissue.m <- acast(rna.tissue, kinase ~ tissue)
rna.tissue.sel <- apply(rna.tissue.m, 1, skewness, na.rm=TRUE)
max.sel <- max(rna.tissue.sel, na.rm=TRUE)
min.sel <- min(rna.tissue.sel, na.rm=TRUE)
rna.tissue.sel <- (rna.tissue.sel-min.sel)/(max.sel-min.sel)

kinases <- names(rna.tissue.sel)
kin.pairs <- expand.grid(kinases, kinases, stringsAsFactors=FALSE)
tissue.sel.df <- cbind(kin.pairs,
                       apply(kin.pairs, 1, function(pair) rna.tissue.sel[pair[1]]),
                       apply(kin.pairs, 1, function(pair) rna.tissue.sel[pair[2]]))
names(tissue.sel.df) <- c("node1", "node2", "tissue.sel1", "tissue.sel2")

write.table(tissue.sel.df, "out/rna-tissue-selectivity.tsv", sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)

rna.celline <- read.delim("data/rna_celline.tsv", as.is=TRUE)
rna.celline.m <- acast(rna.celline, kinase ~ celline)
rna.celline.sel <- apply(rna.celline.m, 1, skewness, na.rm=TRUE)
max.sel <- max(rna.celline.sel, na.rm=TRUE)
min.sel <- min(rna.celline.sel, na.rm=TRUE)
rna.celline.sel <- (rna.celline.sel-min.sel)/(max.sel-min.sel)

kinases <- names(rna.celline.sel)
kin.pairs <- expand.grid(kinases, kinases, stringsAsFactors=FALSE)
celline.sel.df <- cbind(kin.pairs,
                       apply(kin.pairs, 1, function(pair) rna.celline.sel[pair[1]]),
                       apply(kin.pairs, 1, function(pair) rna.celline.sel[pair[2]]))
names(celline.sel.df) <- c("node1", "node2", "celline.sel1", "celline.sel2")

write.table(celline.sel.df, "out/rna-celline-selectivity.tsv", sep="\t", quote=FALSE,
            row.names=FALSE, col.names=TRUE)
