library(reshape2)
library(preprocessCore)

raw.dat <- read.delim("data/gtex-tissue-expression.tsv", as.is=TRUE)

kinome <- read.table("data/human-kinome.txt", as.is=TRUE)$V1
raw.dat <- subset(raw.dat, Gene.Name %in% kinome)

expr <- raw.dat[3:ncol(raw.dat)]

coexpr <- cor(t(expr), method="spearman", use="pairwise.complete.obs")

rownames(coexpr) <- raw.dat$Gene.Name
colnames(coexpr) <- raw.dat$Gene.Name

coexpr.df <- melt(coexpr)
colnames(coexpr.df) <- c("node1", "node2", "tissue.coexpr")

min.cor <- min(coexpr.df$tissue.coexpr, na.rm=TRUE)
max.cor <- max(coexpr.df$tissue.coexpr, na.rm=TRUE)
coexpr.df$tissue.coexpr <- (coexpr.df$tissue.coexpr-min.cor)/(max.cor-min.cor)

write.table(coexpr.df, "out/gtex-tissue-coexpr.tsv", sep="\t",
            row.names=FALSE, col.names=TRUE, quote=FALSE)
