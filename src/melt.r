library(reshape2)

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> TABLE OUT_FILE")
}

data.tbl.file <- argv[1]
out.file <- argv[2]

data.tbl <- read.delim(data.tbl.file, as.is=TRUE)
data.tbl.m <- as.matrix(data.tbl[,2:ncol(data.tbl)])
colnames(data.tbl.m) <- colnames(data.tbl)[2:ncol(data.tbl)]
rownames(data.tbl.m) <- data.tbl[,1]

data.tbl.df <- melt(data.tbl.m)

write.table(data.tbl.df[order(data.tbl.df[,1], data.tbl.df[,2]),],
            out.file, col.names=FALSE, row.names=FALSE,
            quote=FALSE, sep="\t")
