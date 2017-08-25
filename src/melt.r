library(reshape2)
source("src/optimize-activity-table.r")

argv <- commandArgs(TRUE)
if (!(length(argv) %in% c(2, 3))){
    stop("USAGE: <script> TABLE OUT_FILE [NORMALIZE]")
}

data.tbl.file <- argv[1]
out.file <- argv[2]
if (length(argv) == 3){
    library(preprocessCore)
    norm <- TRUE
}else{
    norm <- FALSE
}

data.tbl <- read.delim(data.tbl.file, as.is=TRUE)
if (length(data.tbl) == 1)
    data.tbl <- read.delim(data.tbl.file, as.is=TRUE, sep=",")
data.tbl.m <- as.matrix(data.tbl[,2:ncol(data.tbl)])
colnames(data.tbl.m) <- colnames(data.tbl)[2:ncol(data.tbl)]
rownames(data.tbl.m) <- data.tbl[,1]

if (norm){
    data.tbl.m.norm <- normalize.quantiles(data.tbl.m)
    ## data.tbl.m.norm <- apply(data.tbl.m, 2, function(col) col-median(col, na.rm=TRUE))
    colnames(data.tbl.m.norm) <- colnames(data.tbl.m)
    rownames(data.tbl.m.norm) <- rownames(data.tbl.m)
    data.tbl.m <- data.tbl.m.norm
}

data.tbl.df <- melt(data.tbl.m)

write.table(data.tbl.df[order(data.tbl.df[,1], data.tbl.df[,2]),],
            out.file, col.names=FALSE, row.names=FALSE,
            quote=FALSE, sep="\t")
