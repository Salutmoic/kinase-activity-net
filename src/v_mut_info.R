library('infotheo')
library(reshape2)

do.mut.info <- function(kin1, kin2, kin.act){
    if (kin1 == kin2)
        return(NA)
    r <- mutinformation(kin.act[kin1,], kin.act[kin2,])
    
    return(r)
}


v.do.mut.info <- Vectorize(do.mut.info, vectorize.args=c("kin1", "kin2"))

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> DATA_FILE OUT_FILE")
}

data.file <- argv[1]
out.file <- argv[2]

kin.act.tbl <- read.table(data.file, as.is=TRUE)
names(kin.act.tbl) <- c("kinase", "condition",  "discr.act")
kin.act.m <- acast(kin.act.tbl, kinase ~ condition)

mut.info <- outer(rownames(kin.act.m), rownames(kin.act.m),
                 v.do.mut.info, kin.act.m)
rownames(mut.info) <- rownames(kin.act.m)
colnames(mut.info) <- rownames(kin.act.m)

mut.tbl <- melt(mut.info)
mut.tbl <- mut.tbl[order(mut.tbl[,1], mut.tbl[,2]),]
write.table(mut.tbl, out.file, row.names=FALSE, col.names=FALSE,
            sep="\t", quote=FALSE)
