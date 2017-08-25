library(FunChisq)
library(reshape2)

do.nfchisq <- function(kin1, kin2, kin.act){
    if (kin1 == kin2)
        return(NA)
    cont.tbl <- table(kin.act[kin1,], kin.act[kin2,])
    r <- fun.chisq.test(cont.tbl, method="nfchisq")
    return(r$statistic)
}

v.do.nfchisq <- Vectorize(do.nfchisq, vectorize.args=c("kin1", "kin2"))

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> DATA_FILE OUT_FILE")
}

data.file <- argv[1]
out.file <- argv[2]

kin.act.tbl <- read.table(data.file, as.is=TRUE)
names(kin.act.tbl) <- c("kinase", "condition",  "discr.act")
kin.act.m <- acast(kin.act.tbl, kinase ~ condition)

nfchisq <- outer(rownames(kin.act.m), rownames(kin.act.m),
                 v.do.nfchisq, kin.act.m)
rownames(nfchisq) <- rownames(kin.act.m)
colnames(nfchisq) <- rownames(kin.act.m)

nfchisq.tbl <- melt(nfchisq)
nfchisq.tbl <- nfchisq.tbl[order(nfchisq.tbl[,1], nfchisq.tbl[,2]),]
write.table(nfchisq.tbl, out.file, row.names=FALSE, col.names=FALSE,
            sep="\t", quote=FALSE)
