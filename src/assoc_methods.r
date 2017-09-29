library(ppcor)
library(FNN)
library("FunChisq")
library('infotheo')
library(reshape2)

do.fnn.mut.info <- function(kin1, kin2, kin.act){
    if (kin1 == kin2)
        return(NA)
    r <- mutinfo(kin.act[kin1,], kin.act[kin2,])
    
    return(r)
}


do.nfchisq <- function(kin1, kin2, kin.act){
    if (kin1 == kin2)
        return(NA)
    cont.tbl <- table(kin.act[kin1,], kin.act[kin2,])
    r <- fun.chisq.test(cont.tbl, method="nfchisq")
    return(r$statistic)
}


do.mut.info <- function(kin1, kin2, kin.act){
    if (kin1 == kin2)
        return(NA)
    r <- mutinformation(kin.act[kin1,], kin.act[kin2,])
    
    return(r)
}


argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> DATA_FILE OUT_FILE METHOD")
}

data.file <- argv[1]
out.file <- argv[2]
method <- argv[3]


kin.act.tbl <- read.table(data.file, as.is=TRUE)
names(kin.act.tbl) <- c("kinase", "condition",  "discr.act")
kin.act.m <- acast(kin.act.tbl, kinase ~ condition)



if(method == "mut_info"){
v.do.mut.info <- Vectorize(do.mut.info, vectorize.args=c("kin1", "kin2"))


res.tbl <- outer(rownames(kin.act.m), rownames(kin.act.m),v.do.mut.info, kin.act.m)


}


if(method == "fnn_mut_info"){
v.do.fnn.mut.info <- Vectorize(do.fnn.mut.info, vectorize.args=c("kin1", "kin2"))
res.tbl <- outer(rownames(kin.act.m), rownames(kin.act.m),v.do.fnn.mut.info, kin.act.m)
}

if(method == "nfchisq"){
v.do.nfchisq <- Vectorize(do.nfchisq, vectorize.args=c("kin1", "kin2"))
res.tbl<- outer(rownames(kin.act.m), rownames(kin.act.m),v.do.nfchisq, kin.act.m)

}

if(method == "partcor"){
res.tbl = pcor(t(kin.act.m),method = "s")$estimate

}


rownames(res.tbl) <- rownames(kin.act.m)
colnames(res.tbl) <- rownames(kin.act.m)

tbl <- melt(res.tbl)
tbl <- tbl[order(tbl[,1], tbl[,2]),]
write.table(tbl, out.file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)



