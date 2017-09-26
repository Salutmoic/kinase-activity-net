library(reshape2)

mclust.whole <- function(kin.act.m){
## kin.num stores kinase activities
##  All values fed to Mclust as a one dimensional vector
    suppressMessages(library(mclust))
    kin.act.v <- as.vector(kin.act.m)
    clust <- Mclust(kin.act.v)
    discrete = matrix(clust$classification, nrow=nrow(kin.act.m),
                      ncol=ncol(kin.act.m))
    rownames(discrete) <- rownames(kin.act.m)
    colnames(discrete) <- colnames(kin.act.m)
    return(discrete)
}

mclust.by.row <- function(kin.act.m){
## kin.num stores the kinase activites acorss conditions
## discrete contains the discretized values
    suppressMessages(library(mclust))
    disc <- c()
    for(i in 1:nrow(kin.act.m)){
        res <- Mclust(kin.act.m[i,], G=3)$classification
        disc <- rbind(disc, res)
    }
    discrete <- as.matrix(disc)
    rownames(discrete) <- rownames(kin.act.m)
    colnames(discrete) <- colnames(kin.act.m)
    return(discrete)
}

discr.manual <- function(kin.act.m){
## kin.act.m stores kinase activites acorss conditions
    disc <- c()
    for(i in 1:nrow(kin.act.m)){
        row <- rep(0,ncol(kin.act.m))
        p <- kin.act.m[i,]
        row[which(p <= 1 & p >= -1)] = 2
        row[which(p < -1)] = 1
        row[which(p > 1)] = 3
        if(0 %in% row){
            break
        }
        disc <- rbind(disc, row)
    }
    discrete <- as.matrix(disc)
    rownames(discrete) <- rownames(kin.act.m)
    colnames(discrete) <- colnames(kin.act.m)
    return(discrete)
}

discr.trunc <- function(kin.act.m){
## kin.act.m stores kinase activites acorss conditions
    disc <- c()
    for(i in 1:nrow(kin.act.m)){
        p <- kin.act.m[i,]
        row <- trunc(p)
        disc <- rbind(disc, row)
    }
    discrete <- as.matrix(disc)
    rownames(discrete) <- rownames(kin.act.m)
    colnames(discrete) <- colnames(kin.act.m)
    return(discrete)
}

discretize <- function(kin.act.m, method){
    switch(method,
        mclust.whole = mclust.whole(kin.act.m),
        mclust.by.row = mclust.by.row(kin.act.m),
        manual = discr.manual(kin.act.m),
        trunc = discr.trunc(kin.act.m))
}   

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> METHOD DATA_FILE OUT_FILE")
}

method <- argv[1]
data.file <- argv[2]
out.file <- argv[3]

kin.act.tbl <- read.table(data.file, as.is=TRUE)
names(kin.act.tbl) <- c("kinase", "condition", "activity")

kin.act.m <- acast(kin.act.tbl, kinase ~ condition)
discrete <- discretize(kin.act.m, method)

discrete.tbl <- melt(discrete)

write.table(discrete.tbl, out.file, row.names=FALSE,  col.names=FALSE,
            sep="\t", quote=FALSE)
