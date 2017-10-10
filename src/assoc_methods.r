library(ppcor)
library(FNN)
library("FunChisq")
library('infotheo')
library(reshape2)

average.over.methods <- function(tbls){
    ## here I average over the score of each kinase pair
    ## Funchisq is assymetrical so I average over m[i,j] m[j,i]
    ## for now

    avg.scores = c()

    for(i in 1:(nrow(tbls[[1]])-1)){
        for(j in (i+1):nrow(tbls[[2]])){

            v = 0
            for(tbl in tbls){
                tbl[ is.na(tbl) ] <- NA

                v = v +median(as.numeric(tbl[i,j]),as.numeric(tbl[j,i]),na.rm = T)
            }
            v = v/3

            avg.scores = rbind(avg.scores,c(rownames(tbls[[1]])[i],rownames(tbls[[1]])[j],v))

        }
    }

    return(avg.scores) 
}

normalize.by.max <- function(table){
    ## Normalization method for methods that only return postivie values
    max = max(unlist(table))

    table = table/max

    return(table)
}


normalize.by.range <- function(table){
    ## Normalization method for methods that return postivie and negative values values
    t = unlist(table)
    t = table[!is.na(t)]
    range = range(t)
    min = min(t)
    table = (table-min)/(range[2] - range[1]) 
    return(table)
}

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

if(method == "pcor"){
    res.tbl = pcor(t(kin.act.m),method = "s")$estimate

}

if(method == "all"){

    ## here I do all the methods and normalize them 
    
    methods = list(funch =Vectorize(do.nfchisq, vectorize.args=c("kin1", "kin2")),
                   mutinf = Vectorize(do.fnn.mut.info, vectorize.args=c("kin1", "kin2")))


    tbls = list(abs(pcor(t(kin.act.m),method = "s")$estimate))



    for(i in 1:length(methods)){

        tbls[[i+1]] =  outer(rownames(kin.act.m), rownames(kin.act.m),methods[[i]], kin.act.m)
    }

    for(i in 1:length(tbls)){
        if(i == 1){
            tbls[[i]] = normalize.by.max(tbls[[i]])
        }
        if(i !=1 ){
            tbls[[i]] = normalize.by.range(tbls[[i]])
        }

    }

    names(tbls) = c("partcor","nfchisq","fnn")

    ## here I average over all possible connections 

    rownames(tbls[[1]]) <- rownames(kin.act.m)
    colnames(tbls[[1]]) <- rownames(kin.act.m)

    tbl = average.over.methods(tbls)
    write.table(tbl, out.file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

}

if(method != "all"){
    rownames(res.tbl) <- rownames(kin.act.m)
    colnames(res.tbl) <- rownames(kin.act.m)

    res.tbl[which(is.na(res.tbl))] <- 0.0
    tbl <- melt(res.tbl)
    names(tbl) <- c("node1", "node2", "assoc")
    tbl <- tbl[order(tbl[,1], tbl[,2]),]
    write.table(tbl, out.file, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

}
