suppressMessages(library(Biobase))
library(preprocessCore)
source("src/optimize-activity-table.r")

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> KINACT_PREDS NETWORK CLUSTS")
}
kinact.file <- argv[1]
network.file <- argv[2]
clusts.file <- argv[3]

load("data/external/esetNR.Rdata")
cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals <- exprs(esetNR)

cptac.exps <-  c("272")
cptac.conds <- rownames(subset(cond.anno, experiment_id %in% cptac.exps))
phospho.vals <- phospho.vals[, which(colnames(phospho.vals) %in% cptac.conds)]

load(kinact.file)
kin.act <- log.wKSEA.kinase_condition.clean
kin.act <- kin.act[, which(colnames(kin.act) %in% cptac.conds)]
kin.act.norm <- normalize.quantiles(kin.act)
rownames(kin.act.norm) <- rownames(kin.act)
colnames(kin.act.norm) <- colnames(kin.act)
kin.act.norm <- kin.act

network <- read.delim(network.file, as.is=TRUE)
rownames(network) <- paste(network$prot1, network$prot2, sep="-")
clusts <- read.table(clusts.file, as.is=TRUE)
names(clusts) <- c("prot", "clust")

good.kins <- rownames(kin.act)[unique(unlist(apply(kin.act, 2,
                                                   function(col){
                                                       which(!is.na(col))
                                                   })))]

## cluster 2
clust.kins <- subset(clusts, prot %in% good.kins & clust==2)$prot
clust.kin.act.full <- kin.act[clust.kins,]
clust.kin.act <- make.max.row.table(clust.kin.act.full, 0.25, 1)
clust.kins <- rownames(clust.kin.act)
clust.net <- subset(network, prot1 %in% clust.kins & prot2 %in% clust.kins)

calc.coherence <- function(kinact1, kinact2, rel.sign){
    coh <- abs(kinact1-kinact2)*(kinact1*rel.sign)/kinact2
    return(coh)
}

calc.avg.coherence <- function(kin.act, prot1, prot2, rel.sign){
    all.cohs <- c()
    for (i in 1:ncol(kin.act)){
        kinact1 <- kin.act[prot1, i]
        kinact2 <- kin.act[prot2, i]
        all.cohs <- c(all.cohs, calc.coherence(kinact1, kinact2, rel.sign))
    }
    return(median(all.cohs))
}

gen.rand.subnet <- function(network, kin.act){
    prot1s <- c()
    prot2s <- c()
    probs <- c()
    signs <- c()
    cohs <- c()
    for (n in 1:nrow(network)){
        if (network$prot1[n] == network$prot2[n])
            next
        has.edge <- runif(1) < network$bart.pred.mean[n]
        if (!has.edge)
            next
        rel.sign <- as.integer(runif(1) < network$prob.act.mean[n])*2 - 1
        prot1 <- network$prot1[n]
        prot2 <- network$prot2[n]
        prob <- network$bart.pred.mean[n]
        coh <- calc.avg.coherence(kin.act, prot1, prot2, rel.sign)
        prot1s <- c(prot1s, prot1)
        prot2s <- c(prot2s, prot2)
        probs <- c(probs, prob)
        signs <- c(signs, rel.sign)
        cohs <- c(cohs, coh)
    }
    new.net <- data.frame(prot1=prot1s, prot2=prot2s, bart.pred.mean=probs,
                          sign=signs, coherence=cohs)
    rownames(new.net) <- paste(new.net$prot1, new.net$prot2, sep="-")
    return(new.net)
}

add.edge <- function(full.network, test.network, kin.act, row){
    prot1 <- full.network[row, "prot1"]
    prot2 <- full.network[row, "prot2"]
    prob <- full.network[row, "bart.pred.mean"]
    rel.sign <- as.integer(runif(1) < network$prob.act.mean[n])*2 - 1
    coh <- calc.avg.coherence(kin.act, prot1, prot2, rel.sign)
    new.net <- rbind(test.network, data.frame(prot1=prot1, prot2=prot2,
                                              bart.pred.mean=prob,
                                              sign=rel.sign, coherence=coh))
    rownames(new.net)[nrow(new.net)] <- row
    return(new.net)
}

change.sign <- function(test.network, kin.act, row){
    new.net <- test.network
    prot1 <- new.net[row, "prot1"]
    prot2 <- new.net[row, "prot2"]
    rel.sign <- new.net[row, "sign"] * -1
    new.net[row, "sign"] <- rel.sign
    coh <- calc.avg.coherence(kin.act, prot1, prot2, rel.sign)
    new.net[row, "coherence"] <- coh
    return(new.net)
}

remove.edge <- function(test.network, row){
    new.net <- test.network
    keep.edges <- setdiff(rownames(new.net), row)
    return(new.net[keep.edges,])
}

calc.perf <- function(network, full.network){
    coherence <- max(c(0, sum(network$coherence)))
    prob <- prod(network$bart.pred.mean)
    n <- nrow(network)
    return((coherence*prob*n)^(1/3))
}

best.net <- gen.rand.subnet(clust.net, clust.kin.act)
best.perf <- calc.perf(best.net, clust.net)
all.perfs <- c(best.perf)
op <- NULL
for (i in 1:100000){
    rand.pick <- runif(1)
    if (rand.pick < 0.20){
        new.net <- gen.rand.subnet(clust.net, clust.kin.act)
        op <- "newnet"
   }else{
        if (rand.pick < 0.4){
            row <- sample(rownames(clust.net), 1, prob=clust.net$bart.pred.mean)
            op <- "rand.row"
        }else{
            row <- rownames(best.net)[which.min(best.net$coherence)]
            op <- "worst.row"
        }
        if (!(row %in% rownames(best.net))){
            new.net <- add.edge(clust.net, best.net, clust.kin.act, row)
            op <- paste(op, "newedge")
        }else{
            rand.pick2 <- runif(1)
            prob.act <- clust.net[row, "prob.act.mean"]
            if ((best.net[row, "sign"] == 1 && rand.pick2 > prob.act) |
                (best.net[row, "sign"] == -1 && rand.pick2 < prob.act)){
                new.net <- change.sign(best.net, clust.kin.act, row)
                op <- paste(op, "sign")
            }else{
                new.net <- remove.edge(best.net, row)
                op <- paste(op, "deledge")
            }
        }
   }
    if (nrow(new.net) == 0){
        break
    }
    new.perf <- calc.perf(new.net, clust.net)
    if (new.perf > best.perf){
        best.perf <- new.perf
        best.net <- new.net
    }
    all.perfs <- c(all.perfs, best.perf)
}
print(best.perf)
best.net <- merge(best.net, clust.net, by=c("prot1", "prot2"))
rownames(best.net) <- paste(best.net$prot1, best.net$prot2, sep="-")
message("Original network:")
filt.net <- subset(clust.net, bart.pred.mean>0.95)
table(filt.net$assoc.is.known)
table(filt.net$assoc.is.known)/nrow(filt.net)
message("Optimised network:")
table(best.net$assoc.is.known)
table(best.net$assoc.is.known)/nrow(best.net)
pdf("tmp/coherence.pdf")
plot(all.perfs, type="l", ylab="performance", xlab="iteration")
dev.off()
 
