library(reshape2)
suppressMessages(library(Biobase))
source("src/optimize-activity-table.r")

argv <- commandArgs(TRUE)
if (length(argv) != 1){
    stop("USAGE: <script> KINACT_PREDS")
}
kinact.file <- argv[1]

load("data/external/esetNR.Rdata")
load(kinact.file)
kin.act <- log.wKSEA.kinase_condition.clean
cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals <- exprs(esetNR)
kin.sub.tbl <- read.delim("data/psiteplus-kinase-substrates.tsv", as.is=TRUE,
                          sep="\t")
kin.sub.tbl$SUB_MOD_RSD <- sub("[STY]", "", kin.sub.tbl$SUB_MOD_RSD)
kin.overlap <- read.delim("data/psiteplus-kinase-substrate-overlap.tsv",
                          as.is=TRUE)

## Please note also we added the breast cancer CPTAC data from two different
## sources: The CPTAC database as well as from the supplementary material of
## the paper (which got published later). I would recommend you to use the
## ones from the supplementary material of the paper marked as ³supp² and
## avoid the ones from CPTAC database (³1604_290² to ³1711_290²).
cptac.start <- 1604
cptac.end <- 1711
cptac.range <- 1604:1711
cptac.conds <- paste(as.character(cptac.range), "290", sep="_")
phospho.vals <- phospho.vals[, -which(colnames(phospho.vals) %in% cptac.conds)]

min.sites <- 3

## TODO: remove?
## Only choose kinases that have at least 10 substrates
## kin.overlap.sub <- subset(kin.overlap, no.substrates1 >= min.sites & no.substrates2 >= min.sites)
kin.overlap.sub <- kin.overlap
good.kins <- unique(c(kin.overlap.sub$kinase1, kin.overlap.sub$kinase2))

## Remove redundant kinases
redundant.kins <- get.redundant.kins(kin.overlap.sub, kin.act)

kin.act <- kin.act[setdiff(good.kins, redundant.kins),]

## Filter out kinase activity predictions that were not based on
## enough sites
kin.act <- filter.act.preds(kin.act, phospho.anno, phospho.vals, kin.sub.tbl,
                            min.sites)

## ## We look for "empty" rows here...not really empty, but anything
## ## that's practically empty (>90% missing values).
## empty.rows <- rownames(kin.act)[which(apply(kin.act, 1,
##                                             function(tbl.row){
##                                                 percent.na(tbl.row)>0.9
##                                             }))]
## kin.act <- kin.act[setdiff(rownames(kin.act), empty.rows),]
## empty.cols <- colnames(kin.act)[which(apply(kin.act, 2,
##                                             function(tbl.col){
##                                                 percent.na(tbl.col)>0.9
##                                             }))]
## kin.act <- kin.act[,setdiff(colnames(kin.act), empty.cols)]

## Filter out conditions where a supposedly inhibited kinase has
## positive or unchanged activity
bad.inhib.conds <- get.bad.inhib.conds(kin.act)
kin.act <- kin.act[,setdiff(colnames(kin.act), bad.inhib.conds)]

## Filter out kinases that have overall low activity predictions.
kin.act.filt <- filter.low.activity(kin.act)

kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                as.is=TRUE)

all.kins <- rownames(kin.act)
kin.pairs <- expand.grid(all.kins, all.kins, stringsAsFactors=FALSE)

kin1s <- c()
kin2s <- c()
p.vals <- c()
cors <- c()
for (i in 1:nrow(kin.pairs)){
    kin1 <- kin.pairs[i, 1]
    kin2 <- kin.pairs[i, 2]
    if (kin1 == kin2 || !(kin1 %in% rownames(kin.act)) ||
        !(kin2 %in% rownames(kin.act))){
        kin1s <- c(kin1s, kin1)
        kin2s <- c(kin2s, kin2)
        p.vals <- c(p.vals, NA)
        cors <- c(cors, NA)
        next
    }
    kin2.conds <- kinase.conditions[which(kinase.conditions$Kinase==kin2),
                                    "Condition"]
    kin2.conds <- intersect(kin2.conds, colnames(kin.act))
    kin1.act <- kin.act[kin1,]
    kin2.act <- kin.act[kin2,]
    ## Remove conditions where kinase 2 was inhibited or
    ## activated, so we don't calculate correlations on them
    kin2.act[kin2.conds] <- NA
    kin1.meas.conds <- which(!is.na(kin1.act))
    kin2.meas.conds <- which(!is.na(kin2.act))
    shared.conds <- intersect(kin1.meas.conds, kin2.meas.conds)
    if (length(shared.conds) < 5){
        p.value <- NA
        cor.est <- NA
    }else{
        ct <- cor.test(kin1.act, kin2.act, method="spearman",
                       use="pairwise.complete.obs")
        p.value <- -log10(ct$p.value)
        ## z-transformation of Spearman's rho
        cor.est <- sqrt((length(shared.conds)-3)/1.06)*atanh(ct$estimate)
        if (is.infinite(p.value)){
            p.value <- NA
            cor.est <- NA
        }
    }
    kin1s <- c(kin1s, kin1)
    kin2s <- c(kin2s, kin2)
    p.vals <- c(p.vals, p.value)
    cors <- c(cors, cor.est)
}

cors[is.infinite(cors) & cors < 0] <- min(cors[!is.infinite(cors) & !is.na(cors)])
cors[is.infinite(cors) & cors > 0] <- max(cors[!is.infinite(cors) & !is.na(cors)])

max.cor.est <- max(cors, na.rm=TRUE)
min.cor.est <- min(cors, na.rm=TRUE)
print(c(min.cor.est, max.cor.est))
print(length(which(!is.na(cors))))
cors <- (cors-min.cor.est)/(max.cor.est-min.cor.est)

results <- data.frame(node1=kin1s, node2=kin2s, kinact.cor.p=p.vals)
results$kinact.cor.p <- results$kinact.cor.p/max(results$kinact.cor.p, na.rm=TRUE)
write.table(results[order(results$node1),], "out/kinact-full-cor.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

cor.est.results <- data.frame(node1=kin1s, node2=kin2s, kinact.cor.est=cors)
write.table(cor.est.results[order(cor.est.results$node1),],
            "out/kinact-full-cor-sign.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

## ## kin.cor <- acast(results, node1 ~ node2)
## kin.cor <- matrix(p.vals*sign(cors), nrow=length(all.kins))
## rownames(kin.cor) <- all.kins
## colnames(kin.cor) <- all.kins

## kin1s <- c()
## kin2s <- c()
## p.vals <- c()
## for (i in 1:nrow(kin.pairs)){
##     kin1 <- kin.pairs[i, 1]
##     kin2 <- kin.pairs[i, 2]
##     kin1.cor.pairs <- which(!is.na(kin.cor[kin1,]))
##     kin2.cor.pairs <- which(!is.na(kin.cor[kin2,]))
##     shared.pairs <- intersect(kin1.cor.pairs, kin2.cor.pairs)
##     if (length(shared.pairs) < 3){
##         p.value <- NA
##     }else{
##         kin1.cor <- kin.cor[kin1,]
##         kin2.cor <- kin.cor[kin2,]
##         ct <- cor.test(kin1.cor, kin2.cor, method="spearman",
##                        use="pairwise.complete.obs")
##         p.value <- -log10(ct$p.value)
##         if (is.infinite(p.value))
##             p.value <- NA
##     }
##     kin1s <- c(kin1s, kin1)
##     kin2s <- c(kin2s, kin2)
##     p.vals <- c(p.vals, p.value)
## }

## results <- data.frame(node1=kin1s, node2=kin2s, kinact.cor2.p=p.vals)
## results$kinact.cor2.p <- results$kinact.cor2.p/max(results$kinact.cor2.p, na.rm=TRUE)

## write.table(results[order(results$node1),], "out/kinact-full-cor2.tsv", quote=FALSE,
##             row.names=FALSE, col.names=TRUE, sep="\t")
