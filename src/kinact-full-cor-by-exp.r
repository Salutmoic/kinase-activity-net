library(reshape2)
suppressMessages(library(Biobase))
source("src/optimize-activity-table.r")

normalize.p.vals <- function(cor.p){
    min.p <- min(cor.p[which(cor.p>0)], na.rm=TRUE)
    cor.p[which(cor.p==0)] <- min.p
    cor.p <- log(cor.p)
    min.p <- min(cor.p, na.rm=TRUE)
    max.p <- max(cor.p, na.rm=TRUE)
    cor.p <- (cor.p-min.p)/(max.p-min.p)
    return (cor.p)
}

min.sites <- 3
min.conds <- 5

load("data/external/esetNR.Rdata")
cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals <- exprs(esetNR)

## Please note also we added the breast cancer CPTAC data from two different
## sources: The CPTAC database as well as from the supplementary material of
## the paper (which got published later). I would recommend you to use the
## ones from the supplementary material of the paper marked as ³supp² and
## avoid the ones from CPTAC database (³1604_290² to ³1711_290²).
cptac.exps <-  c("290", "292")
cptac.conds <- rownames(subset(cond.anno, experiment_id %in% cptac.exps))
phospho.vals <- phospho.vals[, -which(colnames(phospho.vals) %in% cptac.conds)]

exp.cond.table <- table(cond.anno$experiment_id)
good.exps <- sort(names(exp.cond.table)[which(exp.cond.table>=min.conds)])
good.exps <- setdiff(good.exps, cptac.exps)

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop(paste0("USAGE: <script> KINACT_PREDS EXP_NUM (EXP_NUM=[1-",
                length(good.exps), "])"))
}
kinact.file <- argv[1]
exp.num <- as.integer(argv[2])
if (exp.num <= 0 || exp.num > length(good.exps))
    stop(paste0("USAGE: <script> KINACT_PREDS EXP_NUM (EXP_NUM=[1-",
                length(good.exps), "])"))

load(kinact.file)
kin.act <- log.wKSEA.kinase_condition.clean
kin.sub.tbl <- read.delim("data/psiteplus-kinase-substrates.tsv", as.is=TRUE,
                          sep="\t")
kin.sub.tbl$SUB_MOD_RSD <- sub("[STY]", "", kin.sub.tbl$SUB_MOD_RSD)
kin.overlap <- read.delim("data/psiteplus-kinase-substrate-overlap.tsv",
                          as.is=TRUE)

experiment.id <- good.exps[exp.num]
exp.conds <- subset(cond.anno, experiment_id==experiment.id)
good.conds <- intersect(rownames(exp.conds), colnames(kin.act))
kin.act.exp <- kin.act[, good.conds]
phospho.vals.exp <- phospho.vals[, good.conds]

## TODO: remove?
## Only choose kinases that have at least 10 substrates
## kin.overlap.sub <- subset(kin.overlap, no.substrates1 >= min.sites & no.substrates2 >= min.sites)
kin.overlap.sub <- kin.overlap
good.kins <- unique(c(kin.overlap.sub$kinase1, kin.overlap.sub$kinase2))

## Remove redundant kinases
redundant.kins <- get.redundant.kins(kin.overlap.sub, kin.act.exp)

kin.act.exp <- kin.act.exp[setdiff(good.kins, redundant.kins),]

## Filter out kinase activity predictions that were not based on
## enough sites
kin.act.exp <- filter.act.preds(kin.act.exp, phospho.anno, phospho.vals.exp,
                                kin.sub.tbl, min.sites)

## Filter out conditions where a supposedly inhibited kinase has
## positive or unchanged activity
bad.inhib.conds <- get.bad.inhib.conds(kin.act.exp)
kin.act.exp <- kin.act.exp[,setdiff(colnames(kin.act.exp), bad.inhib.conds)]

## Filter out kinases that have overall low activity predictions.
## kin.act.filt <- filter.low.activity(kin.act.exp)
kin.act.filt <- kin.act.exp

kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                as.is=TRUE)

all.kins <- rownames(kin.act.exp)
kin.pairs <- expand.grid(all.kins, all.kins, stringsAsFactors=FALSE)

kin1s <- c()
kin2s <- c()
p.vals <- c()
cors <- c()
for (i in 1:nrow(kin.pairs)){
    kin1 <- kin.pairs[i, 1]
    kin2 <- kin.pairs[i, 2]
    if (kin1 == kin2 || !(kin1 %in% rownames(kin.act.filt)) ||
        !(kin2 %in% rownames(kin.act.filt))){
        kin1s <- c(kin1s, kin1)
        kin2s <- c(kin2s, kin2)
        p.vals <- c(p.vals, NA)
        cors <- c(cors, NA)
        next
    }
    kin2.conds <- kinase.conditions[which(kinase.conditions$Kinase==kin2),
                                    "Condition"]
    kin2.conds <- intersect(kin2.conds, colnames(kin.act.filt))
    kin1.act <- kin.act.filt[kin1,]
    kin2.act <- kin.act.filt[kin2,]
    ## Remove conditions where kinase 2 was inhibited or
    ## activated, so we don't calculate correlations on them
    kin2.act[kin2.conds] <- NA
    kin1.meas.conds <- which(!is.na(kin1.act))
    kin2.meas.conds <- which(!is.na(kin2.act))
    shared.conds <- intersect(kin1.meas.conds, kin2.meas.conds)
    if (length(shared.conds) < min.conds){
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

results <- data.frame(node1=kin1s, node2=kin2s, kinact.cor.p=p.vals)
results$kinact.cor.p <- normalize.p.vals(results$kinact.cor.p)
if (all(is.na(results$kinact.cor.p))){
    stop("No results")
}
colnames(results)[ncol(results)] <- paste0("kinact.cor.p.", experiment.id)
write.table(results[order(results$node1),],
            paste0("out/kinact-full-cor-by-exp/kinact-full-cor-exp-",
                   experiment.id, ".tsv"),
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

cor.est.results <- data.frame(node1=kin1s, node2=kin2s, kinact.cor.est=cors)
colnames(cor.est.results)[ncol(cor.est.results)] <- paste0("kinact.cor.est.",
                                                           experiment.id)
write.table(cor.est.results[order(cor.est.results$node1),],
            paste0("out/kinact-full-cor-by-exp/kinact-full-cor-sign-exp-",
                   experiment.id, ".tsv"),
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
