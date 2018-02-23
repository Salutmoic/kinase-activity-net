suppressMessages(library(Biobase))
library(entropy)
## for the get.bad.inhib.conds function
source("src/optimize-activity-table.r")

argv <- commandArgs(TRUE)
if (length(argv) != 4){
    stop("USAGE: <script> KINACT_PREDS KINASE_LIST MIN_SITES OUT_FILE")
}

kinact.file <- argv[1]
kinase.list.file <- argv[2]
min.sites <- as.integer(argv[3])
out.file <- argv[4]

kinase.list <- read.table(kinase.list.file, as.is=TRUE)$V1

kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                as.is=TRUE)
kinase.conditions <- subset(kinase.conditions, Regulation=="down")

## Load data
load(kinact.file)
kin.act <- log.wKSEA.kinase_condition.clean
load("data/external/esetNR.Rdata")
cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals <- exprs(esetNR)
kin.sub.tbl <- read.delim("data/psiteplus-kinase-substrates.tsv", as.is=TRUE, sep="\t")
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

## Filter out conditions where a supposedly inhibited kinase has
## positive or unchanged activity
bad.inhib.conds <- get.bad.inhib.conds(kin.act)
kin.act <- kin.act[,setdiff(colnames(kin.act), bad.inhib.conds)]

## Filter out kinases that have overall low activity predictions.
kin.act.filt <- filter.low.activity(kin.act)

## kin.act[kin.act>-1 & kin.act<1] <- 1

kin1s <- c()
kin2s <- c()
scores <- c()

for (kin1 in unique(kinase.conditions$Kinase)){
    if (!(kin1 %in% row.names(kin.act.filt)))
        next
    kin1.conds <- subset(kinase.conditions, Kinase==kin1)$Condition
    kin1.conds <- kin1.conds[which(kin1.conds %in% colnames(kin.act.filt))]
    kin1.acts <- kin.act.filt[kin1, kin1.conds]
    kin1.sign <- sign(kin1.acts)
    for (kin2 in kinase.list){
        if (!(kin2 %in% rownames(kin.act.filt))){
            next
        }
        kin2.conds <- subset(kinase.conditions, Kinase==kin2)$Condition
        kin1.conds.tmp <- setdiff(kin1.conds, kin2.conds)
        if (length(kin1.conds.tmp) == 0){
            next
        }
        kin2.acts <- kin.act.filt[kin2, kin1.conds.tmp]
        kin2.sign <- sign(kin2.acts)
        if (!all(is.na(kin1.acts)) && !all(is.na(kin2.acts))){
            ## score <-  (cor(kin1.acts, kin2.acts, use="pairwise.complete.obs", method="spearman")+1)/2
            matches <- as.integer(kin2.sign == kin1.sign)
            matches <- na.omit(matches)
            matches <- sort(matches, decreasing=TRUE)
            score <- sum(matches/log2((1:length(matches))+1))
        }else{
            score <- NA
        }
        scores <- c(scores, score)
        kin1s <- c(kin1s, kin1)
        kin2s <- c(kin2s, kin2)
    }
}

scores <-  scores/max(scores, na.rm=TRUE)

results <- data.frame(node1=kin1s, node2=kin2s, score=scores)
results <- results[order(results$node1),]

write.table(results, out.file, sep="\t", quote=FALSE, col.names=TRUE,
            row.names=FALSE)
