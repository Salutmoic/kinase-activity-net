suppressMessages(library(Biobase))
source("src/optimize-activity-table.r")

argv <- commandArgs(TRUE)
if (length(argv) != 4){
    stop("USAGE: <script> STRATEGY MIN_SITES ENTROPY_FILTER NA_THRESHOLD")
}
strategy <- argv[1]
min.sites <- as.integer(argv[2])
entropy.filter.stren <- as.numeric(argv[3])
na.threshold <- as.numeric(argv[4])

## General parameters
entropy.bins <- 10

## Load data
load("data/external/esetNR.Rdata")
load("data/log.wKSEA.kinase_condition.clean.Rdata")
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

## TODO: remove?
## Only choose kinases that have at least 10 substrates
## kin.overlap.sub <- subset(kin.overlap, no.substrates1 >= min.sites & no.substrates2 >= min.sites)
kin.overlap.sub <- kin.overlap
good.kins <- unique(c(kin.overlap.sub$kinase1, kin.overlap.sub$kinase2))

## Remove redundant kinases
redundant.kins <- get.redundant.kins(kin.overlap.sub)

kin.act <- kin.act[setdiff(good.kins, redundant.kins),]

## Filter out kinase activity predictions that were not based on
## enough sites
kin.act <- filter.act.preds(kin.act, phospho.anno, phospho.vals, kin.sub.tbl,
                            min.sites)

## We look for "empty" rows here...not really empty, but anything
## that's practically empty (>90% missing values).
empty.rows <- rownames(kin.act)[which(apply(kin.act, 1,
                                            function(tbl.row){
                                                percent.na(tbl.row)>0.9
                                            }))]
kin.act <- kin.act[setdiff(rownames(kin.act), empty.rows),]
empty.cols <- colnames(kin.act)[which(apply(kin.act, 2,
                                            function(tbl.col){
                                                percent.na(tbl.col)>0.9
                                            }))]
kin.act <- kin.act[,setdiff(colnames(kin.act), empty.cols)]

## Filter out conditions where a supposedly inhibited kinase has
## positive or unchanged activity
bad.inhib.conds <- get.bad.inhib.conds(kin.act)
kin.act <- kin.act[,setdiff(colnames(kin.act), bad.inhib.conds)]

## Filter out kinases that have overall low activity predictions.
kin.act.filt <- filter.low.activity(kin.act)

## "CDK2", "CDK4",
egf.kinases <- c("AKT1", "AKT2", "ARAF", "BRAF", "CHUK", "GSK3A",
                 "GSK3B", "JAK2", "MAP2K1", "MAP2K2", "MAPK1",
                 "MAPK3", "MTOR", "PAK1", "PDPK1", "PRKAA2", "PRKACA",
                 "PTK2", "RAF1", "RPS6KB1", "RPS6KB2", "SGK1", "SRC")
##### Only kinases w/ perturbations
## egf.kinases <- c("AKT1", "BRAF", "MAP2K1", "MAPK1", "MTOR", "PDPK1", "PRKACA",
##                  "RPS6KB1", "SRC")
kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                as.is=TRUE)

## kinase.invivoconditions <- read.delim("data/external/kinase_invivoconditions.csv",
##                                       as.is=TRUE, sep=",")
## names(kinase.invivoconditions) <- c("Condition", "Kinase", "Description",
##                                     "invivoAct", "pAct")

kinase.cond <- subset(kinase.conditions, Kinase %in% egf.kinases &
                                         Kinase %in% rownames(kin.act.filt) &
                                         Condition %in% colnames(kin.act.filt) &
                                         Control == "Control" &
                                         Regulation=="down")
kinase.cond$activity <- apply(kinase.cond, 1,
                              function(row){
                                  kinase <- row["Kinase"]
                                  cond <- row["Condition"]
                                  return(kin.act.filt[kinase, cond])
                              })
kinase.cond <- subset(kinase.cond, activity < 0)

## kinase.cond.best <- c()
## for (kin in unique(kinase.cond$Kinase)){
##     kinase.cond.sub <- subset(kinase.cond, Kinase==kin)
##     act.order <- order(order(kinase.cond.sub$activity))
##     cond.nas <- c()
##     for (cond in kinase.cond.sub$Condition){
##         cond.na <- percent.na(kin.act.filt[,cond])
##         cond.nas <- c(cond.nas, cond.na)
##     }
##     na.order <- order(order(cond.nas))
##     best.cond <- which.min(abs(act.order-na.order))
##     kinase.cond.best <- rbind(kinase.cond.best, kinase.cond.sub[best.cond,])
## }
## names(kinase.cond.best) <- names(kinase.cond)
## kinase.cond <- kinase.cond.best

## kinase.invivocond <- subset(kinase.invivoconditions, Kinase %in% egf.kinases)

egf.kin.act <- kin.act.filt[rownames(kin.act.filt) %in% kinase.cond$Kinase, ]


egf.kin.act.pert.full <- egf.kin.act[, unique(kinase.cond$Condition)]
egf.kin.act.pert <- egf.kin.act.pert.full[which(apply(egf.kin.act.pert.full,
                                                      1, percent.na) < 0.5),]
egf.kin.act.pert <- egf.kin.act.pert[,which(apply(egf.kin.act.pert,
                                                  2, percent.na) < 1/3)]
## lost.conds <- which(!(kinase.cond$Condition %in% colnames(egf.kin.act.pert)))
## lost.kins <- unique(kinase.cond$Kinase[lost.conds])
## kept.kins <- setdiff(rownames(egf.kin.act.pert), lost.kins)
## egf.kin.act.pert <- egf.kin.act.pert[kept.kins,]

egf.kin.act.all <- egf.kin.act[, which(apply(egf.kin.act, 2, percent.na) < 1/3)]

non.pert.conds <- setdiff(colnames(kin.act.filt), unique(kinase.cond$Condition))

## Build the optimized table according to the desired strategy
if (strategy=="balanced"){
    max.unpert.table <- make.balanced.table(subset(kin.act.filt,
                                                   select=non.pert.conds),
                                            na.threshold,
                                            precious.rows=rownames(egf.kin.act.pert))
}else if (strategy=="max-rows"){
    max.unpert.table <- make.max.row.table(subset(kin.act.filt,
                                                  select=non.pert.conds),
                                           na.threshold, 50,
                                           precious.rows=rownames(egf.kin.act.pert))
}else if (strategy=="max-cols"){
    max.unpert.table <- make.max.col.table(subset(kin.act.filt,
                                                  select=non.pert.conds),
                                           na.threshold, 20,
                                           precious.rows=rownames(egf.kin.act.pert))
}

if (is.null(max.unpert.table))
    stop("Failed to optimize table.")

max.unpert.table.filt <- filter.low.activity(max.unpert.table)

## Remove any kinases or conditions that do not show enough
## variability in their activity predictions.
## max.unpert.table.filt <- filter.by.entropy(max.unpert.table, entropy.filter.stren,
##                                            entropy.bins)


## Output the tble
print.table.info(egf.kin.act.pert)
write.activity.table(egf.kin.act.pert, "data/kinact-egf-pert.tsv")
print.table.info(egf.kin.act.all)
write.activity.table(egf.kin.act.all, "data/kinact-egf-all.tsv")
print.table.info(max.unpert.table.filt)
write.activity.table(max.unpert.table.filt,
                     paste0("data/kinact-egf-unpert-", strategy, ".tsv"))

## Do imputation and output the imputed table
egf.kin.act.pert.imp <- impute.on.table(egf.kin.act.pert)
egf.kin.act.pert.imp.df <- melt(egf.kin.act.pert.imp)
names(egf.kin.act.pert.imp.df) <- c("kinase", "condition", "activity")
kinase.cond.pairs <- expand.grid(rownames(egf.kin.act.pert.imp),
                                 colnames(egf.kin.act.pert.imp)[1:(ncol(egf.kin.act.pert.imp)-1)])
egf.kin.act.pert.imp.df$is.pert <- unlist(apply(kinase.cond.pairs, 1,
                                              function(pair){
                                                  cond.kinases <- subset(kinase.cond,
                                                                         Condition==pair[2])$Kinase
                                                  return (pair[1] %in% cond.kinases)
                                              }))
egf.kin.act.pert.imp.df <- egf.kin.act.pert.imp.df[order(egf.kin.act.pert.imp.df$kinase,
                                                     egf.kin.act.pert.imp.df$condition),]
write.table(egf.kin.act.pert.imp.df, "data/kinact-egf-pert-imp.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

max.unpert.table.imp <- impute.on.table(max.unpert.table.filt)
write.activity.table(max.unpert.table.imp,
                     paste0("data/kinact-egf-unpert-", strategy, "-imp.tsv"))

egf.kin.act.all.imp <- impute.on.table(egf.kin.act.all)
write.activity.table(egf.kin.act.all.imp, "data/kinact-egf-all-imp.tsv")
