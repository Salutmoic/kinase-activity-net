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
kin.sub.tbl <- read.delim("data/human_kinase_table.tsv", as.is=TRUE, sep="\t")
kin.sub.tbl$SUB_MOD_RSD <- sub("[STY]", "", kin.sub.tbl$SUB_MOD_RSD)
kin.overlap <- read.delim("data/kinase_substrate_overlap.tsv", as.is=TRUE)

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

## Build the optimized table according to the desired strategy
if (strategy=="balanced"){
    max.table <- make.balanced.table(kin.act.filt, na.threshold)
}else if (strategy=="max-rows"){
    max.table <- make.max.row.table(kin.act.filt, na.threshold, 50)
}else if (strategy=="max-cols"){
    max.table <- make.max.col.table(kin.act.filt, na.threshold, 20)
}
if (is.null(max.table))
    stop("Failed to optimize table.")

## Remove any kinases or conditions that do not show enough
## variability in their activity predictions.
max.table.filt <- filter.by.entropy(max.table, entropy.filter.stren,
                                    entropy.bins)

## Output the tble
print.table.info(max.table.filt)
write.activity.table(max.table.filt, paste0("data/kinact-", strategy, ".tsv"))

## Do imputation and output the imputed table
max.table.imp <- impute.on.table(max.table.filt)
write.activity.table(max.table.imp, paste0("data/kinact-", strategy, "-imp.tsv"))
