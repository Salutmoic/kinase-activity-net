suppressMessages(library(Biobase))
suppressMessages(library(reshape2))
suppressMessages(library(preprocessCore))

load("data/external/esetNR.Rdata")
load("data/ks-kinact.Rdata")
kin.act.full <- log.wKSEA.kinase_condition.clean
cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals.full <- exprs(esetNR)
kin.sub.tbl <- read.delim("data/psiteplus-kinase-substrates.tsv", as.is=TRUE,
                          sep="\t")
## kin.sub.tbl$SUB_MOD_RSD <- sub("[STY]", "", kin.sub.tbl$SUB_MOD_RSD)
kin.overlap <- read.delim("data/psiteplus-kinase-substrate-overlap.tsv",
                          as.is=TRUE)

wilkes.conds <- rownames(subset(cond.anno, experiment_id == "272"))
phospho.vals.sub <- phospho.vals.full[, which(colnames(phospho.vals.full) %in% wilkes.conds)]
kin.act <- kin.act.full[, which(colnames(kin.act.full) %in% wilkes.conds)]
