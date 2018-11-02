suppressMessages(library(Biobase))
suppressMessages(library(reshape2))
suppressMessages(library(preprocessCore))
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
phospho.vals.full <- exprs(esetNR)
kin.sub.tbl <- read.delim("data/psiteplus-kinase-substrates.tsv", as.is=TRUE,
                          sep="\t")
## kin.sub.tbl$SUB_MOD_RSD <- sub("[STY]", "", kin.sub.tbl$SUB_MOD_RSD)
kin.overlap <- read.delim("data/psiteplus-kinase-substrate-overlap.tsv",
                          as.is=TRUE)

## cptac.290.range <- 1604:1711
## cptac.292.range <- 1781:1860
## cptac.290.conds <- paste(as.character(cptac.290.range), "290", sep="_")
## cptac.292.conds <- paste(as.character(cptac.292.range), "292", sep="_")
## cptac.conds <- c(cptac.290.conds, cptac.292.conds)
## cptac.pubs <- c("176", "177")
## cptac.conds <- rownames(subset(cond.anno, publication %in% cptac.pubs & experiment_id != "291"))
cptac.conds <- rownames(subset(cond.anno, experiment_id == "293"))
phospho.vals.sub <- phospho.vals.full[, which(colnames(phospho.vals.full) %in% cptac.conds)]

all.kins <- read.table("data/human-kinome.txt", as.is=TRUE)[,1]

psites <- read.table("data/pride-phosphosites.tsv", as.is=TRUE, sep="\t")
rownames(psites) <- paste(psites[,1], psites[,2], sep="_")

ensp.id.map.full <- read.table("data/ensembl-id-map.tsv", as.is=TRUE)
names(ensp.id.map.full) <- c("ensembl", "gene.name")
ensp.id.map <- subset(ensp.id.map.full, gene.name %in% all.kins)
rownames(ensp.id.map) <- ensp.id.map$ensembl

phospho.vals.ensp <- unlist(lapply(strsplit(rownames(phospho.vals.sub), split="_"),
                                   function(foo) foo[1]))
phospho.vals <- phospho.vals.sub[which(phospho.vals.ensp %in% ensp.id.map$ensembl),]
phospho.vals.ensp <- unlist(lapply(strsplit(rownames(phospho.vals), split="_"),
                                   function(foo) foo[1]))
phospho.vals.pos <- unlist(lapply(strsplit(rownames(phospho.vals), split="_"),
                                  function(foo) foo[2]))
phospho.vals.gene <- ensp.id.map[phospho.vals.ensp, "gene.name"]

rownames(phospho.vals) <- paste(phospho.vals.gene, phospho.vals.pos, sep="_")

good.sites <- which(rownames(phospho.vals) %in% rownames(psites))
phospho.vals <- phospho.vals[good.sites,]
phospho.vals.gene <- phospho.vals.gene[good.sites]

phospho.vals.norm <- normalize.quantiles(phospho.vals)
rownames(phospho.vals.norm) <- rownames(phospho.vals)
colnames(phospho.vals.norm) <- colnames(phospho.vals)
phospho.vals <- phospho.vals.norm

## reg.sites <- read.delim("data/psiteplus-reg-sites.tsv", as.is=TRUE)
## reg.sites$MOD_RSD <- sub("^[STY]", "", reg.sites$MOD_RSD)
## reg.sites$MOD_RSD <- sub("-p$", "", reg.sites$MOD_RSD)
## rownames(reg.sites) <- paste(reg.sites$GENE, reg.sites$MOD_RSD, sep="_")

phosfun <- read.table("data/phosfun.tsv", as.is=TRUE)
names(phosfun) <- c("prot", "pos", "score")
rownames(phosfun) <- paste(phosfun$prot, phosfun$pos, sep="_")
## for (site in rownames(phosfun)){
##     if (site %in% rownames(reg.sites))
##         phosfun[site, "score"] <- 1.0
## }
## med.phosfun <- median(phosfun$score)

kin.overlap.sub <- kin.overlap
good.kins <- unique(c(kin.overlap.sub$kinase1, kin.overlap.sub$kinase2))

## Remove redundant kinases
redundant.kins <- get.redundant.kins(kin.overlap.sub, kin.act)

kin.act <- kin.act[setdiff(good.kins, redundant.kins),]

## Filter out kinase activity predictions that were not based on
## enough sites
kin.act <- filter.act.preds(kin.act, phospho.anno, phospho.vals, kin.sub.tbl,
                            min.sites)

## We look for "empty" rows here...not really empty, but anything
## that's practically empty (>90% missing values).
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

## Remove experimental conditions
kin.act.filt <- kin.act.filt[,setdiff(colnames(kin.act.filt), kinase.conditions$Condition)]
kin.act.filt <- kin.act.filt[,intersect(colnames(phospho.vals), colnames(kin.act.filt))]
phospho.vals <- phospho.vals[,colnames(kin.act.filt)]

p.vals <- c()
cors <- c()
for (kin in all.kins){
    if (!(kin %in% rownames(kin.act.filt)) || !(kin %in% phospho.vals.gene)){
        p.vals <- c(p.vals, NA)
        cors <- c(cors, NA)
        next
    }
    kin.sites <- which(phospho.vals.gene == kin)
    kin.sites.phosfun <- phosfun[kin.sites,]
    kin.acts <- kin.act.filt[kin,]
    kin.acts.conds <- which(!is.na(kin.acts))
    site.p.vals <- c()
    site.cors <- c()
    for (kin.site in kin.sites){
        site.phospho <- phospho.vals[kin.site,]
        site.phospho.conds <- which(!is.na(site.phospho))
        shared.conds <- intersect(site.phospho.conds, kin.acts.conds)
        if (length(shared.conds) < 5){
            site.p.vals <- c(site.p.vals, NA)
            site.cors <- c(site.cors, NA)
            next
        }
        ## Get the functional score prediction for both sites
        kin.phosfun.i <- rownames(phospho.vals)[kin.site]
        if (!(kin.phosfun.i %in% rownames(kin.sites.phosfun) ||
              is.na(kin.sites.phosfun[kin.phosfun.i,]))){
            kin.phosfun <- med.phosfun
        }else{
            kin.phosfun <- kin.sites.phosfun[kin.phosfun.i, "score"]
        }
        ct <- cor.test(site.phospho, kin.acts, method="spearman",
                       use="pairwise.complete.obs")
        p.value <- -log10(ct$p.value)
        cor.est <- sqrt((length(shared.conds)-3)/1.06)*atanh(ct$estimate)
        if (!is.infinite(p.value)){ ## && p.value > best.p.value){
            site.p.vals <- c(site.p.vals, p.value)
            ## z-transformation of Spearman's rho
            site.cors <- c(site.cors, cor.est)
        }else{
            site.p.vals <- c(site.p.vals, NA)
            site.cors <- c(site.cors, NA)
        }
    }
    if (length(site.p.vals) == 0){
        p.val <- NA
        cor.est <- NA
    }else{
        p.val <- weighted.mean(site.p.vals, w=kin.sites.phosfun$score, na.rm=TRUE)
        cor.est <- weighted.mean(site.cors, w=kin.sites.phosfun$score, na.rm=TRUE)
    }
    if (is.null(p.val)){
        p.val <- NA
        cor.est <- NA
    }
    p.vals <- c(p.vals, p.val)
    cors <- c(cors, cor.est)
}

sign.guess.results <- data.frame(kinase=all.kins, sign.guess=cors)
write.table(sign.guess.results[order(sign.guess.results$kinase),],
            "out/kinase-reg-sign-guess.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")
