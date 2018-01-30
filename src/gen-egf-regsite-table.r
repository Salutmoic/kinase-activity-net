suppressMessages(library(Biobase))
source("src/optimize-activity-table.r")

load("data/external/esetNR.Rdata")
cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals <- exprs(esetNR)
cptac.start <- 1604
cptac.end <- 1711
cptac.range <- 1604:1711
cptac.conds <- paste(as.character(cptac.range), "290", sep="_")
phospho.vals.sub <- phospho.vals[, -which(colnames(phospho.vals) %in% cptac.conds)]

all.kins <- read.table("data/human-kinome.txt", as.is=TRUE)[,1]

ensp.id.map.full <- read.table("data/ensembl-id-map.tsv", as.is=TRUE)
names(ensp.id.map.full) <- c("ensembl", "gene.name")
ensp.id.map <- subset(ensp.id.map.full, gene.name %in% all.kins)
rownames(ensp.id.map) <- ensp.id.map$ensembl

phospho.vals.sub.ensp <- unlist(lapply(strsplit(rownames(phospho.vals.sub), split="_"),
                                       function(foo) foo[1]))
phospho.vals.sub.gene <- ensp.id.map[phospho.vals.sub.ensp, "gene.name"]

phospho.vals <- phospho.vals.sub[which(phospho.vals.sub.ensp %in% ensp.id.map$ensembl),]
phospho.vals.ensp <- unlist(lapply(strsplit(rownames(phospho.vals), split="_"),
                                   function(foo) foo[1]))
phospho.vals.pos <- unlist(lapply(strsplit(rownames(phospho.vals), split="_"),
                                  function(foo) foo[2]))
phospho.vals.gene <- ensp.id.map[phospho.vals.ensp, "gene.name"]

rownames(phospho.vals) <- paste(phospho.vals.gene, phospho.vals.pos, sep="_")

reg.sites <- read.delim("data/reg_sites.tsv", as.is=TRUE)
reg.sites$MOD_RSD <- sub("^[STY]", "", reg.sites$MOD_RSD)
rownames(reg.sites) <- paste(reg.sites$GENE, reg.sites$MOD_RSD, sep="_")

phosfun <- read.delim("data/phosfun.tsv", as.is=TRUE)
names(phosfun) <- c("prot", "pos", "score")
rownames(phosfun) <- paste(phosfun$prot, phosfun$pos, sep="_")
for (site in rownames(phosfun)){
    if (site %in% rownames(reg.sites))
        phosfun[site, "score"] <- 1.0
}

kin.subs.full <- read.delim("data/reduced_kinase_table.tsv", as.is=TRUE)
kin.subs.full$SUB_MOD_RSD <- sub("^[STY]", "", kin.subs.full$SUB_MOD_RSD)

## egf.kinases <- c("AKT1", "AKT2", "ARAF", "BRAF", "CHUK", "GSK3A", "GSK3B",
##                  "JAK1", "JAK2", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "MTOR",
##                  "PAK1", "PDPK1", "PRKAA2", "PRKACA", "PTK2", "RAF1",
##                  "RPS6KB1", "RPS6KB2", "SGK1", "SRC")
## egf.kinases <- c("AKT1", "ARAF", "BRAF", "CHUK", "GSK3A", "GSK3B",
##                  "JAK1", "JAK2", "MAP2K1", "MAPK1", "MTOR",
##                  "PAK1", "PDPK1", "PRKAA2", "PRKACA", "PTK2", "RPS6KA1"
##                  "RPS6KB1", "RPS6KB2", "SGK1", "SRC")
## erk.kinases <- c("ARAF", "BRAF", "MAP2K1", "MAPK1", "RPS6KB1", "SRC")
erk.kinases <- c("ARAF", "BRAF", "MAP2K1", "MAPK1", "RPS6KB1")

kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                as.is=TRUE)

kinase.cond <- subset(kinase.conditions, Kinase %in% erk.kinases &
                                         Kinase %in% phospho.vals.gene &
                                         Condition %in% colnames(phospho.vals) &
                                         Regulation=="down")

## Remove sites with multiple kinases
kin.subs.sub <- subset(kin.subs.full, GENE %in% kinase.cond$Kinase)
all.subs <- paste(kin.subs.sub$SUB_GENE, kin.subs.sub$SUB_MOD_RSD, sep="_")
dup.subs <- all.subs[which(duplicated(all.subs))]
good.subs <- which(!(all.subs %in% dup.subs))
kin.subs <- kin.subs.sub[good.subs,]

sub.sites.all <- paste(kin.subs$SUB_GENE, kin.subs$SUB_MOD_RSD, sep="_")
sub.sites <- intersect(sub.sites.all, rownames(phospho.vals))
phospho.vals.inhib <- phospho.vals[sub.sites, kinase.cond$Condition]
phospho.vals.inhib <- phospho.vals.inhib[which(apply(phospho.vals.inhib,1,percent.na)<1),
                                         which(apply(phospho.vals.inhib,2,percent.na)<1)]

best.sites <- list()
for (kin in unique(kinase.cond$Kinase)){
    print(kin)
    subs <- subset(kin.subs, GENE==kin)
    sub.sites.all <- paste(subs$SUB_GENE, subs$SUB_MOD_RSD, sep="_")
    sub.sites <- intersect(sub.sites.all, rownames(phospho.vals.inhib))
    if (length(sub.sites) > 1){
        site.na <- apply(phospho.vals.inhib[sub.sites,], 1, percent.na)
        best <- which.min(site.na)
        best.sites[[kin]] <- sub.sites[best]
    }else{
        best.sites[[kin]] <- sub.sites
    }
}
best.sites <- unlist(best.sites)

phospho.vals.opt <- phospho.vals.inhib[best.sites,]

best.conds <- c()
for (n in 1:nrow(phospho.vals.opt)){
    kin <- names(best.sites[n])
    print(kin)
    kin.conds <- subset(kinase.cond, Kinase==kin)$Condition
    print(kin.conds)
    if (length(kin.conds) > 1){
        cond.na <- apply(phospho.vals.opt[,kin.conds], 2, percent.na)
        best <- which.min(cond.na)
        best.conds <- c(best.conds, names(cond.na)[best])
    }else{
        best.conds <- c(best.conds, kin.conds[1])
    }
}

phospho.vals.opt2 <- phospho.vals.opt[,best.conds]
colnames(phospho.vals.opt2) <- sapply(colnames(phospho.vals.opt2), function(x) kinase.cond[which(kinase.cond$Condition==x), "Kinase"])
phospho.vals.opt2
