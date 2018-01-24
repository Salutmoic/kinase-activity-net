suppressMessages(library(Biobase))
suppressMessages(library(reshape2))
load("data/external/esetNR.Rdata")

cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals.full <- exprs(esetNR)

cptac.start <- 1604
cptac.end <- 1711
cptac.range <- 1604:1711
cptac.conds <- paste(as.character(cptac.range), "290", sep="_")
phospho.vals.sub <- phospho.vals.full[, -which(colnames(phospho.vals.full) %in% cptac.conds)]

all.kins <- read.table("data/human-kinome.txt", as.is=TRUE)[,1]

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

reg.sites <- read.delim("data/reg_sites.tsv", as.is=TRUE)
reg.sites$MOD_RSD <- sub("^[STY]", "", reg.sites$MOD_RSD)
## reg.sites <- reg.sites.full[which(grepl("enzymatic activity, induced",
##                                         reg.sites.full$ON_FUNCTION) |
##                                   grepl("enzymatic activity, inhibited",
##                                         reg.sites.full$ON_FUNCTION)),]
rownames(reg.sites) <- paste(reg.sites$GENE, reg.sites$MOD_RSD, sep="_")

## activity induction => 1, activity inhibition => -1
## reg.sites$direction <- as.integer(grepl("induced", reg.sites$ON_FUNCTION))*2 - 1


phosfun <- read.delim("data/phosfun.tsv", as.is=TRUE)
names(phosfun) <- c("prot", "pos", "score")
rownames(phosfun) <- paste(phosfun$prot, phosfun$pos, sep="_")
for (site in rownames(phosfun)){
    if (site %in% rownames(reg.sites))
        phosfun[site, "score"] <- 1.0
}
med.phosfun <- median(phosfun$score)

kin.pairs <- expand.grid(all.kins, all.kins, stringsAsFactors=FALSE)
kin1s <- c()
kin2s <- c()
p.vals <- c()
cors <- c()
for (i in 1:nrow(kin.pairs)){
    kin1 <- kin.pairs[i, 1]
    kin2 <- kin.pairs[i, 2]
    if (kin1 == kin2 || !(kin1 %in% phospho.vals.gene) ||
        !(kin2 %in% phospho.vals.gene)){
        kin1s <- c(kin1s, kin1)
        kin2s <- c(kin2s, kin2)
        p.vals <- c(p.vals, NA)
        next
    }
    kin1.sites <- which(phospho.vals.gene == kin1)
    kin2.sites <- which(phospho.vals.gene == kin2)
    pair.p.vals <- c()
    ## best.p.value <- -Inf
    for (kin1.site in kin1.sites){
        for (kin2.site in kin2.sites){
            kin1.meas.conds <- which(!is.na(phospho.vals[kin1.site,]))
            kin2.meas.conds <- which(!is.na(phospho.vals[kin2.site,]))
            shared.conds <- intersect(kin1.meas.conds, kin2.meas.conds)
            if (length(shared.conds) < 5)
                next
            kin1.phosfun.i <- rownames(phospho.vals)[kin1.site]
            kin2.phosfun.i <- rownames(phospho.vals)[kin2.site]
            if (!(kin1.phosfun.i %in% rownames(phosfun))){
                kin1.phosfun <- med.phosfun
            }else{
                kin1.phosfun <- phosfun[kin1.phosfun.i, "score"]
            }
            if (!(kin2.phosfun.i %in% rownames(phosfun))){
                kin2.phosfun <- med.phosfun
            }else{
                kin2.phosfun <- phosfun[kin2.phosfun.i, "score"]
            }
            kin1.phospho <- phospho.vals[kin1.site,]
            kin2.phospho <- phospho.vals[kin2.site,]
            ct <- cor.test(kin1.phospho, kin2.phospho, method="spearman",
                           use="pairwise.complete.obs")
            ## p.value <- -log10(ct$p.value/(kin1.phosfun*kin2.phosfun))
            p.value <- -log10(ct$p.value)*kin1.phosfun*kin2.phosfun
            if (!is.infinite(p.value)){ ## && p.value > best.p.value){
                ## best.p.value <- p.value
                ## best.cor.est <- ct$estimate
                pair.p.vals <- c(pair.p.vals, p.value)
            }
        }
    }
    if (length(pair.p.vals) == 0){
        p.val <- NA
    }else{
        p.val <- max(pair.p.vals, na.rm=TRUE)
    }
    if (is.null(p.val))
        p.val <- NA
    kin1s <- c(kin1s, kin1)
    kin2s <- c(kin2s, kin2)
    p.vals <- c(p.vals, p.val)
}

results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor.p=p.vals)
results$reg.site.cor.p <- results$reg.site.cor.p/max(results$reg.site.cor.p, na.rm=TRUE)

write.table(results[order(results$node1),], "out/reg-site-cor.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

message("Done with cor, starting cor^2")

## kin.cor <- acast(results, node1 ~ node2)
kin.cor <- matrix(p.vals, nrow=length(all.kins))
rownames(kin.cor) <- all.kins
colnames(kin.cor) <- all.kins

kin1s <- c()
kin2s <- c()
p.vals <- c()
for (i in 1:nrow(kin.pairs)){
    kin1 <- kin.pairs[i, 1]
    kin2 <- kin.pairs[i, 2]
    if (kin1 == kin2){
        kin1s <- c(kin1s, kin1)
        kin2s <- c(kin2s, kin2)
        p.vals <- c(p.vals, NA)
        next
    }
    kin1.cor.pairs <- which(!is.na(kin.cor[kin1,]))
    kin2.cor.pairs <- which(!is.na(kin.cor[kin2,]))
    shared.pairs <- intersect(kin1.cor.pairs, kin2.cor.pairs)
    if (length(shared.pairs) < 3){
        p.value <- NA
    }else{
        kin1.cor <- kin.cor[kin1,]
        kin2.cor <- kin.cor[kin2,]
        ct <- cor.test(kin1.cor, kin2.cor, method="spearman",
                       use="pairwise.complete.obs")
        p.value <- -log10(ct$p.value)
        if (is.infinite(p.value))
            p.value <- NA
    }
    kin1s <- c(kin1s, kin1)
    kin2s <- c(kin2s, kin2)
    p.vals <- c(p.vals, p.value)
}

results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor2.p=p.vals)
results$reg.site.cor2.p <- results$reg.site.cor2.p/max(results$reg.site.cor2.p, na.rm=TRUE)

write.table(results[order(results$node1),], "out/reg-site-cor2.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")
