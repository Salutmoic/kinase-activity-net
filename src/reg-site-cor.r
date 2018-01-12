suppressMessages(library(Biobase))
suppressMessages(library(reshape2))
load("data/external/esetNR.Rdata")

cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals <- exprs(esetNR)

cptac.start <- 1604
cptac.end <- 1711
cptac.range <- 1604:1711
cptac.conds <- paste(as.character(cptac.range), "290", sep="_")
phospho.vals <- phospho.vals[, -which(colnames(phospho.vals) %in% cptac.conds)]

reg.sites.full <- read.delim("data/reg_sites.tsv", as.is=TRUE)
reg.sites.full$MOD_RSD <- sub("^[STY]", "", reg.sites.full$MOD_RSD)
reg.sites <- reg.sites.full[which(grepl("enzymatic activity, induced",
                                        reg.sites.full$ON_FUNCTION) |
                                  grepl("enzymatic activity, inhibited",
                                        reg.sites.full$ON_FUNCTION)),]

## activity induction => 1, activity inhibition => -1
reg.sites$direction <- as.integer(grepl("induced", reg.sites$ON_FUNCTION))*2 - 1


reg.sites.ensp.full <- merge(phospho.anno, reg.sites,
                             by.x=c("gene_name", "positions"),
                             by.y=c("GENE", "MOD_RSD"))
rownames(reg.sites.ensp.full) <- paste(reg.sites.ensp.full$ensp,
                                       reg.sites.ensp.full$positions,
                                       "P", sep="_")

reg.sites.i <- which(rownames(reg.sites.ensp.full) %in% rownames(phospho.vals))
reg.sites.ensp <- reg.sites.ensp.full[reg.sites.i,]

all.kins <- unique(reg.sites.ensp$gene_name)
kin.pairs <- expand.grid(all.kins, all.kins, stringsAsFactors=FALSE)
kin1s <- c()
kin2s <- c()
p.vals <- c()
cors <- c()
for (i in 1:nrow(kin.pairs)){
    kin1 <- kin.pairs[i, 1]
    kin2 <- kin.pairs[i, 2]
    kin1.reg.sites <- subset(reg.sites.ensp, gene_name==kin1)
    kin2.reg.sites <- subset(reg.sites.ensp, gene_name==kin2)
    pair.p.vals <- c()
    ## best.p.value <- -Inf
    for (kin1.ensp in rownames(kin1.reg.sites)){
        for (kin2.ensp in rownames(kin2.reg.sites)){
            kin1.ensp.dir <- reg.sites.ensp[kin1.ensp, "direction"]
            kin2.ensp.dir <- reg.sites.ensp[kin2.ensp, "direction"]
            kin1.meas.conds <- which(!is.na(phospho.vals[kin1.ensp,]))
            kin2.meas.conds <- which(!is.na(phospho.vals[kin2.ensp,]))
            shared.conds <- intersect(kin1.meas.conds, kin2.meas.conds)
            if (length(shared.conds) < 3)
                next
            kin1.ensp.phospho <- phospho.vals[kin1.ensp,] * kin1.ensp.dir
            kin2.ensp.phospho <- phospho.vals[kin2.ensp,] * kin2.ensp.dir
            ct <- cor.test(kin1.ensp.phospho, kin2.ensp.phospho, method="spearman",
                           use="pairwise.complete.obs")
            p.value <- -log10(ct$p.value)
            if (!is.infinite(p.value)){ ## && p.value > best.p.value){
                ## best.p.value <- p.value
                ## best.cor.est <- ct$estimate
                pair.p.vals <- c(pair.p.vals, p.value)
            }
        }
    }
    ## if (is.infinite(best.p.value)){
    ##     best.p.value <- NA
    ##     best.cor.est <- NA
    ## }
    if (length(pair.p.vals) == 0){
        p.val <- NA
    }else{
        p.val <- max(pair.p.vals, na.rm=TRUE)
    }
    if (is.null(p.val))
        p.val <- NA
    kin1s <- c(kin1s, kin1)
    kin2s <- c(kin2s, kin2)
    ## p.vals <- c(p.vals, best.p.value)
    p.vals <- c(p.vals, p.val)
    ## cors <- c(cors, best.cor.est)
}

results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor.p=p.vals)
results$reg.site.cor.p <- results$reg.site.cor.p/max(results$reg.site.cor.p, na.rm=TRUE)

write.table(results[order(results$node1),], "out/reg-site-cor.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

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
