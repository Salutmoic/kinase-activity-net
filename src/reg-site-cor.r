suppressMessages(library(Biobase))
suppressMessages(library(reshape2))
load("data/external/esetNR.Rdata")

calc.dcg <- function(scores){
    dcg <- sum(sapply(1:length(scores),
                      function(i){
                          scores[i]/log2(i+1)
                      }), na.rm=TRUE)
    return (dcg)
}

cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals.full <- exprs(esetNR)

cptac.start <- 1604
cptac.end <- 1711
cptac.range <- 1604:1711
cptac.conds <- paste(as.character(cptac.range), "290", sep="_")
phospho.vals.sub <- phospho.vals.full[, -which(colnames(phospho.vals.full) %in% cptac.conds)]

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


## reg.sites <- read.delim("data/psiteplus-reg-sites.tsv", as.is=TRUE)
## reg.sites$MOD_RSD <- sub("^[STY]", "", reg.sites$MOD_RSD)
## reg.sites$MOD_RSD <- sub("-p$", "", reg.sites$MOD_RSD)
## rownames(reg.sites) <- paste(reg.sites$GENE, reg.sites$MOD_RSD, sep="_")

phosfun <- read.delim("data/phosfun.tsv", as.is=TRUE)
names(phosfun) <- c("prot", "pos", "score")
rownames(phosfun) <- paste(phosfun$prot, phosfun$pos, sep="_")
## for (site in rownames(phosfun)){
##     if (site %in% rownames(reg.sites))
##         phosfun[site, "score"] <- 1.0
## }
## med.phosfun <- median(phosfun$score)

kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                as.is=TRUE)

kin.pairs <- expand.grid(all.kins, all.kins, stringsAsFactors=FALSE)

kin1s <- c()
kin2s <- c()
p.vals <- c()
cors <- c()
mean.cors <- c()
dcgs <- c()
signed.dcgs <- c()
for (i in 1:nrow(kin.pairs)){
    kin1 <- kin.pairs[i, 1]
    kin2 <- kin.pairs[i, 2]
    if (kin1 == kin2 || !(kin1 %in% phospho.vals.gene) ||
        !(kin2 %in% phospho.vals.gene)){
        kin1s <- c(kin1s, kin1)
        kin2s <- c(kin2s, kin2)
        p.vals <- c(p.vals, NA)
        cors <- c(cors, NA)
        mean.cors <- c(mean.cors, NA)
        dcgs <- c(dcgs, NA)
        signed.dcgs <- c(signed.dcgs, NA)
        next
    }
    kin1.sites <- which(phospho.vals.gene == kin1)
    kin2.sites <- which(phospho.vals.gene == kin2)
    kin1.conds <- kinase.conditions[which(kinase.conditions$Kinase==kin1),
                                    "Condition"]
    kin2.conds <- kinase.conditions[which(kinase.conditions$Kinase==kin2),
                                    "Condition"]
    pair.p.vals <- c()
    pair.cors <- c()
    pair.weights <- c()
    ## best.p.value <- -Inf
    for (kin1.site in kin1.sites){
        for (kin2.site in kin2.sites){
            kin1.phospho <- phospho.vals[kin1.site,]
            kin2.phospho <- phospho.vals[kin2.site,]
            ## Remove conditions where the kinases were inhibited or
            ## activated, because the phosphorylation state of a site
            ## might not reflect the chemically induced activity of
            ## the kinase (e.g. if the site is regulated by a
            ## different kinase, the site might still get
            ## phosphorylated, even if the substrate kinase has been
            ## inhibited.  So we don't calculate correlations on these
            ## conditions.
            kin1.phospho[kin1.conds] <- NA
            kin2.phospho[kin2.conds] <- NA
            kin1.meas.conds <- which(!is.na(kin1.phospho))
            kin2.meas.conds <- which(!is.na(kin2.phospho))
            shared.conds <- intersect(kin1.meas.conds, kin2.meas.conds)
            if (length(shared.conds) < 5)
                next
            ## Get the functional score prediction for both sites
            kin1.phosfun.i <- rownames(phospho.vals)[kin1.site]
            kin2.phosfun.i <- rownames(phospho.vals)[kin2.site]
            kin1.phosfun <- phosfun[kin1.phosfun.i, "score"]
            kin2.phosfun <- phosfun[kin2.phosfun.i, "score"]
            ct <- cor.test(kin1.phospho, kin2.phospho, method="spearman",
                           use="pairwise.complete.obs")
            p.value <- -log10(ct$p.value)
            ## z-transformation of Spearman's rho
            cor.est <- sqrt((length(shared.conds)-3)/1.06)*atanh(ct$estimate)
            if (!is.infinite(p.value)){ ## && p.value > best.p.value){
                pair.p.vals <- c(pair.p.vals, p.value)
                pair.cors <- c(pair.cors, cor.est)
                pair.weights <- c(pair.weights, kin1.phosfun*kin2.phosfun)
            }
        }
    }
    if (length(pair.p.vals) == 0){
        p.val <- NA
        cor.est <- NA
        mean.cor <- NA
        idcg <- NA
        signed.idcg <- NA
    }else{
        p.val <- max(pair.p.vals*pair.weights, na.rm=TRUE)
        j <- which.max(pair.p.vals*pair.weights)
        cor.est <- pair.cors[j]*pair.weights[j]
        mean.cor <- weighted.mean(pair.cors, w=pair.weights, na.rm=TRUE)
        ## DCG
        pair.p.vals <- pair.p.vals/max(pair.p.vals, na.rm=TRUE)
        dcg <- calc.dcg(pair.weights[order(pair.p.vals, decreasing=TRUE)])
        max.dcg <- calc.dcg(pair.weights[order(pair.weights, decreasing=TRUE)])
        min.dcg <- calc.dcg(pair.weights[order(pair.weights, decreasing=FALSE)])
        idcg <- (dcg-min.dcg)/(max.dcg-min.dcg)
        idcg <- idcg * max(pair.p.vals, na.rm=TRUE) * max(pair.weights, na.rm=TRUE)
        signed.weights <- pair.cors * pair.weights
        signed.dcg <- calc.dcg(signed.weights[order(pair.p.vals, decreasing=TRUE)])
        if (is.na(signed.dcg)){
            signed.idcg <- NA
        }else if (signed.dcg < 0){
            best.signed.dcg <- calc.dcg(signed.weights[order(signed.weights)])
            worst.signed.dcg <- calc.dcg(signed.weights[order(signed.weights, decreasing=TRUE)])
            signed.idcg <- (signed.dcg-worst.signed.dcg)/(best.signed.dcg-worst.signed.dcg)
            signed.idcg <- signed.idcg * min(signed.weights, na.rm=TRUE) * max(pair.p.vals, na.rm=TRUE)
        }else{
            best.signed.dcg <- calc.dcg(signed.weights[order(signed.weights, decreasing=TRUE)])
            worst.signed.dcg <- calc.dcg(signed.weights[order(signed.weights)])
            signed.idcg <- (signed.dcg-worst.signed.dcg)/(best.signed.dcg-worst.signed.dcg)
            signed.idcg <- signed.idcg * max(signed.weights, na.rm=TRUE) * max(pair.p.vals, na.rm=TRUE)
        }
    }
    if (is.null(p.val)){
        p.val <- NA
        cor.est <- NA
        mean.cor <- NA
        idcg <- NA
        signed.idcg <- NA
    }
    kin1s <- c(kin1s, kin1)
    kin2s <- c(kin2s, kin2)
    p.vals <- c(p.vals, p.val)
    cors <- c(cors, cor.est)
    mean.cors <- c(mean.cors, mean.cor)
    dcgs <- c(dcgs, idcg)
    signed.dcgs <- c(signed.dcgs, signed.idcg)
}

results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor.p=p.vals)
results$reg.site.cor.p <- results$reg.site.cor.p/max(results$reg.site.cor.p, na.rm=TRUE)
write.table(results[order(results$node1),], "out/reg-site-cor.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

cor.est.results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor.est=cors)
write.table(cor.est.results[order(cor.est.results$node1),],
            "out/reg-site-cor-sign.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

cor.mean.results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor.mean=mean.cors)
write.table(cor.mean.results[order(cor.mean.results$node1),],
            "out/reg-site-mean-cor-sign.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

cor.dcg.results <- data.frame(node1=kin1s, node2=kin2s, reg.site.dcg=dcgs)
write.table(cor.dcg.results[order(cor.dcg.results$node1),],
            "out/reg-site-cor-dcg.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

cor.signed.dcg.results <- data.frame(node1=kin1s, node2=kin2s,
                                     reg.site.signed.dcg=dcgs)
write.table(cor.signed.dcg.results[order(cor.signed.dcg.results$node1),],
            "out/reg-site-cor-dcg-signed.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")


## message("Done with cor, starting cor^2")

## ## kin.cor <- acast(results, node1 ~ node2)
## kin.cor <- matrix(p.vals, nrow=length(all.kins))
## rownames(kin.cor) <- all.kins
## colnames(kin.cor) <- all.kins

## kin1s <- c()
## kin2s <- c()
## p.vals <- c()
## for (i in 1:nrow(kin.pairs)){
##     kin1 <- kin.pairs[i, 1]
##     kin2 <- kin.pairs[i, 2]
##     if (kin1 == kin2){
##         kin1s <- c(kin1s, kin1)
##         kin2s <- c(kin2s, kin2)
##         p.vals <- c(p.vals, NA)
##         next
##     }
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

## results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor2.p=p.vals)
## results$reg.site.cor2.p <- results$reg.site.cor2.p/max(results$reg.site.cor2.p, na.rm=TRUE)

## write.table(results[order(results$node1),], "out/reg-site-cor2.tsv", quote=FALSE,
##             row.names=FALSE, col.names=TRUE, sep="\t")
