suppressMessages(library(Biobase))
suppressMessages(library(reshape2))
suppressMessages(library(preprocessCore))
load("data/external/esetNR.Rdata")

calc.dcg <- function(scores){
    dcg <- sum(sapply(1:length(scores),
                      function(i){
                          scores[i]/log2(i+1)
                      }), na.rm=TRUE)
    return (dcg)
}

normalize.p.vals <- function(cor.p){
    min.p <- min(cor.p[which(cor.p>0)], na.rm=TRUE)
    cor.p[which(cor.p==0)] <- min.p
    cor.p <- log(cor.p)
    min.p <- min(cor.p, na.rm=TRUE)
    max.p <- max(cor.p, na.rm=TRUE)
    cor.p <- (cor.p-min.p)/(max.p-min.p)
    return (cor.p)
}

min.conds <- 5

cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals.full <- exprs(esetNR)

cptac.start <- 1604
cptac.end <- 1711
cptac.range <- 1604:1711
cptac.conds <- paste(as.character(cptac.range), "290", sep="_")
phospho.vals.sub <- phospho.vals.full[, -which(colnames(phospho.vals.full) %in% cptac.conds)]

exp.cond.table <- table(cond.anno$experiment_id)
## good.exps <- sort(names(exp.cond.table)[which(exp.cond.table>=min.conds)])
good.exps <- c("272", "293")

argv <- commandArgs(TRUE)
if (length(argv) != 1){
    stop(paste0("USAGE: <script> EXP_NUM (EXP_NUM=[1-", length(good.exps), "])"))
}
exp.num <- as.integer(argv[1])
if (exp.num <= 0 || exp.num > length(good.exps))
    stop(paste0("USAGE: <script> EXP_NUM (EXP_NUM=[1-", length(good.exps), "])"))

experiment.id <- good.exps[exp.num]
exp.conds <- subset(cond.anno, experiment_id==experiment.id)
good.conds <- intersect(rownames(exp.conds), colnames(phospho.vals.sub))
phospho.vals.exp <- phospho.vals.sub[, good.conds]

all.kins <- read.table("data/human-kinome.txt", as.is=TRUE)[,1]

psites <- read.table("data/pride-phosphosites.tsv", as.is=TRUE, sep="\t")
rownames(psites) <- paste(psites[,1], psites[,2], sep="_")

ensp.id.map.full <- read.table("data/ensembl-id-map.tsv", as.is=TRUE)
names(ensp.id.map.full) <- c("ensembl", "gene.name")
ensp.id.map <- subset(ensp.id.map.full, gene.name %in% all.kins)
rownames(ensp.id.map) <- ensp.id.map$ensembl

phospho.vals.ensp <- unlist(lapply(strsplit(rownames(phospho.vals.exp), split="_"),
                                   function(foo) foo[1]))
phospho.vals <- phospho.vals.exp[which(phospho.vals.ensp %in% ensp.id.map$ensembl),]
phospho.vals.ensp <- unlist(lapply(strsplit(rownames(phospho.vals), split="_"),
                                   function(foo) foo[1]))
phospho.vals.pos <- unlist(lapply(strsplit(rownames(phospho.vals), split="_"),
                                  function(foo) foo[2]))
phospho.vals.gene <- ensp.id.map[phospho.vals.ensp, "gene.name"]
quant.prots <- unique(phospho.vals.gene)

rownames(phospho.vals) <- paste(phospho.vals.gene, phospho.vals.pos, sep="_")

good.sites <- which(rownames(phospho.vals) %in% rownames(psites))
phospho.vals <- phospho.vals[good.sites,]
phospho.vals.gene <- phospho.vals.gene[good.sites]
phospho.vals.pos <- phospho.vals.pos[good.sites]

phospho.vals.norm <- normalize.quantiles(phospho.vals)
rownames(phospho.vals.norm) <- rownames(phospho.vals)
colnames(phospho.vals.norm) <- colnames(phospho.vals)
phospho.vals <- phospho.vals.norm

phosfun <- read.table("data/phosfun.tsv", as.is=TRUE)
names(phosfun) <- c("prot", "pos", "score")
rownames(phosfun) <- paste(phosfun$prot, phosfun$pos, sep="_")

site.signs <- read.delim("out/reg-site-bart-sign-preds-rescaled.tsv", as.is=TRUE)
names(site.signs) <- c("prot", "pos", "score")
rownames(site.signs) <- paste(site.signs$prot, site.signs$pos, sep="_")

kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                as.is=TRUE)
kinase.conditions <- subset(kinase.conditions,
                            Condition %in% colnames(phospho.vals))

kin.pairs <- expand.grid(all.kins, all.kins, stringsAsFactors=FALSE)

num.quant <- function(v) length(which(!is.na(v)))

kin1s <- c()
kin2s <- c()
results <- NULL
p.vals <- c()
all.p.vals <- c()
all.weights1 <- c()
all.weights2 <- c()
all.signs1 <- c()
all.signs2 <- c()
cors <- c()
all.cors <- c()
all.pos1 <- c()
all.pos2 <- c()
mean.cors <- c()
for (i in 1:nrow(kin.pairs)){
    kin1 <- kin.pairs[i, 1]
    kin2 <- kin.pairs[i, 2]
    if (kin1 == kin2 || !(kin1 %in% phospho.vals.gene) ||
        !(kin2 %in% phospho.vals.gene)){
        kin1s <- c(kin1s, kin1)
        kin2s <- c(kin2s, kin2)
        p.vals <- c(p.vals, NA)
        all.p.vals <- c(all.p.vals, NA)
        all.weights1 <- c(all.weights1, NA)
        all.weights2 <- c(all.weights2, NA)
        all.signs1 <- c(all.signs1, NA)
        all.signs2 <- c(all.signs2, NA)
        cors <- c(cors, NA)
        all.cors <- c(all.cors, NA)
        mean.cors <- c(mean.cors, NA)
        all.pos1 <- c(all.pos1, NA)
        all.pos2 <- c(all.pos2, NA)
        next
    }
    kin1.all.sites <- which(phospho.vals.gene == kin1)
    kin1.good.sites <- sapply(kin1.all.sites,
                              function(site){
                                  site.quant <- num.quant(phospho.vals[site,])
                                  return(site.quant>=min.conds)
                              })
    kin1.sites <- kin1.all.sites[which(kin1.good.sites)]
    kin2.all.sites <- which(phospho.vals.gene == kin2)
    kin2.good.sites <- sapply(kin2.all.sites,
                              function(site){
                                  site.quant <- num.quant(phospho.vals[site,])
                                  return(site.quant>=min.conds)
                              })
    if (!any(kin1.good.sites) || !any(kin2.good.sites)){
        kin1s <- c(kin1s, kin1)
        kin2s <- c(kin2s, kin2)
        p.vals <- c(p.vals, NA)
        all.p.vals <- c(all.p.vals, NA)
        all.weights1 <- c(all.weights1, NA)
        all.weights2 <- c(all.weights2, NA)
        all.signs1 <- c(all.signs1, NA)
        all.signs2 <- c(all.signs2, NA)
        cors <- c(cors, NA)
        all.cors <- c(all.cors, NA)
        mean.cors <- c(mean.cors, NA)
        all.pos1 <- c(all.pos1, NA)
        all.pos2 <- c(all.pos2, NA)
        next
    }
    kin2.sites <- kin2.all.sites[which(kin2.good.sites)]

    kin1.conds <- kinase.conditions[which(kinase.conditions$Kinase==kin1),
                                    "Condition"]
    kin1.conds <- intersect(kin1.conds, colnames(phospho.vals))
    kin2.conds <- kinase.conditions[which(kinase.conditions$Kinase==kin2),
                                    "Condition"]
    kin2.conds <- intersect(kin2.conds, colnames(phospho.vals))
    pair.p.vals <- c()
    pair.cors <- c()
    pair.weights <- c()
    pair.signs <- c()
    weights1 <- c()
    weights2 <- c()
    signs1 <- c()
    signs2 <- c()
    sites1 <- c()
    sites2 <- c()
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
            ## Get the functional score prediction for both sites
            kin1.phosfun.i <- rownames(phospho.vals)[kin1.site]
            kin2.phosfun.i <- rownames(phospho.vals)[kin2.site]
            kin1.phosfun <- phosfun[kin1.phosfun.i, "score"]
            kin2.phosfun <- phosfun[kin2.phosfun.i, "score"]
            kin1.sign <- site.signs[kin1.phosfun.i, "score"]
            kin2.sign <- site.signs[kin2.phosfun.i, "score"]
            kin1.phospho[kin1.conds] <- NA
            kin2.phospho[kin2.conds] <- NA
            if (all(is.na(kin1.phospho)) || all(is.na(kin2.phospho))){
                ## pair.p.vals <- c(pair.p.vals, NA)
                ## pair.cors <- c(pair.cors, NA)
                ## pair.weights <- c(pair.weights, kin1.phosfun*kin2.phosfun)
                ## weights1 <- c(weights1, kin1.phosfun)
                ## weights2 <- c(weights2, kin2.phosfun)
                ## signs1 <- c(signs1, kin1.sign)
                ## signs2 <- c(signs2, kin2.sign)
                ## sites1 <- c(sites1, phospho.vals.pos[kin1.site])
                ## sites2 <- c(sites2, phospho.vals.pos[kin2.site])
                next
            }
            kin1.meas.conds <- which(!is.na(kin1.phospho))
            kin2.meas.conds <- which(!is.na(kin2.phospho))
            shared.conds <- intersect(kin1.meas.conds, kin2.meas.conds)
            if (length(shared.conds) < min.conds){
                pair.p.vals <- c(pair.p.vals, NA)
                pair.cors <- c(pair.cors, NA)
                pair.weights <- c(pair.weights, kin1.phosfun*kin2.phosfun)
                weights1 <- c(weights1, kin1.phosfun)
                weights2 <- c(weights2, kin2.phosfun)
                signs1 <- c(signs1, kin1.sign)
                signs2 <- c(signs2, kin2.sign)
                sites1 <- c(sites1, phospho.vals.pos[kin1.site])
                sites2 <- c(sites2, phospho.vals.pos[kin2.site])
                pair.signs <- c(pair.signs, sign(kin1.sign)*sign(kin2.sign))
                next
            }
            ct <- cor.test(kin1.phospho, kin2.phospho, method="spearman",
                           use="pairwise.complete.obs")
            if (ct$p.value==0){
                p.value <- 6
            }else{
                p.value <- -log10(ct$p.value)
            }
            ## z-transformation of Spearman's rho
            if (ct$estimate == 1){
                estimate <- 1 - 1e-6
            }else if (ct$estimate == -1){
                estimate <- -1 + 1e-6
            }else{
                estimate <- ct$estimate
            }
            cor.est <- sqrt((length(shared.conds)-3)/1.06)*atanh(estimate)
            if (sign(cor.est) != sign(kin1.sign))
                cor.est <- NA
            pair.p.vals <- c(pair.p.vals, p.value)
            pair.cors <- c(pair.cors, cor.est)
            pair.weights <- c(pair.weights, kin1.phosfun*kin2.phosfun)
            weights1 <- c(weights1, kin1.phosfun)
            weights2 <- c(weights2, kin2.phosfun)
            signs1 <- c(signs1, kin1.sign)
            signs2 <- c(signs2, kin2.sign)
            sites1 <- c(sites1, phospho.vals.pos[kin1.site])
            sites2 <- c(sites2, phospho.vals.pos[kin2.site])
            pair.signs <- c(pair.signs, sign(kin1.sign)*sign(kin2.sign))
        }
    }
    if (length(pair.p.vals) == 0){
        p.val <- NA
        cor.est <- NA
        mean.cor <- NA
        cor.list <- NA
        p.val.list <- NA
        weights1.list <- NA
        weights2.list <- NA
        signs1.list <- NA
        signs2.list <- NA
        pos1.list <- NA
        pos2.list <- NA
    }else{
        p.val <- max(pair.p.vals*pair.weights, na.rm=TRUE)
        cor.ests <- pair.cors * pair.weights * pair.signs
        min.cor.est <- min(cor.ests, na.rm=TRUE)
        max.cor.est <- max(cor.ests, na.rm=TRUE)
        if (abs(min.cor.est) > max.cor.est){
            cor.est <- min.cor.est
        }else{
            cor.est <- max.cor.est
        }
        ## j <- which.max(pair.p.vals*pair.weights)
        ## cor.est <- pair.cors[j]*pair.weights[j]*pair.signs[j]
        mean.cor <- weighted.mean(pair.cors, w=pair.weights, na.rm=TRUE)
        p.val.list <- paste(pair.p.vals, collapse=",")
        weights1.list <- paste(weights1, collapse=",")
        weights2.list <- paste(weights2, collapse=",")
        signs1.list <- paste(signs1, collapse=",")
        signs2.list <- paste(signs2, collapse=",")
        cor.list <- paste(pair.cors, collapse=",")
        pos1.list <- paste(sites1, collapse=",")
        pos2.list <- paste(sites2, collapse=",")
    }
    if (is.null(p.val)){
        p.val <- NA
        cor.est <- NA
        mean.cor <- NA
        p.val.list <- NA
        weights1.list <- NA
        weights2.list <- NA
        signs1.list <- NA
        signs2.list <- NA
        cor.list <- NA
        pos1.list <- NA
        pos2.list <- NA
    }
    kin1s <- c(kin1s, kin1)
    kin2s <- c(kin2s, kin2)
    p.vals <- c(p.vals, p.val)
    all.p.vals <- c(all.p.vals, p.val.list)
    all.weights1 <- c(all.weights1, weights1.list)
    all.weights2 <- c(all.weights2, weights2.list)
    all.signs1 <- c(all.signs1, signs1.list)
    all.signs2 <- c(all.signs2, signs2.list)
    cors <- c(cors, cor.est)
    all.cors <- c(all.cors, cor.list)
    mean.cors <- c(mean.cors, mean.cor)
    all.pos1 <- c(all.pos1, pos1.list)
    all.pos2 <- c(all.pos2, pos2.list)
}
results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor.p=p.vals,
                      all.p.vals=all.p.vals, all.weights1=all.weights1,
                      all.weights2=all.weights2, kin1.pos=all.pos1,
                      kin2.pos=all.pos2)
results$reg.site.cor.p.norm <- normalize.p.vals(results$reg.site.cor.p)
if (all(is.na(results$reg.site.cor.p.norm))){
    stop("No results")
}
colnames(results)[ncol(results)] <- paste0("cor.p.", experiment.id)
write.table(results[order(results$node1),],
            paste0("out/reg-site-cor-by-exp/reg-site-cor-exp-", experiment.id, "-annot.tsv"),
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

cor.est.results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor.est=cors,
                              all.cor.ests=all.cors, all.weights1=all.weights1,
                              all.weights2=all.weights2, all.signs1=all.signs1,
                              all.signs2=all.signs2, kin1.pos=all.pos1,
                              kin2.pos=all.pos2)
colnames(cor.est.results)[3] <- paste0("cor.est.", experiment.id)
write.table(cor.est.results[order(cor.est.results$node1),],
            paste0("out/reg-site-cor-by-exp/reg-site-cor-sign-exp-",
                   experiment.id, "-annot.tsv"),
            quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

## cor.mean.results <- data.frame(node1=kin1s, node2=kin2s, reg.site.cor.mean=mean.cors)
## colnames(cor.mean.results)[ncol(cor.mean.results)] <- paste0("cor.mean.", experiment.id)
## write.table(cor.mean.results[order(cor.mean.results$node1),],
##             paste0("out/reg-site-cor-by-exp/reg-site-mean-cor-sign-exp-",
##                    experiment.id, ".tsv"),
##             quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
