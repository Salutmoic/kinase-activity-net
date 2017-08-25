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

cptac.pubs <- c("176", "177")
## cptac.exps <-  c("291", "292")
cptac.conds <- rownames(subset(cond.anno, publication %in% cptac.pubs & experiment_id != "291"))
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
good.sites.i <- which(rownames(phospho.vals) %in% rownames(psites))
good.sites <- sort(rownames(phospho.vals)[good.sites.i])
phospho.vals <- phospho.vals[good.sites,]
phospho.vals.gene <- sort(phospho.vals.gene[good.sites.i])

phosfun <- read.table("data/phosfun.tsv", as.is=TRUE)
names(phosfun) <- c("prot", "pos", "score")
rownames(phosfun) <- paste(phosfun$prot, phosfun$pos, sep="_")
phosfun <- phosfun[good.sites,]

site.signs <- read.delim("out/reg-site-bart-sign-preds.tsv", as.is=TRUE)
names(site.signs) <- c("prot", "pos", "score")
rownames(site.signs) <- paste(site.signs$prot, site.signs$pos, sep="_")
site.signs <- site.signs[good.sites,]

kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                as.is=TRUE)
## inhib.conds <- subset(kinase.conditions, Regulation=="down")$Condition
disrupt.conds <- kinase.conditions$Condition
good.conds <- setdiff(colnames(phospho.vals), disrupt.conds)
phospho.vals <- phospho.vals[,good.conds]


site.signs$score[site.signs$score<0] <- -site.signs$score[site.signs$score<0]/min(site.signs$score, na.rm=TRUE)
site.signs$score[site.signs$score>0] <- site.signs$score[site.signs$score>0]/max(site.signs$score, na.rm=TRUE)

weights <- phosfun$score * site.signs$score

weighted.phospho.vals <- apply(phospho.vals, 2, `*`, weights)

## phospho.vals.adj <- apply(phospho.vals, 2,
##                           function(col){
##                               return (col-median(col, na.rm=TRUE))
##                           })

all.kins <- unique(phospho.vals.gene)
reg.site.kin.act <- t(sapply(all.kins,
                    function(kin){
                        i <- which(phospho.vals.gene==kin)
                       if (length(i)==1)
                           return (weighted.phospho.vals[i,])
                        kin.sites <- weighted.phospho.vals[i,]
                        kin.phosfun <- phosfun[i,"score"]
                        acts <- apply(kin.sites, 2,
                                      function(col){
                                          if (all(is.na(col)))
                                              return (NA)
                                          return(sum(col, na.rm=TRUE)/length(which(!is.na(col))))
                                          ## col <- col[order(kin.phosfun, decreasing=TRUE)]
                                          ## dcg <- calc.dcg(col)
                                          ## if (dcg < 0){
                                          ##     col <- sort(col, decreasing=FALSE, na.last=TRUE)
                                          ##     best.dcg <- calc.dcg(col)
                                          ##     col <- sort(col, decreasing=TRUE, na.last=FALSE)
                                          ##     worst.dcg <- calc.dcg(col)
                                          ##     sign <- -1
                                          ## }else{
                                          ##     col <- sort(col, decreasing=TRUE, na.last=TRUE)
                                          ##     best.dcg <- calc.dcg(col)
                                          ##     col <- sort(col, decreasing=FALSE, na.last=FALSE)
                                          ##     worst.dcg <- calc.dcg(col)
                                          ##     sign <- 1
                                          ## }
                                          ## if (best.dcg == worst.dcg)
                                          ##     return (sign(dcg))
                                          ## return (sign*(dcg-worst.dcg)/(best.dcg-worst.dcg))
                                      })
                        return (unlist(acts))
                    }))

kin.pairs <- expand.grid(all.kins, all.kins, stringsAsFactors=FALSE)

kin1s <- c()
kin2s <- c()
p.vals <- c()
cors <- c()
for (i in 1:nrow(kin.pairs)){
    kin1 <- kin.pairs[i, 1]
    kin2 <- kin.pairs[i, 2]
    if (kin1 == kin2){
        kin1s <- c(kin1s, kin1)
        kin2s <- c(kin2s, kin2)
        p.vals <- c(p.vals, NA)
        cors <- c(cors, NA)
        next
    }
    kin1.act <- kin.act[kin1,]
    kin2.act <- kin.act[kin2,]
    kin1.meas.conds <- which(!is.na(kin1.act))
    kin2.meas.conds <- which(!is.na(kin2.act))
    shared.conds <- intersect(kin1.meas.conds, kin2.meas.conds)
    if (length(shared.conds) < 5){
        p.value <- NA
        cor.est <- NA
    }else{
        ct <- cor.test(kin1.act, kin2.act, method="spearman",
                       use="pairwise.complete.obs")
        p.value <- -log10(ct$p.value)
        ## z-transformation of Spearman's rho
        cor.est <- sqrt((length(shared.conds)-3)/1.06)*atanh(ct$estimate)
        if (is.infinite(p.value)){
            p.value <- NA
            cor.est <- NA
        }
    }
    kin1s <- c(kin1s, kin1)
    kin2s <- c(kin2s, kin2)
    p.vals <- c(p.vals, p.value)
    cors <- c(cors, cor.est)
}

cors[is.infinite(cors) & cors < 0] <- min(cors[!is.infinite(cors) & !is.na(cors)])
cors[is.infinite(cors) & cors > 0] <- max(cors[!is.infinite(cors) & !is.na(cors)])

max.cor.est <- max(cors, na.rm=TRUE)
min.cor.est <- min(cors, na.rm=TRUE)
print(c(min.cor.est, max.cor.est))
print(length(which(!is.na(cors))))
cors <- (cors-min.cor.est)/(max.cor.est-min.cor.est)

results <- data.frame(node1=kin1s, node2=kin2s, kinact.cor.p=p.vals)
results$kinact.cor.p <- results$kinact.cor.p/max(results$kinact.cor.p, na.rm=TRUE)
write.table(results[order(results$node1),], "out/reg-site-kinact-full-cor.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

cor.est.results <- data.frame(node1=kin1s, node2=kin2s, kinact.cor.est=cors)
write.table(cor.est.results[order(cor.est.results$node1),],
            "out/reg-site-kinact-full-cor-sign.tsv", quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

