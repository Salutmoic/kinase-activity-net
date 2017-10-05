library(dplyr)
library(reshape2)

cor.prior <- function(cor.vals, is.reg, pscore){
    ## Fisher's r-to-Z
    F.r <- atanh(cor.vals)
    ## Inverse-Gamma hyperparameters for the unknown variance.  Scale
    ## ('b') is larger when it's a regulatory relationship.
    ig.a <- 9
    pscore.mod <- 1/(1+exp(-0.1*pscore))
    if (is.reg){
        ig.b <- 9*pscore.mod
    }else{
        ig.b <- 0.1*pscore.mod
    }
    if (any(ig.b < 0))
        print(ig.b)
    prior <- ((gamma(ig.a+0.5)/gamma(ig.a))*
              ((2*pi*ig.b)^-0.5)*
              ((1+(F.r^2)/(2*ig.b))^(-(ig.a+0.5))))
    return (prior)
}

nfchisq.prior <- function(nfchisq.vals, is.reg, pscore){
    ## Normalized functional Chi^2 statistics are normally distributed
    ig.a <- 9
    pscore.mod <- 1/(1+exp(-0.1*pscore))
    if (is.reg){
        ig.b <- 9*pscore.mod
    }else{
        ig.b <- 0.1*pscore.mod
    }
    prior <- ((gamma(ig.a+0.5)/gamma(ig.a))*
              ((2*pi*ig.b)^-0.5)*
              ((1+(nfchisq.vals^2)/(2*ig.b))^(-(ig.a+0.5))))
    return (prior)
}

is.reg.prior <- function(kin.sub, reg.sites){
    ## Calculate uniform prior probability of a regulatory relationship
    kinases <- read.table("data/external/human-kinome.txt", as.is=TRUE)
    kin.sub.hsap <- subset(kin.sub, KIN_ORGANISM=="human" & SUB_ORGANISM=="human" & SUB_GENE!="" & SUB_GENE %in% kinases$V1,
                           select=c("GENE", "SUB_GENE", "SUB_MOD_RSD"))
    reg.sites.hsap <- subset(reg.sites, ORGANISM=="human" & GENE %in% kinases$V1 & endsWith(MOD_RSD, "-p"),
                             select=c("GENE", "MOD_RSD"))
    colnames(reg.sites.hsap) <- c("SUB_GENE", "SUB_MOD_RSD")
    reg.sites.hsap$SUB_MOD_RSD <- sub("-p", "", reg.sites.hsap$SUB_MOD_RSD)
    reg.kin.sub <- unique(merge(kin.sub.hsap, reg.sites.hsap, all=FALSE))
    reg.kin.sub <- unique(subset(reg.kin.sub, GENE!=SUB_GENE, select=c("GENE", "SUB_GENE")))

    reg.kin.sub.by.sub <- group_by(reg.kin.sub, SUB_GENE)
    reg.kin.per.sub <- summarise(reg.kin.sub.by.sub, num.reg.kin=n())
    mean.regulators <- round(mean(reg.kin.per.sub$num.reg.kin))
    num.char.kinases <- length(unique(c(reg.kin.sub$GENE, reg.kin.sub$SUB_GENE)))

    ## What's the probability, in a pure combinatorics sense, that kinase
    ## A regulates kinase B, given that on average k (mean.regulators) out
    ## of n (num.char.kinases) kinases regulate any given kinase.  So, out
    ## of n_C_k combinations of regulators, what's the probability that
    ## the combination includes kinase A?  (1_C_1 * (n-1)_C_(k-i))/n_C_k.
    ## This simplifies to k/n, btw...
    p.reg <- choose(num.char.kinases-1, mean.regulators-1)/choose(num.char.kinases, mean.regulators)
    return (p.reg)
}

argv <- commandArgs(TRUE)
if (length(argv) != 7){
    stop("USAGE: <script> ASSOC_TYPE ASSOC_FILE PHOS_SCORE_FILE PHOS_BACK_FILE KIN_SUB_FILE REG_SITES_FILE OUT_FILE")
}

assoc.type <- argv[1]
assoc.file <- argv[2]
phos.score.file <- argv[3]
phos.back.file <- argv[4]
kin.sub.file <- argv[5]
reg.sites.file <- argv[6]
out.file <- argv[7]

assoc.tbl <- read.delim(assoc.file, as.is=TRUE)
names(assoc.tbl) <- c("kin1", "kin2", "assoc")

kin.sub <- read.delim(kin.sub.file, as.is=TRUE)
reg.sites <- read.delim(reg.sites.file, as.is=TRUE)

phos.score.tbl <- read.table(phos.score.file, as.is=TRUE)
names(phos.score.tbl) <- c("kin1", "kin2", "phos.score")
phos.score.mean <- mean(phos.score.tbl$phos.score)
phos.score.sd <- sd(phos.score.tbl$phos.score)

dat.tbl <- merge(assoc.tbl, phos.score.tbl)

phos.back.tbl.tmp <- read.table(phos.back.file, as.is=TRUE)
phos.back.tbl <- phos.back.tbl.tmp[2:ncol(phos.back.tbl.tmp)]
phos.back.mean <- mean(as.matrix(phos.back.tbl))
phos.back.sd <- sd(as.matrix(phos.back.tbl))

dat.tbl[which(dat.tbl$kin1==dat.tbl$kin2), "assoc"] <- NA

if (assoc.type %in% c("pcor", "pcor-filter", "scor", "scor-filter")){
    p.assoc.given.reg.pscore <- cor.prior(dat.tbl$assoc, TRUE,
                                          dat.tbl$phos.score)
    p.assoc.given.nreg.pscore <- cor.prior(dat.tbl$assoc, FALSE,
                                           dat.tbl$phos.score)
}else if (assoc.type == "nfchisq"){
    p.assoc.given.reg.pscore <- nfchisq.prior(dat.tbl$assoc, TRUE,
                                              dat.tbl$phos.score)
    p.assoc.given.nreg.pscore <- nfchisq.prior(dat.tbl$assoc, FALSE,
                                               dat.tbl$phos.score)
}

p.reg <- is.reg.prior(kin.sub, reg.sites)

## A bit of a cheat here.  We want the probability of seeing the PSSM
## score for A->B given that A regulates B (or that A doesn't regulate
## B).  However, since "A regulates B" implies "A phosphorylates B",
## and since we don't expect regulatory sites to have significantly
## different PSSM scores from nonregulatory sites, we can just use the
## probability of seeing the A->B PSSM score given that A
## phosphorylates B (or that A doesn't phosphorylate B).
p.pscore.given.reg <- pnorm(dat.tbl$phos.score, mean=phos.score.mean,
                            sd=phos.score.sd)
p.pscore.given.nreg <- pnorm(dat.tbl$phos.score, mean=phos.back.mean,
                                sd=phos.back.sd)

p.cor.pscore.given.reg <- p.assoc.given.reg.pscore * p.pscore.given.reg
p.cor.pscore.given.nreg <- p.assoc.given.nreg.pscore * p.pscore.given.nreg

posterior <- (p.cor.pscore.given.reg*p.reg)/(p.cor.pscore.given.reg*p.reg+p.cor.pscore.given.nreg*(1-p.reg))

posterior.tbl <- data.frame(kin1=dat.tbl$kin1, kin2=dat.tbl$kin2,
                            posterior=posterior)
posterior.tbl$posterior[which(posterior.tbl$kin1==posterior.tbl$kin2)] <- 0.00

write.table(posterior.tbl, out.file, sep="\t", row.names=FALSE,
            col.names=FALSE, quote=FALSE)
