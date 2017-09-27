library(dplyr)
library(reshape2)

cor.prior <- function(cor.vals, is.reg){
    ## Fisher's r-to-Z
    F.r <- atanh(cor.vals)
    ## Inverse-Gamma hyperparameters for the unknown variance.  Scale
    ## ('b') is larger when it's a regulatory relationship.
    ig.a <- 9
    if (is.reg){
        ig.b <- 9
     }else{
         ig.b <- 0.1
    }
    prior <- ((gamma(ig.a+0.5)/gamma(ig.a))*
              ((2*pi*ig.b)^-0.5)*
              ((1+(F.r^2)/(2*ig.b))^(-(ig.a+0.5))))
    return (prior)
}

nfchisq.prior <- function(nfchisq.vals, is.reg){
    ## Normalized functional Chi^2 statistics are normally distributed
    if (is.reg){
        ig.b <- 9
     }else{
         ig.b <- 0.1
    }
    prior <- ((gamma(ig.a+0.5)/gamma(ig.a))*
              ((2*pi*ig.b)^-0.5)*
              ((1+(nfchisq.vals^2)/(2*ig.b))^(-(ig.a+0.5))))
    return (prior)
}

is.reg.prior <- function(){
    ## Calculate uniform prior probability of a regulatory relationship
    kinases <- read.table("~/nfs-home/data/annotation/human-kinome.txt", as.is=TRUE)
    kin.sub <- read.delim("~/nfs-home/data/proteomes/PTMs/phosphosite-plus/Kinase_Substrate_Dataset-2016-04.txt", as.is=TRUE)
    reg.sites <- read.delim("~/nfs-home/data/proteomes/PTMs/phosphosite-plus/Regulatory_sites-2016-02-27.tsv", as.is=TRUE)
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
if (length(argv) != 3){
    stop("USAGE: <script> PRIOR_TYPE DATA_FILE OUT_FILE")
}

prior.type <- argv[1]
data.file <- argv[2]
out.file <- argv[3]

dat.tbl <- read.delim(data.file, as.is=TRUE)
names(dat.tbl) <- c("node1", "node2", "assoc")

dat.tbl[which(dat.tbl$node1==dat.tbl$node2), "assoc"] <- NA

if (prior.type %in% c("pcor", "pcor-filter", "scor", "scor-filter")){
    dat.prior.given.reg <- cor.prior(dat.tbl$assoc, TRUE)
    dat.prior.given.nreg <- cor.prior(dat.tbl$assoc, FALSE)
}else if (prior.type == "nfchisq"){
    dat.prior.given.reg <- nfchisq.prior(dat.tbl$assoc, TRUE)
    dat.prior.given.nreg <- nfchisq.prior(dat.tbl$assoc, FALSE)
}

p.reg <- is.reg.prior()

posterior <- (dat.prior.given.reg*p.reg)/(dat.prior.given.reg*p.reg+dat.prior.given.nreg*(1-p.reg))

posterior.tbl <- data.frame(node1=dat.tbl$node1, node2=dat.tbl$node2,
                            posterior=posterior)
posterior.tbl$posterior[which(posterior.tbl$node1==posterior.tbl$node2)] <- 0.00

write.table(posterior.tbl, out.file, sep="\t", row.names=FALSE,
            col.names=FALSE, quote=FALSE)
