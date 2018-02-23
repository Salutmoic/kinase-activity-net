suppressMessages(library(Biobase))

percent.na <- function(x){
    if (is.null(x) || is.na(length(x)))
        return(0.0)
    return(length(which(is.na(x)))/length(x))
}

choose.redundant.kin <- function(row, kinase.sites, phospho.data){
    kin1 <- row[["kinase1"]]
    no.subs1 <- row[["no.substrates1"]]
    kin2 <- row[["kinase2"]]
    no.subs2 <- row[["no.substrates2"]]
    prop.shared <- row[["proportion.of.substrates.shared"]]
    ## Calculate the difference in amount of missing values between
    ## the two kinases' activity predictions
    kin1.sites <- kinase.sites[[kin1]][which(kinase.sites[[kin1]] %in% rownames(phospho.data))]
    kin2.sites <- kinase.sites[[kin2]][which(kinase.sites[[kin2]] %in% rownames(phospho.data))]
    kin1.data <- phospho.data[kin1.sites,]
    kin2.data <- phospho.data[kin2.sites,]
    kin1.perc.na <- percent.na(kin1.data)
    kin2.perc.na <- percent.na(kin2.data)
    perc.na.diff <- abs(kin1.perc.na-kin2.perc.na)
    if (prop.shared > 0.5){
        ## For small differences in missingness, the kinase with fewer
        ## substrates is redundant
        if (perc.na.diff < 0.1){
            if (no.subs1 == no.subs2){
                ## (unless the #substrates is equal, in which case,
                ## choose by %NA again)
                if (kin1.perc.na > kin2.perc.na)
                    return(kin1)
                else
                    return(kin2)
            }
            if (no.subs1 > no.subs2){
                return(kin2)
            }else{
                return(kin1)
            }
        }else{
            ## otherwise, the kinase with more missing data is
            ## redundant
            if (kin1.perc.na > kin2.perc.na){
                return(kin1)
            }else{
                return(kin2)
            }
        }
    }
}

get.redundant.kins <- function(kin.overlap, kinase.sites, phospho.data,
                               phospho.anno){
    kin.overlap <- kin.overlap[order(kin.overlap$proportion.of.substrates.shared,
                                     decreasing=TRUE),]
    redundant.kins <- c()
    for (i in 1:nrow(kin.overlap)){
        row <- kin.overlap[i,]
        if (row[["kinase1"]] %in% redundant.kins || row[["kinase2"]] %in% redundant.kins)
            next
        redundant.kin <- choose.redundant.kin(row, kinase.sites, phospho.data)
        redundant.kins <- c(redundant.kins, redundant.kin)
    }
    return(redundant.kins)
}

get.kinase.sites <- function(kin.sub.tbl, phospho.anno){
    kinases <- unique(kin.sub.tbl$GENE)
    kinase.sites <- list()
    for (kinase in kinases){
        kin.subs <- subset(kin.sub.tbl, GENE == kinase)
        ## Merge with the phospho data, mainly to convert between HGNC
        ## gene name and Ensembl Protein ID
        kin.subs.ensp <- merge(phospho.anno, kin.subs,
                               by.x=c("gene_name", "positions"),
                               by.y=c("SUB_GENE", "SUB_MOD_RSD"))
        if (is.null(kin.subs.ensp) || nrow(kin.subs.ensp) == 0){
            kinase.sites[[kinase]] <- NA
            next
        }
        ## Construct site IDs as they appear in the phospho data.
        kin.ensp.sites <- paste(kin.subs.ensp$ensp,
                                kin.subs.ensp$positions, "P", sep="_")
        kinase.sites[[kinase]] <- kin.ensp.sites
    }
    return(kinase.sites)
}

cond.do.ks <- function(phospho.data, kinase.sites){
    ensp.sites <- rownames(phospho.data)
    kinases.ks <- c()
    for (kinase in names(kinase.sites)){
        if (is.na(kinase.sites[[kinase]]) ||
            length(which(ensp.sites %in% kinase.sites[[kinase]]))==0){
            kinases.ks <- c(kinases.ks, NA)
            next
        }
        kinase.subs <- phospho.data[which(rownames(phospho.data) %in% kinase.sites[[kinase]])]
        if (length(which(!is.na(kinase.subs))) < 3){
            kinases.ks <- c(kinases.ks, NA)
            next
        }
        other.subs <- phospho.data[which(!(rownames(phospho.data) %in% kinase.sites[[kinase]]))]
        ks <- ks.test(kinase.subs, other.subs, alternative="two.sided")
        median.fc <- median(kinase.subs, na.rm=TRUE)
        median.other.fc <- median(other.subs, na.rm=TRUE)
        fc.sign <- sign(median.fc-median.other.fc)
        kinases.ks <- c(kinases.ks, -log10(ks$p.value)*fc.sign)
    }
    names(kinases.ks) <- names(kinase.sites)
    return(kinases.ks)
}

argv <- commandArgs(TRUE)
if (length(argv) != 2){
    stop("USAGE: <script> COND_NUM USE_AUTOPHOS")
}
## This is a bit of wierdness to accommodate LSF job arrays
cond.num <- as.integer(sub("^\\\\", "", argv[1]))
use.autophos <- as.logical(argv[2])

## argv <- commandArgs(TRUE)
## if (length(argv) != 1)
##     stop("USAGE: <script> USE_AUTOPHOS")
## }
## ## This is a bit of wierdness to accommodate LSF job arrays
## use.autophos <- as.logical(argv[1])

## Load phospho data
load("data/external/esetNR.Rdata")
cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals <- exprs(esetNR)

## Remove CPTAC datasets
cptac.start <- 1604
cptac.end <- 1711
cptac.range <- 1604:1711
cptac.conds <- paste(as.character(cptac.range), "290", sep="_")
phospho.vals <- phospho.vals[, -which(colnames(phospho.vals) %in% cptac.conds)]

## all.kins <- read.table("data/human-kinome.txt", as.is=TRUE)[,1]
## ensp.id.map.full <- read.table("data/ensembl-id-map.tsv", as.is=TRUE)
## names(ensp.id.map.full) <- c("ensembl", "gene.name")
## ensp.id.map <- subset(ensp.id.map.full, gene.name %in% all.kins)
## rownames(ensp.id.map) <- ensp.id.map$ensembl

## phospho.vals.sub.ensp <- unlist(lapply(strsplit(rownames(phospho.vals.sub), split="_"),
##                                        function(foo) foo[1]))
## phospho.vals.sub.gene <- ensp.id.map[phospho.vals.sub.ensp, "gene.name"]
## phospho.vals <- phospho.vals.sub[which(phospho.vals.sub.ensp %in% ensp.id.map$ensembl),]
## phospho.vals.ensp <- unlist(lapply(strsplit(rownames(phospho.vals), split="_"),
##                                    function(foo) foo[1]))
## phospho.vals.pos <- unlist(lapply(strsplit(rownames(phospho.vals), split="_"),
##                                   function(foo) foo[2]))
## phospho.vals.gene <- ensp.id.map[phospho.vals.ensp, "gene.name"]
## rownames(phospho.vals) <- paste(phospho.vals.gene, phospho.vals.pos, sep="_")

if (cond.num > ncol(phospho.vals)){
    stop(paste("Invalid condition number.  Max # =", ncol(phospho.vals)))
}

## Load kinase-substrate relationships
kin.sub.tbl <- read.delim("data/psiteplus-kinase-substrates.tsv", as.is=TRUE,
                          sep="\t")
kin.sub.tbl$SUB_MOD_RSD <- sub("[STY]", "", kin.sub.tbl$SUB_MOD_RSD)

## Remove all autophosphorylations, as in Ochoa et al 2016.  The
## suspicion is that autophosphorylation tends to be nonspecific and
## thus will be a source of noise.  We can benchmark how much of an
## effect this has for us.
if (!use.autophos){
    kin.sub.tbl <- subset(kin.sub.tbl, GENE != SUB_GENE)
}

## ## Get all the sites annotated as substrates for each kinase.
kinase.sites <- get.kinase.sites(kin.sub.tbl, phospho.anno)

## kin.overlap <- read.delim("data/psiteplus-kinase-substrate-overlap.tsv", as.is=TRUE)
## ## Remove redundant kinases
## redundant.kins <- get.redundant.kins(kin.overlap, kinase.sites, phospho.vals,
##                                      phospho.anno)
## good.kins <- setdiff(kin.sub.tbl$GENE, redundant.kins)
## kin.sub.tbl <- subset(kin.sub.tbl, !(GENE %in% redundant.kins))
## ## Remove ambiguous sites
## all.subs <- paste(kin.sub.tbl$SUB_GENE, kin.sub.tbl$SUB_MOD_RSD, sep="_")
## dup.subs <- all.subs[which(duplicated(all.subs))]
## good.subs <- which(!(all.subs %in% dup.subs))
## print(length(good.subs))
## kin.sub.tbl <- kin.sub.tbl[good.subs,]

## ## Get all the sites annotated as substrates for each kinase, again.
## kinase.sites <- get.kinase.sites(kin.sub.tbl, phospho.anno)

## Estimate kinase activities for this condition.
kin.act <- cond.do.ks(subset(phospho.vals,
                                select=colnames(phospho.vals)[cond.num]),
                      kinase.sites)

## Reformat as a narrow table
kin.act.tbl <- data.frame(kinase=names(kin.act),
                          cond=rep(colnames(phospho.vals)[cond.num], length(kin.act)),
                          kin.act=kin.act)

if (use.autophos){
    out.dir <- "tmp/ks-autophos/"
}else{
    out.dir <- "tmp/ks/"
}

write.table(kin.act.tbl,
            paste0(out.dir, "cond-", cond.num, ".tsv"),
            sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
