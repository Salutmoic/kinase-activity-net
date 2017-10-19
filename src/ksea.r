library(ksea)
suppressMessages(library(Biobase))

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
        if (is.null(kin.subs.ensp) || nrow(kin.subs.ensp) == 0)
            next
        ## Construct site IDs as they appear in the phospho data.
        kin.ensp.sites <- paste(kin.subs.ensp$ensp, kin.subs.ensp$positions,
                                "P", sep="_")
        kinase.sites[[kinase]] <- kin.ensp.sites
    }
    return(kinase.sites)
}

cond.run.ksea <- function(phospho.data, kinase.sites, num.trials){
    phospho.data <- na.omit(phospho.data[,1])
    phospho.data <- sort(phospho.data, decreasing=TRUE)
    ensp.sites <- names(phospho.data)
    kinases.ksea <- ksea_batchKinases(ensp.sites, phospho.data,
                                      kinase.sites, trial=num.trials)
    for (kinase in names(kinase.sites)){
        if (length(which(ensp.sites %in% kinase.sites[[kinase]]))==0)
            next
        kinase.subs <- phospho.data[which(names(phospho.data) %in% kinase.sites[[kinase]])]
        median.fc <- median(kinase.subs)
        fc.sign <- sign(median.fc)
        ksea.i <- which(grepl(paste0("^", kinase, "\\."), names(kinases.ksea)))
        if (length(ksea.i)==0)
            next
        names(kinases.ksea)[ksea.i] <- kinase
        if (is.na(kinases.ksea[ksea.i])){
            next
        }
        if (kinases.ksea[ksea.i] == 0){
            kinases.ksea[ksea.i] <- 1/num.trials
        }
        kinases.ksea[ksea.i] <- -log10(kinases.ksea[ksea.i])*fc.sign
    }
    return(kinases.ksea)
}

argv <- commandArgs(TRUE)
if (!(length(argv) %in% c(2, 3))){
    stop("USAGE: <script> COND_NUM NUM_TRIALS USE_AUTOPHOS")
}
## This is a bit of wierdness to accommodate LSF job arrays
cond.num <- as.integer(sub("^\\\\", "", argv[1]))
num.trials <- as.integer(argv[2])
if (length(argv) == 3){
    use.autophos <- as.logical(argv[3])
    if (is.na(use.autophos))
        use.autophos <- FALSE
}else{
    use.autophos <- FALSE
}

## Setup multicore forking for mclapply based in bsub available cores
hosts <- Sys.getenv("LSB_MCPU_HOSTS")
if(length(hosts) >0){
    ncores <- unlist(strsplit(hosts," "))[2]
}else{
    ncores <- 1
}
options("mc.cores" = ncores)

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

if (cond.num > ncol(phospho.vals)){
    stop(paste("Invalid condition number.  Max # =", ncol(phospho.vals)))
}

## Load kinase-substrate relationships
kin.sub.tbl <- read.delim("data/human_kinase_table.tsv", as.is=TRUE, sep="\t")
kin.sub.tbl$SUB_MOD_RSD <- sub("[STY]", "", kin.sub.tbl$SUB_MOD_RSD)

## Remove all autophosphorylations, as in Ochoa et al 2016.  The
## suspicion is that autophosphorylation tends to be nonspecific and
## thus will be a source of noise.  We can benchmark how much of an
## effect this has for us.
if (!use.autophos){
    kin.sub.tbl <- subset(kin.sub.tbl, GENE != SUB_GENE)
}

## Get all the sites annotated as substrates for each kinase.
kinase.sites <- get.kinase.sites(kin.sub.tbl, phospho.anno)

## Estimate kinase activities for this condition.
kin.act <- cond.run.ksea(subset(phospho.vals, select=colnames(phospho.vals)[cond.num]),
                         kinase.sites, num.trials)

## Reformat as a narrow table
kin.act.tbl <- data.frame(kinase=names(kin.act),
                          cond=rep(colnames(phospho.vals)[cond.num], length(kin.act)),
                          kin.act=kin.act)

if (use.autophos){
    out.dir <- "tmp/ksea-autophos/"
}else{
    out.dir <- "tmp/ksea/"
}

write.table(kin.act.tbl,
            paste0(out.dir, "cond-", cond.num, ".tsv"),
            sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
