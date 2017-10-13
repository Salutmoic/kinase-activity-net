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

cond.run.ksea <- function(phospho.data, kinase.sites){
    phospho.data <- na.omit(phospho.data)
    phospho.data <- phospho.data[order(phospho.data, decreasing=TRUE)]
    kinases.ksea <- ksea_batchKinases(names(phospho.data), phospho.data,
                                      kinase.sites, trial=1000)
    return(-log10(as.vector(kinases.ksea)))
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

## Load kinase-substrate relationships
kin.sub.tbl <- read.delim("data/human_kinase_table.tsv", as.is=TRUE, sep="\t")
kin.sub.tbl$SUB_MOD_RSD <- sub("[STY]", "", kin.sub.tbl$SUB_MOD_RSD)

## Remove all autophosphorylations, as in Ochoa et al 2016.  The
## suspicion is that autophosphorylation tends to be nonspecific and
## thus will be a source of noise.  We can benchmark how much of an
## effect this has for us.
kin.sub.tbl <- subset(kin.sub.tbl, GENE != SUB_GENE)

kinase.sites <- get.kinase.sites(kin.sub.tbl, phospho.anno)
kin.act <- as.data.frame(cbind(apply(phospho.vals, 2, cond.run.ksea,
                                     kinase.sites)))
rownames(kin.act) <- names(kinase.sites)
colnames(kin.act) <- colnames(phospho.vals)

save(kin.act, "data/log.wGSEA.kinase_condition.clean.Rdata")
