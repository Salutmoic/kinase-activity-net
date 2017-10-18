suppressMessages(library(reshape2))
suppressMessages(library(VIM))
suppressMessages(library(entropy))
## library(FNN)
suppressMessages(library(Biobase))

filter.act.preds <- function(kin.act, phospho.anno, phospho.vals, min.sites){
    for (kinase in rownames(kin.act)){
        kin.subs <- subset(kin.sub.tbl, GENE==kinase)
        kin.subs.ensp <- merge(phospho.anno, kin.subs,
                               by.x=c("gene_name", "positions"),
                               by.y=c("SUB_GENE", "SUB_MOD_RSD"))
        ensp.sites <- paste(kin.subs.ensp$ensp, kin.subs.ensp$positions, "P",
                            sep="_")
        ensp.sites <- intersect(rownames(phospho.vals), ensp.sites)
        if (length(ensp.sites) == 0){
            next
        }
        if (length(ensp.sites) == 1){
            kin.act[kinase, 1:ncol(kin.act)] <- NA
            next
        }
        sites.detected <- apply(phospho.vals[ensp.sites,], 2,
                                function(cond) length(which(!is.na(cond))))
        bad.conds <- setdiff(colnames(phospho.vals)[which(sites.detected<min.sites)],
                             colnames(kin.act)[which(is.na(kin.act[kinase,]))])
        kin.act[kinase, bad.conds] <- NA
        ## if (length(bad.conds) > 0){
        ##     message(paste(length(bad.conds), "conditions omitted for", kinase,
        ##                   "due to lack of evidence."))
        ## }
    }
    return(kin.act)
}

choose.redundant.kin <- function(row){
    kin1 <- row[1]
    no.subs1 <- as.integer(row[2])
    kin2 <- row[3]
    no.subs2 <- as.integer(row[4])
    prop.shared <- as.double(row[6])
    ## Calculate the difference in amount of missing values between
    ## the two kinases' activity predictions
    kin1.perc.na <- percent.na(kin.act[kin1,])
    kin2.perc.na <- percent.na(kin.act[kin2,])
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

get.redundant.kins <- function(kin.overlap){
    kin.overlap <- kin.overlap[order(kin.overlap$proportion.of.substrates.shared,
                                     decreasing=TRUE),]
    redundant.kins <- c()
    for (i in 1:nrow(kin.overlap)){
        row <- unlist(kin.overlap[i,])
        if (row[1] %in% redundant.kins || row[2] %in% redundant.kins)
            next
        redundant.kin <- choose.redundant.kin(row)
        redundant.kins <- c(redundant.kins, redundant.kin)
    }
    return(redundant.kins)
}

percent.na <- function(x){
    if (is.null(x) || is.na(length(x)))
        return(0.0)
    return(length(which(is.na(x)))/length(x))
}

make.balanced.table <- function(full.tbl, na.threshold){
    ## Iteratively remove the worst row or column until all rows and
    ## cols have less than the threshold %NA.
    while (nrow(full.tbl) > 1 && ncol(full.tbl) > 1){
        max.row.na <- max(apply(full.tbl, 1, percent.na))
        max.col.na <- max(apply(full.tbl, 2, percent.na))
        if (max.row.na > max.col.na){
            max.na <- which.max(apply(full.tbl, 1, percent.na))
            worst.row <- rownames(full.tbl)[max.na]
            new.rows <- setdiff(rownames(full.tbl), worst.row)
            full.tbl <- full.tbl[new.rows,]
        }else{
            max.na <- which.max(apply(full.tbl, 2, percent.na))
            worst.col <- colnames(full.tbl)[max.na]
            new.cols <- setdiff(colnames(full.tbl), worst.col)
            full.tbl <- subset(full.tbl, select=new.cols)
        }
        if (all(apply(full.tbl, 1, function(tbl.row) percent.na(tbl.row) < na.threshold)) &&
            all(apply(full.tbl, 2, function(tbl.col) percent.na(tbl.col) < na.threshold))){
            return(full.tbl)
        }
    }
    return(NULL)
}

make.max.row.table <- function(full.tbl, na.threshold, min.cols){
    ## 1) Get rid of all columns that have more than the threshold %NA
    ## 2) Try removing more columns until all of the rows have less than the threshold
    ## 3) If all the rows pass and we haven't removed too many columns, return the table
    ## 4) Otherwise, remove the worst row and repeat from #2
    bad.cols <- colnames(full.tbl)[which(apply(full.tbl, 2,
                                               function(tbl.col){
                                                   percent.na(tbl.col) > na.threshold
                                               }))]
    good.cols <- setdiff(colnames(full.tbl), bad.cols)
    sub.tbl <- subset(full.tbl, select=good.cols)
    while (nrow(sub.tbl) > 1){
        tmp.tbl <- sub.tbl
        fail <- FALSE
        while (any(apply(tmp.tbl, 1, function(tbl.row) percent.na(tbl.row) > na.threshold))){
            max.na.col <- which.max(apply(tmp.tbl, 2, percent.na))
            good.cols <- setdiff(colnames(tmp.tbl), colnames(tmp.tbl)[max.na.col])
            if (is.null(good.cols) || length(good.cols) <= 1){
                fail <- TRUE
                break
            }
            tmp.tbl <- subset(tmp.tbl, select=good.cols)
        }
        if (!fail && ncol(tmp.tbl) >= min.cols)
            return(tmp.tbl)
        max.na.row <- which.max(apply(sub.tbl, 1, percent.na))
        good.rows <- setdiff(rownames(sub.tbl), rownames(sub.tbl)[max.na.row])
        if (is.null(good.rows) || length(good.rows) == 0)
            return (NULL)
        sub.tbl <- sub.tbl[good.rows,]
    }
    return (NULL)
}

make.max.col.table <- function(full.tbl, na.threshold, min.rows){
    ## 1) Get rid of all rows that have more than the threshold %NA
    ## 2) Try removing more rows until all of the columns have less than the threshold
    ## 3) If all the columns pass and we haven't removed too many rows, return the table
    ## 4) Otherwise, remove the worst column and repeat from #2
    bad.rows <- rownames(full.tbl)[which(apply(full.tbl, 1,
                                               function(tbl.row){
                                                   percent.na(tbl.row) > na.threshold
                                               }))]
    good.rows <- setdiff(rownames(full.tbl), bad.rows)
    sub.tbl <- full.tbl[good.rows,]

    while (ncol(sub.tbl) > 1){
        tmp.tbl <- sub.tbl
        fail <- FALSE
        while (any(apply(tmp.tbl, 2, percent.na) > na.threshold)){
            max.na.row <- which.max(apply(tmp.tbl, 1, percent.na))
            good.rows <- setdiff(rownames(tmp.tbl), rownames(tmp.tbl)[max.na.row])
            if (is.null(good.rows) || length(good.rows) == 0){
                fail <- TRUE
                break
            }
            tmp.tbl <- tmp.tbl[good.rows,]
        }
        if (!fail && nrow(tmp.tbl) >= min.rows)
            return(tmp.tbl)
        max.na.col <- which.max(apply(sub.tbl, 2, percent.na))
        good.cols <- setdiff(colnames(sub.tbl), colnames(sub.tbl)[max.na.col])
        if (is.null(good.cols) || length(good.cols) == 0)
            return (NULL)
        sub.tbl <- subset(sub.tbl, select=good.cols)
    }
    return (NULL)}

filter.low.activity <- function(full.tbl){
    bad.rows <- which(apply(full.tbl, 1, function(tbl.row) all(abs(tbl.row) < 1.0)))
    if (length(bad.rows) == 0){
        return(full.tbl)
    }
    good.rows <- setdiff(rownames(full.tbl, bad.rows))
    if (length(good.rows) < nrow(full.tbl)){
        message(paste("Filtering out", nrow(full.tbl)-length(good.rows), "rows due to low activity"))
    }
    filtered.tbl <- full.tbl[good.rows,]
    return(filtered.tbl)
}

filter.by.entropy <- function(full.tbl, entropy.filter.stren, num.bins){
    tbl.min <- min(full.tbl, na.rm=TRUE)
    tbl.max <- max(full.tbl, na.rm=TRUE)
    col.entropy <- apply(full.tbl, 2,
                          function(col){
                              entropy(discretize(na.omit(col),
                                                 numBins=num.bins,
                                                 r=range(tbl.min, tbl.max)))})
    good.cols <- which(col.entropy > entropy.filter.stren*max(col.entropy, na.rm=TRUE))
    if (length(good.cols) < ncol(full.tbl)){
        message(paste("Filtering out", ncol(full.tbl)-length(good.cols), "columns due to low entropy"))
    }
    filtered.tbl <- subset(full.tbl, select=good.cols)
    row.entropy <- apply(filtered.tbl, 1,
                         function(row){
                             entropy(discretize(na.omit(row),
                                                numBins=num.bins,
                                                r=range(tbl.min, tbl.max)))})
    good.rows <- which(row.entropy > entropy.filter.stren*max(row.entropy, na.rm=TRUE))
    if (length(good.rows) < nrow(full.tbl)){
        message(paste("Filtering out", nrow(full.tbl)-length(good.rows), "rows due to low entropy"))
    }
    filtered.tbl <- filtered.tbl[good.rows,]
    return(filtered.tbl)
}

write.activity.table <- function(tbl, file.name){
    tbl.df <- melt(tbl)
    ## lame hack....
    if (ncol(tbl.df)!=3){
        tbl.tmp <- tbl
        tbl.tmp$kinase <- rownames(tbl)
        tbl.df <- melt(tbl.tmp)
    }
    names(tbl.df) <- c("kinase", "condition", "activity")
    tbl.df <- tbl.df[order(tbl.df$kinase, tbl.df$condition),]
    write.table(tbl.df, file.name, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

impute.on.table <- function(tbl){
    tbl.imp <- data.frame(kNN(tbl, imp_var=FALSE))
    rownames(tbl.imp) <- rownames(tbl)
    colnames(tbl.imp) <- colnames(tbl)
    tbl.imp$kinase <- rownames(tbl.imp)
    return(tbl.imp)
}

print.table.info <- function(tbl){
    message(paste0("dims=", dim(tbl)[1], " x ", dim(tbl)[2], "\n",
                   "%NA=", percent.na(as.matrix(tbl))))
}

argv <- commandArgs(TRUE)
if (length(argv) != 4){
    stop("USAGE: <script> STRATEGY MIN_SITES ENTROPY_FILTER NA_THRESHOLD")
}
strategy <- argv[1]
min.sites <- as.integer(argv[2])
entropy.filter.stren <- as.numeric(argv[3])
na.threshold <- as.numeric(argv[4])

## General parameters
entropy.bins <- 10

## Load data
load("data/external/esetNR.Rdata")
load("data/log.wKSEA.kinase_condition.clean.Rdata")
kin.act <- log.wKSEA.kinase_condition.clean
cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals <- exprs(esetNR)
kin.sub.tbl <- read.delim("data/human_kinase_table.tsv", as.is=TRUE, sep="\t")
kin.sub.tbl$SUB_MOD_RSD <- sub("[STY]", "", kin.sub.tbl$SUB_MOD_RSD)
kin.overlap <- read.delim("data/kinase_substrate_overlap.tsv", as.is=TRUE)

## Please note also we added the breast cancer CPTAC data from two different
## sources: The CPTAC database as well as from the supplementary material of
## the paper (which got published later). I would recommend you to use the
## ones from the supplementary material of the paper marked as ³supp² and
## avoid the ones from CPTAC database (³1604_290² to ³1711_290²).
cptac.start <- 1604
cptac.end <- 1711
cptac.range <- 1604:1711
cptac.conds <- paste(as.character(cptac.range), "290", sep="_")
phospho.vals <- phospho.vals[, -which(colnames(phospho.vals) %in% cptac.conds)]

## TODO: remove?
## Only choose kinases that have at least 10 substrates
## kin.overlap.sub <- subset(kin.overlap, no.substrates1 >= min.sites & no.substrates2 >= min.sites)
kin.overlap.sub <- kin.overlap
good.kins <- unique(c(kin.overlap.sub$kinase1, kin.overlap.sub$kinase2))

## Remove redundant kinases
redundant.kins <- get.redundant.kins(kin.overlap.sub)

kin.act <- kin.act[setdiff(good.kins, redundant.kins),]

## Filter out kinase activity predictions that were not based on
## enough sites
kin.act <- filter.act.preds(kin.act, phospho.anno, phospho.vals, min.sites)

## We look for "empty" rows here...not really empty, but anything
## that's practically empty (>90% missing values).
empty.rows <- rownames(kin.act)[which(apply(kin.act, 1,
                                            function(tbl.row){
                                                percent.na(tbl.row)>0.9
                                            }))]
kin.act <- kin.act[setdiff(rownames(kin.act), empty.rows),]
empty.cols <- colnames(kin.act)[which(apply(kin.act, 2,
                                            function(tbl.col){
                                                percent.na(tbl.col)>0.9
                                            }))]
kin.act <- kin.act[,setdiff(colnames(kin.act), empty.cols)]

## Filter out kinases that have overall low activity predictions.
kin.act.filt <- filter.low.activity(kin.act)

## Build the optimized table according to the desired strategy
if (strategy=="balanced"){
    max.table <- make.balanced.table(kin.act.filt, na.threshold)
}else if (strategy=="max-rows"){
    max.table <- make.max.row.table(kin.act.filt, na.threshold, 50)
}else if (strategy=="max-cols"){
    max.table <- make.max.col.table(kin.act.filt, na.threshold, 20)
}
if (is.null(max.table))
    stop("Failed to optimize table.")

## Remove any kinases or conditions that do not show enough
## variability in their activity predictions.
max.table.filt <- filter.by.entropy(max.table, entropy.filter.stren,
                                    entropy.bins)

## Output the tble
print.table.info(max.table.filt)
write.activity.table(max.table.filt, paste0("data/kinact-", strategy, ".tsv"))

## Do imputation and output the imputed table
max.table.imp <- impute.on.table(max.table.filt)
write.activity.table(max.table.imp, paste0("data/kinact-", strategy, "-imp.tsv"))
