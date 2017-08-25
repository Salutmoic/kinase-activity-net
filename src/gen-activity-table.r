library(reshape2)
library(VIM)
library(entropy)
## library(FNN)
suppressMessages(library(Biobase))
load("data/esetNR.Rdata")
load("data/log.wGSEA.kinase_condition.clean.Rdata")

kin.act <- log.wGSEA.kinase_condition.clean

cond.anno <- pData(esetNR)
phospho.anno <- fData(esetNR)
phospho.vals <- exprs(esetNR)

## Please note also we added the breast cancer CPTAC data from two different
## sources: The CPTAC database as well as from the supplementary material of
## the paper (which got published later). I would recommend you to use the
## ones from the supplementary material of the paper marked as ³supp² and
## avoid the ones from CPTAC database (³1604_290² to ³1711_290²).

## cptac.start <- which(colnames(phospho.vals)=="1604_290")
## cptac.end <- which(colnames(phospho.vals)=="1711_290")
cptac.start <- 1604
cptac.end <- 1711
cptac.range <- 1604:1711
cptac.conds <- paste(as.character(cptac.range), "290", sep="_")
phospho.vals <- phospho.vals[, -which(colnames(phospho.vals) %in% cptac.conds)]

cptac.start <- which(colnames(phospho.vals)=="1604_290")
cptac.end <- which(colnames(phospho.vals)=="1711_290")
kin.act <- kin.act[,-which(colnames(kin.act) %in% cptac.conds)]

kin.overlap <- read.delim("data/kinase_substrate_overlap.tsv", as.is=TRUE)

kin.subs.a <- kin.overlap[,c("kinase1", "no.substrates")]
kin.subs.b <- kin.overlap[,c("kinase2", "no.substrates.1")]
colnames(kin.subs.b) <- colnames(kin.subs.a)
kin.subs <- unique(rbind(kin.subs.a, kin.subs.b))
good.kins <- kin.subs$kinase1[which(kin.subs$no.substrates>10)]

kin.overlap.sub <- subset(kin.overlap, kinase1 %in% good.kins & kinase2 %in% good.kins)
redundant.kins <- unique(unlist(apply(kin.overlap.sub, 1,
                                      function(row){
                                          kin1 <- row[1]
                                          no.subs1 <- row[2]
                                          kin2 <- row[3]
                                          no.subs2 <- row[4]
                                          prop.shared <- row[6]
                                          if (prop.shared > 0.5){
                                              if (no.subs1 >= no.subs2){
                                                  return (kin2)
                                              }else{
                                                  return (kin1)
                                              }
                                          }
                                      })))

kin.act <- kin.act[which(rownames(kin.act) %in% good.kins & (!rownames(kin.act) %in% redundant.kins)),]

percent.na <- function(x) return(length(which(is.na(x)))/length(x))

make.max.table <- function(full.tbl, na.threshold){
    ## Sort the columns by number of missing columns.  Then take one
    ## column at a time until one row has >threshold missing data.  At
    ## that point, stop, revert to previous iteration.
    sorted.tbl <- full.tbl[,order(apply(full.tbl, 2,
                                        function(col) {
                                            length(which(!is.na(col)))
                                        }),
                                  decreasing=TRUE)]
    bad.cols <- c()
    for (i in 1:ncol(sorted.tbl)){
        if (percent.na(sorted.tbl[,i]) > na.threshold){
            bad.cols <- c(bad.cols, colnames(sorted.tbl)[i])
        }
    }
    good.cols <- setdiff(colnames(sorted.tbl), bad.cols)
    good.sorted.tbl <- subset(sorted.tbl, select=good.cols)

    num.rows <- nrow(good.sorted.tbl)
    while (num.rows > 1){
        for (i in ncol(good.sorted.tbl):2){
            tmp.tbl <- good.sorted.tbl[,1:i]
            if (all(apply(tmp.tbl, 1, function(row) percent.na(row) < na.threshold))){
                return(tmp.tbl)
            }
        }
        worst.row <- rownames(good.sorted.tbl)[which.max(apply(good.sorted.tbl, 1, percent.na))]
        new.rows <- setdiff(rownames(good.sorted.tbl), worst.row)
        good.sorted.tbl <- good.sorted.tbl[new.rows,]
        num.rows <- num.rows - 1
    }
}

make.max.row.table <- function(full.tbl, na.threshold, min.cols){
    ## Sort the columns by number of missing columns.  Then take one
    ## column at a time until one row has >threshold missing data.  At
    ## that point, stop, revert to previous iteration.
    sorted.tbl <- full.tbl[,order(apply(full.tbl, 2,
                                        function(col) {
                                            length(which(!is.na(col)))
                                        }),
                                  decreasing=TRUE)]
    bad.cols <- c()
    for (i in 1:ncol(sorted.tbl)){
        if (percent.na(sorted.tbl[,i]) > na.threshold){
            bad.cols <- c(bad.cols, colnames(sorted.tbl)[i])
        }
    }
    good.cols <- setdiff(colnames(sorted.tbl), bad.cols)
    good.sorted.tbl <- subset(sorted.tbl, select=good.cols)

    num.rows <- nrow(good.sorted.tbl)
    while(num.rows > 1){
        for (i in 1:ncol(good.sorted.tbl)){
            if (i==1){
                final.tbl <- as.data.frame(good.sorted.tbl[,i],
                                           col.names=colnames(good.sorted.tbl)[i],
                                           row.names=rownames(good.sorted.tbl))
                colnames(final.tbl) <- colnames(good.sorted.tbl)[i]
            }else{
                tmp.tbl <- cbind(final.tbl, good.sorted.tbl[,i])
                if (all(apply(tmp.tbl, 1, function(row) percent.na(row) < na.threshold))){
                    cnames <- colnames(final.tbl)
                    final.tbl <- cbind(final.tbl, good.sorted.tbl[,i])
                    colnames(final.tbl) <- c(cnames, colnames(good.sorted.tbl)[i])
                }else{
                    next
                }
            }
        }
        if (ncol(final.tbl) >= min.cols){
            return(final.tbl)
        }
        worst.row <- rownames(good.sorted.tbl)[which.max(apply(good.sorted.tbl[,1:min.cols], 1, percent.na))]
        new.rows <- setdiff(rownames(good.sorted.tbl), worst.row)
        good.sorted.tbl <- good.sorted.tbl[new.rows,]
        num.rows <- num.rows - 1
    }
}

make.max.col.table <- function(full.tbl, na.threshold, min.rows){
    ## Sort the columns by number of missing columns.  Then take one
    ## column at a time until one row has >threshold missing data.  At
    ## that point, stop, revert to previous iteration.
    sorted.tbl <- full.tbl[order(apply(full.tbl, 1,
                                        function(col) {
                                            length(which(!is.na(col)))
                                        }),
                                 decreasing=TRUE),]
    bad.rows <- c()
    for (i in 1:nrow(sorted.tbl)){
        if (percent.na(sorted.tbl[i,]) > na.threshold){
            bad.rows <- c(bad.rows, rownames(sorted.tbl)[i])
        }
    }
    good.rows <- setdiff(rownames(sorted.tbl), bad.rows)
    good.sorted.tbl <- sorted.tbl[good.rows,]

    num.cols <- ncol(good.sorted.tbl)
    while(num.cols > 1){
        for (i in 1:nrow(good.sorted.tbl)){
            if (i==1){
                final.tbl <- as.data.frame(t(good.sorted.tbl[i,]),
                                           row.names=rownames(good.sorted.tbl)[i],
                                           col.names=colnames(good.sorted.tbl))
            }else{
                tmp.tbl <- rbind(final.tbl, good.sorted.tbl[i,])
                if (all(apply(tmp.tbl, 2, function(col) percent.na(col) < na.threshold))){
                    rnames <- rownames(final.tbl)
                    final.tbl <- rbind(final.tbl, good.sorted.tbl[i,])
                    rownames(final.tbl) <- c(rnames, rownames(good.sorted.tbl)[i])
                }else{
                    next
                }
            }
        }
        if (nrow(final.tbl) >= min.rows){
            return(final.tbl)
        }
        worst.col <- colnames(good.sorted.tbl)[which.max(apply(good.sorted.tbl[1:min.rows,], 2, percent.na))]
        new.cols <- setdiff(colnames(good.sorted.tbl), worst.col)
        good.sorted.tbl <- good.sorted.tbl[,new.cols]
        num.cols <- num.cols - 1
    }
}

filter.low.activity <- function(full.tbl){
    bad.rows <- which(apply(full.tbl, 1, function(row) all(abs(row) < 1.0)))
    if (length(bad.rows) == 0){
        return(full.tbl)
    }
    good.rows <- setdiff(rownames(full.tbl, bad.rows))
    if (length(good.rows) < nrow(full.tbl)){
        message(paste("Filtering out", nrow(full.tbl)-length(good.rows), "rows due to low entropy"))
    }
    filtered.tbl <- full.tbl[good.rows,]
    return(filtered.tbl)
}

filter.by.entropy <- function(full.tbl, entropy.threshold, num.bins){
    tbl.min <- min(full.tbl, na.rm=TRUE)
    tbl.max <- max(full.tbl, na.rm=TRUE)
    good.cols <- which(apply(full.tbl, 2,
                             function(col){
                                 entropy(discretize(na.omit(col),
                                                    numBins=num.bins,
                                                    r=range(tbl.min, tbl.max))) > entropy.threshold
                             }))
    if (length(good.cols) < ncol(full.tbl)){
        message(paste("Filtering out", ncol(full.tbl)-length(good.cols), "columns due to low entropy"))
    }
    filtered.tbl <- subset(full.tbl, select=good.cols)
    good.rows <- rownames(filtered.tbl)[which(apply(filtered.tbl, 1,
                                                    function(row){
                                                        entropy(discretize(na.omit(row),
                                                                           numBins=num.bins,
                                                                           r=range(tbl.min, tbl.max))) > entropy.threshold
                                                    }))]
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

impute.on.table <- function(tbl, filter.entropy=TRUE){
    tbl.imp <- data.frame(kNN(tbl, imp_var=FALSE))
    rownames(tbl.imp) <- rownames(tbl)
    colnames(tbl.imp) <- colnames(tbl)
    tbl.imp$kinase <- rownames(tbl.imp)
    return(tbl.imp)
}

print.table.info <- function(tbl){
    message(paste0("dims=", dim(tbl)[1], " x ", dim(tbl)[2], "\n%NA=", percent.na(as.matrix(tbl))))
}

na.threshold <- 1/3
entropy.threshold <- 1.0
entropy.bins <- 30

kin.act.filt <- filter.low.activity(kin.act)

message("Creating table with maximum rows and columns")
max.table <- make.max.table(kin.act.filt, na.threshold)
max.table.filt <- filter.by.entropy(max.table, entropy.threshold, entropy.bins)
print.table.info(max.table.filt)
write.activity.table(max.table.filt, "kinase-activity-max-rows-max-cols.tsv")
max.table.imp <- impute.on.table(max.table.filt)
write.activity.table(max.table.imp, "kinase-activity-max-rows-max-cols-imp.tsv")

message("Creating table with maximum rows")
max.row.table <- make.max.row.table(kin.act.filt, na.threshold, 20)
max.row.table.filt <- filter.by.entropy(max.row.table, entropy.threshold, entropy.bins)
print.table.info(max.row.table)
write.activity.table(max.row.table, "kinase-activity-max-rows.tsv")
max.row.table.imp <- impute.on.table(max.row.table)
write.activity.table(max.row.table.imp, "kinase-activity-max-rows-imp.tsv")

message("Creating table with maximum columns")
max.col.table <- make.max.col.table(kin.act.filt, na.threshold, 20)
max.col.table.filt <- filter.by.entropy(max.col.table, entropy.threshold, entropy.bins)
print.table.info(max.col.table)
write.activity.table(max.col.table, "kinase-activity-max-cols.tsv")
max.col.table.imp <- impute.on.table(max.col.table)
write.activity.table(max.col.table.imp, "kinase-activity-max-cols-imp.tsv")

## egf.kinases <- c("AKT1", "AKT2", "ARAF", "BRAF", "CDK2", "CDK4",
##                  "CHUK", "GSK3A", "GSK3B", "JAK2", "MAP2K1", "MAP2K2",
##                  "MAPK1", "MAPK3", "MTOR", "PAK1", "PDPK1", "PRKAA2", "PRKACA",
##                  "PTK2", "RAF1", "RPS6KB1", "RPS6KB2", "SGK1", "SRC")
##### Only kinases w/ perturbations
egf.kinases <- c("AKT1", "BRAF", "MAP2K1", "MAPK1", "MTOR", "PDPK1", "PRKACA", "RPS6KB1", "SRC")

kinase.conditions <- read.delim("data/kinase-condition-pairs.tsv", as.is=TRUE)

kinase.invivoconditions <- read.delim("data/kinase_invivoconditions.csv", as.is=TRUE, sep=",")
names(kinase.invivoconditions) <- c("Condition", "Kinase", "Description", "invivoAct", "pAct")

kinase.cond <- subset(kinase.conditions, Kinase %in% egf.kinases)
## kinase.invivocond <- subset(kinase.invivoconditions, Kinase %in% egf.kinases)
egf.kin.act <- kin.act.filt[rownames(kin.act.filt) %in% egf.kinases, ]
egf.kin.act.pert <- make.max.table(egf.kin.act[,unique(kinase.cond$Condition)], na.threshold)
print.table.info(egf.kin.act.pert)
write.activity.table(egf.kin.act.pert, "egf-kinase-activity-perturbed.tsv")
egf.kin.act.pert.imp <- impute.on.table(egf.kin.act.pert)
egf.kin.act.pert.imp.df <- melt(egf.kin.act.pert.imp)
names(egf.kin.act.pert.imp.df) <- c("kinase", "condition", "value")
kinase.cond.pairs <- expand.grid(rownames(egf.kin.act.pert.imp), colnames(egf.kin.act.pert.imp)[1:(ncol(egf.kin.act.pert.imp)-1)])
egf.kin.act.pert.imp.df$is.pert <- unlist(apply(kinase.cond.pairs, 1,
                                                function(pair){
                                                    cond.kinases <- subset(kinase.cond, Condition==pair[2])$Kinase
                                                    return (pair[1] %in% cond.kinases)
                                                }))
egf.kin.act.pert.imp.df <- egf.kin.act.pert.imp.df[order(egf.kin.act.pert.imp.df$kinase,
                                                         egf.kin.act.pert.imp.df$condition),]
write.table(egf.kin.act.pert.imp.df, "egf-kinase-activity-perturbed-imp.tsv", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

non.pert.conds <- setdiff(colnames(egf.kin.act), kinase.cond$Condition)
egf.kin.act.nonpert <- make.max.table(subset(egf.kin.act, select=non.pert.conds), na.threshold)
print.table.info(egf.kin.act.nonpert)
write.activity.table(egf.kin.act.nonpert, "egf-kinase-activity-unperturbed.tsv")
egf.kin.act.nonpert.imp <- impute.on.table(egf.kin.act.nonpert)
egf.kin.act.nonpert.imp.df <- melt(egf.kin.act.nonpert.imp)
names(egf.kin.act.nonpert.imp.df) <- c("kinase", "condition", "value")
kinase.cond.pairs <- expand.grid(rownames(egf.kin.act.nonpert.imp), colnames(egf.kin.act.nonpert.imp)[1:(ncol(egf.kin.act.nonpert.imp)-1)])
egf.kin.act.nonpert.imp.df$is.pert <- unlist(apply(kinase.cond.pairs, 1,
                                                   function(pair){
                                                       cond.kinases <- subset(kinase.cond, Condition==pair[2])$Kinase
                                                       return (pair[1] %in% cond.kinases)
                                                   }))
egf.kin.act.nonpert.imp.df <- egf.kin.act.nonpert.imp.df[order(egf.kin.act.nonpert.imp.df$kinase,
                                                               egf.kin.act.nonpert.imp.df$condition),]
write.table(egf.kin.act.nonpert.imp.df, "egf-kinase-activity-unperturbed-imp.tsv", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

kin.act.egf.nonpert <- make.max.table(subset(kin.act, select=non.pert.conds), na.threshold)
print.table.info(kin.act.egf.nonpert)
write.activity.table(kin.act.egf.nonpert, "kinase-activity-egf-unperturbed.tsv")
kin.act.egf.nonpert.imp <- impute.on.table(kin.act.egf.nonpert)
kin.act.egf.nonpert.imp.df <- melt(kin.act.egf.nonpert.imp)
names(kin.act.egf.nonpert.imp.df) <- c("kinase", "condition", "value")
kinase.cond.pairs <- expand.grid(rownames(kin.act.egf.nonpert.imp), colnames(kin.act.egf.nonpert.imp)[1:(ncol(kin.act.egf.nonpert.imp)-1)])
kin.act.egf.nonpert.imp.df$is.pert <- unlist(apply(kinase.cond.pairs, 1,
                                                   function(pair){
                                                       cond.kinases <- subset(kinase.cond, Condition==pair[2])$Kinase
                                                       return (pair[1] %in% cond.kinases)
                                                   }))
kin.act.egf.nonpert.imp.df <- kin.act.egf.nonpert.imp.df[order(kin.act.egf.nonpert.imp.df$kinase,
                                                               kin.act.egf.nonpert.imp.df$condition),]
write.table(kin.act.egf.nonpert.imp.df, "kinase-activity-egf-unperturbed-imp.tsv", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
