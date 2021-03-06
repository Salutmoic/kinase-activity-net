suppressMessages(library(reshape2))
suppressMessages(library(VIM))
suppressMessages(library(entropy))

filter.act.preds <- function(kin.act, phospho.anno, phospho.vals, kin.sub.tbl,
                             min.sites){
    for (kinase in rownames(kin.act)){
        kin.subs <- subset(kin.sub.tbl, GENE==kinase)
        kin.subs.ensp <- merge(phospho.anno, kin.subs,
                               by.x=c("gene_name", "positions"),
                               by.y=c("SUB_GENE", "position"))
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
        bad.conds <- intersect(bad.conds, colnames(kin.act))
        kin.act[kinase, bad.conds] <- NA
        ## if (length(bad.conds) > 0){
        ##     message(paste(length(bad.conds), "conditions omitted for", kinase,
        ##                   "due to lack of evidence."))
        ## }
    }
    return(kin.act)
}

choose.redundant.kin <- function(row, kin.act){
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

get.redundant.kins <- function(kin.overlap, kin.act){
    kin.overlap <- kin.overlap[order(kin.overlap$proportion.of.substrates.shared,
                                     decreasing=TRUE),]
    redundant.kins <- c()
    for (i in 1:nrow(kin.overlap)){
        row <- unlist(kin.overlap[i,])
        if (row[1] %in% redundant.kins || row[3] %in% redundant.kins)
            next
        if (row[1] %in% rownames(kin.act) && row[3] %in% rownames(kin.act)){
            redundant.kin <- choose.redundant.kin(row, kin.act)
            redundant.kins <- c(redundant.kins, redundant.kin)
        }else{
            if (!(row[1] %in% rownames(kin.act)))
                redundant.kins <- c(redundant.kins, row[1])
            if (!(row[3] %in% rownames(kin.act)))
                redundant.kins <- c(redundant.kins, row[3])
        }
    }
    return(redundant.kins)
}

get.bad.inhib.conds <- function(kin.act){
    kinase.conditions <- read.delim("data/external/kinase-condition-pairs.tsv",
                                    as.is=TRUE)
    conds <- colnames(kin.act)
    bad.inhib.conds <- apply(kinase.conditions, 1,
                             function(row){
                                 cond <- row[1]
                                 kin <- row[6]
                                 direction <- row[7]
                                 if (!(kin %in% rownames(kin.act)) ||
                                     !(cond %in% colnames(kin.act)) ||
                                     is.na(kin.act[kin, cond]) ||
                                     (direction == "down" && kin.act[kin, cond] < 0) ||
                                     (direction == "up" && kin.act[kin, cond] > 0) )
                                     return(FALSE)
                                 return(TRUE)
                             })
    return(unique(kinase.conditions$Condition[which(bad.inhib.conds)]))
}

percent.na <- function(x){
    if (is.null(x) || is.na(length(x)))
        return(0.0)
    return(length(which(is.na(x)))/length(x))
}

make.balanced.table <- function(full.tbl, na.threshold, precious.rows=c()){
    ## Iteratively remove the worst row or column until all rows and
    ## cols have less than the threshold %NA.
    trashable <- setdiff(rownames(full.tbl), precious.rows)
    while (nrow(full.tbl) > 1 && ncol(full.tbl) > 1){
        max.row.na <- max(apply(full.tbl, 1, percent.na))
        max.col.na <- max(apply(full.tbl, 2, percent.na))
        if (max.row.na > max.col.na){
            max.na <- which.max(apply(full.tbl[trashable,], 1, percent.na))
            worst.row <- rownames(full.tbl[trashable,])[max.na]
            new.rows <- setdiff(rownames(full.tbl), worst.row)
            trashable <-  setdiff(trashable, worst.row)
            full.tbl <- full.tbl[new.rows,]
        }else{
            max.na <- which.max(apply(full.tbl, 2, percent.na))
            worst.col <- colnames(full.tbl)[max.na]
            new.cols <- setdiff(colnames(full.tbl), worst.col)
            full.tbl <- subset(full.tbl, select=new.cols)
        }
        if (length(trashable) <= 1 &&
            any(apply(full.tbl, 1, function(tbl.row) percent.na(tbl.row) > na.threshold)))
            return(NULL)
        if (length(full.tbl) == 0 || ncol(full.tbl) == 0 || nrow(full.tbl) == 0)
            return(NULL)
        if (all(apply(full.tbl[trashable,], 1, function(tbl.row) percent.na(tbl.row) < na.threshold)) &&
            all(apply(full.tbl, 2, function(tbl.col) percent.na(tbl.col) < na.threshold))){
            return(full.tbl)
        }
    }
    return(NULL)
}

make.max.row.table <- function(full.tbl, na.threshold, min.cols, precious.rows=c()){
    ## 1) Get rid of all columns that have more than the threshold %NA
    ## 2) Try removing more columns until all of the rows have less than the threshold
    ## 3) If all the rows pass and we haven't removed too many columns, return the table
    ## 4) Otherwise, remove the worst row and repeat from #2
    trashable <- setdiff(rownames(full.tbl), precious.rows)
    bad.cols <- colnames(full.tbl)[which(apply(full.tbl, 2,
                                               function(tbl.col){
                                                   percent.na(tbl.col) > na.threshold
                                               }))]
    good.cols <- setdiff(colnames(full.tbl), bad.cols)
    sub.tbl <- subset(full.tbl, select=good.cols)
    while (!is.null(dim(sub.tbl)) && nrow(sub.tbl) > 1){
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
        max.na.row <- which.max(apply(sub.tbl[trashable,], 1, percent.na))
        good.rows <- setdiff(rownames(sub.tbl), rownames(sub.tbl[trashable,])[max.na.row])
        if (is.null(good.rows) || length(good.rows) == 0)
            return (NULL)
        trashable <- setdiff(trashable, rownames(sub.tbl[trashable,])[max.na.row])
        sub.tbl <- sub.tbl[good.rows,]
    }
    return (NULL)
}

make.max.col.table <- function(full.tbl, na.threshold, min.rows, precious.rows=c()){
    ## 1) Get rid of all rows that have more than the threshold %NA
    ## 2) Try removing more rows until all of the columns have less than the threshold
    ## 3) If all the columns pass and we haven't removed too many rows, return the table
    ## 4) Otherwise, remove the worst column and repeat from #2
    trashable <- setdiff(rownames(full.tbl), precious.rows)
    bad.rows <- rownames(full.tbl[trashable,])[which(apply(full.tbl[trashable,], 1,
                                                           function(tbl.row){
                                                               percent.na(tbl.row) > na.threshold
                                                           }))]
    good.rows <- setdiff(rownames(full.tbl), bad.rows)
    trashable <- setdiff(trashable, bad.rows)
    sub.tbl <- full.tbl[good.rows,]

    while (!is.null(dim(sub.tbl)) && ncol(sub.tbl) > 1){
        tmp.tbl <- sub.tbl
        fail <- FALSE
        while (length(trashable) > 1 &&
               any(apply(tmp.tbl[trashable,], 2, percent.na) > na.threshold)){
            max.na.row <- which.max(apply(tmp.tbl[trashable,], 1, percent.na))
            good.rows <- setdiff(rownames(tmp.tbl), rownames(tmp.tbl[trashable,])[max.na.row])
            if (is.null(good.rows) || length(good.rows) == 0){
                fail <- TRUE
                break
            }
            trashable <- setdiff(trashable, rownames(tmp.tbl[trashable,])[max.na.row])
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
