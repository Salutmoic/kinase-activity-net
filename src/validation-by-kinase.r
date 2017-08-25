suppressMessages(library(ROCR))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(gridExtra))
source("src/rocr-ggplot2.r")

argv <- commandArgs(TRUE)
if (!(length(argv) %in% c(5, 6))){
    stop("USAGE: <script> DATA_FILE VAL_SET NEG_VAL_SET RAND_NEGS DIRECTED [PRED_COLUMN]")
}

pred_score_file <- argv[1]
val_set_file <- argv[2]
neg_val_set_file <- argv[3]
use_rand_negs <- as.logical(argv[4])
directed <- as.logical(argv[5])
if (length(argv) == 6){
    pred_col <- as.integer(argv[6])
}else{
    pred_col <- NULL
}

if (use_rand_negs && directed)
    stop("Select either RAND_NEGS or DIRECTED")

pred_score_full <- read.delim(pred_score_file, as.is = TRUE)
colnames(pred_score_full)[1:2] <- c("prot1", "prot2")
pred_score_full <- pred_score_full[order(pred_score_full$prot1,
                                         pred_score_full$prot2),]
pred_score_full <- subset(pred_score_full, prot1 != prot2)
rownames(pred_score_full) <- paste(pred_score_full$prot1, pred_score_full$prot2,
                                   sep="-")

for (i in 3:ncol(pred_score_full)){
    column <- pred_score_full[,i]
    if (all(!is.numeric(column)))
        next
    col.min <- min(column, na.rm=TRUE)
    col.max <- max(column, na.rm=TRUE)
    col.norm <- (column - col.min)/(col.max - col.min)
    pred_score_full[,i] <- col.norm
    ## We assume that highly specific expression should give lower
    ## chance of a regulatory relationship, so to avoid an inverted
    ## roc/prec-recall curve, convert to expression-generality
    if (colnames(pred_score_full)[i] %in% c("tissue.sel1", "tissue.sel2",
                                            "celline.sel1", "celline.sel2"))
        pred_score_full[,i] <- 1.0 - pred_score_full[,i]
}

true_intxns_tbl <- read.table(val_set_file, as.is = TRUE)

if (is.null(pred_col)){
    pred_cols <- 3:ncol(pred_score_full)
}else{
    pred_cols <- pred_col
}

kinases <- sort(unique(true_intxns_tbl[, 1]))

data_basename <- strsplit(basename(pred_score_file), split = "\\.")[[1]][1]
data_basename <- paste(data_basename, "val-by-kinase", sep="-")
if (use_rand_negs){
    data_basename <- paste(data_basename, "rand-negs", sep="-")
}
if (directed){
    data_basename <- paste(data_basename, "direct", sep="-")
}
out_img <- paste0("img/", data_basename, ".pdf")

pdf(out_img, title=data_basename, width=10.5, height=7)

all.aucs <- NULL
for (kinase in kinases){
    kin_true_intxns_tbl <- subset(true_intxns_tbl, V1==kinase)
    if (nrow(kin_true_intxns_tbl) < 5)
        next
    possible_true_intxns_full <- paste(kin_true_intxns_tbl[, 1],
                                      kin_true_intxns_tbl[, 2],
                                      sep = "-")

    preds <- vector("list", length(pred_cols))
    ## This procedure is sufficient for training-free predictions.
    for (i in pred_cols){
        pred_data <- NULL
        label_data <- NULL
        pred_score_kin <- pred_score_full[which(pred_score_full[,1]==kinase),]
        pred_score <- pred_score_kin[, c(1, 2, i)]
        pred_score <- pred_score[complete.cases(pred_score),]
        if (!all(is.numeric(pred_score[,3]))){
            preds[[i-2]] <- NA
            next
        }
        possible_true_intxns <- intersect(possible_true_intxns_full,
                                          rownames(pred_score))
        if (length(possible_true_intxns)==0){
            preds[[i-2]] <- NA
            next
        }
        if (use_rand_negs){
            possible_false_intxns <- rownames(pred_score)
            n <- 100
        }else if (directed){
            rev_true_intxns <- paste(kin_true_intxns_tbl[, 2],
                                     kin_true_intxns_tbl[, 1], sep = "-")
            possible_false_intxns <- rev_true_intxns
            possible_false_intxns <- intersect(possible_false_intxns,
                                               rownames(pred_score))
            n <- 1
        }else{
            false_intxns_tbl <- read.table(neg_val_set_file, as.is = TRUE)
            possible_false_intxns <- paste(false_intxns_tbl[, 1],
                                           false_intxns_tbl[, 2],
                                           sep = "-")
            possible_false_intxns <- intersect(possible_false_intxns,
                                               rownames(pred_score))
            n <- 1
        }
        possible_false_intxns <- setdiff(possible_false_intxns,
                                         possible_true_intxns)
        sample_size <- round(0.8 * min(length(possible_false_intxns),
                                       length(possible_true_intxns)))
        for (j in 1:n){
            false_intxns <- sample(possible_false_intxns, size=sample_size)
            true_intxns <- sample(possible_true_intxns, size=sample_size)
            ## Put together the tables of predictions and TRUE/FALSE labels
            intxns <- c(true_intxns, false_intxns)
            pred_vals <- pred_score[intxns, 3]
            labels <- c(rep(TRUE, length(true_intxns)),
                        rep(FALSE, length(false_intxns)))
            if (is.null(pred_data)){
                pred_data <- pred_vals
                label_data <- labels
            }else{
                pred_data <- cbind(pred_data, pred_vals)
                label_data <- cbind(label_data, labels)
            }
        }
        pred <- prediction(pred_data, label_data)
        preds[[i-2]] <- pred
    }

    status <- unlist(lapply(preds, function(pred){all(is.na(pred))}))
    if (all(is.na(status))){
        next
    }

    pred_score_melt <- melt(pred_score_kin[,c(1, 2, pred_cols)])[c("prot1", "prot2", "variable", "value")]
    names(pred_score_melt)[3] <- "feature"
    rownames(pred_score_melt) <- rownames(pred_score_kin)
    pred_score_melt$is.known <- rownames(pred_score_melt) %in% possible_true_intxns_full

    feature_names <- list(cor.p.293="CPTAC phospho.",
                          cor.p.272="Wilkes phospho.",
                          prot.atlas.tissue="Prot. Atlas tissue coexpr.",
                          prot.atlas.celline="Prot. Atlas cell-line coexpr.",
                          tissue.coexpr="GTEX coexpr.",
                          max.pssm.score="PSSM",
                          dcg="DCG",
                          max.func.score="Substrate max. func. score",
                          tissue.sel1="Regulator tissue spec.",
                          tissue.sel2="Substrate tissue spec.",
                          celline.sel1="Regulator cell-line spec.",
                          celline.sel2="Substrate cell-line spec.",
                          bart.pred.mean="Mean BART prob. - full",
                          bart.pred.mean.pure="Mean BART prob. - full, pure",
                          assoc.pred.mean="Mean BART prob. - assoc.",
                          direct.pred.mean="Mean BART prob. - direct.")

    pred_score_melt$feature <- factor(pred_score_melt$feature,
                                      levels=levels(pred_score_melt$feature),
                                      labels=feature_names[levels(pred_score_melt$feature)])

    ## Density plot of prediction scores
    if (length(unique(pred_score_melt$feature))>1){
        p1 <- ggplot(pred_score_melt, aes(value, colour=feature, fill=feature)) +
            geom_density(size=1.5, alpha=0.1) +
            xlim(0, 1) +
            theme_bw(base_size=12) +
            labs(title=unique(levels(pred_score_melt$feature))[1])
    }else{
        p1 <- ggplot(pred_score_melt, aes(value, colour=is.known, fill=is.known)) +
            geom_density(size=1.5, alpha=0.1) +
            xlim(0, 1) +
            theme_bw(base_size=12)
    }

    ## ROC curve w/ AUC info
    rocs <- build.perf2d.df(preds, "fpr", "tpr",
                            colnames(pred_score_full)[3:ncol(pred_score_full)],
                            n, TRUE)
    rocs$feature <- factor(rocs$feature,
                           levels=levels(rocs$feature),
                           labels=feature_names[levels(rocs$feature)])
    p2 <- rocr2d.ggplot2(rocs, "false postive rate", "true positive rate", n,
                         c(0, 1), 12)

    aucs <- build.perf1d.df(preds, "auc",
                            colnames(pred_score_full)[3:ncol(pred_score_full)],
                            n, FALSE)
    aucs$feature <- factor(aucs$feature,
                           levels=levels(aucs$feature),
                           labels=feature_names[levels(aucs$feature)])
    if (length(levels(aucs$feature)) > 1){
        abline <- c(0.5, 0)
    }else{
        abline <- NULL
        if (is.null(all.aucs)){
            aucs.tmp <- aucs
            aucs.tmp$kinase <- rep(kinase, nrow(aucs.tmp))
            aucs.tmp$n <- rep(nrow(kin_true_intxns_tbl), nrow(aucs.tmp))
            all.aucs <- aucs.tmp
        }else{
            aucs.tmp <- aucs
            aucs.tmp$kinase <- rep(kinase, nrow(aucs.tmp))
            aucs.tmp$n <- rep(nrow(kin_true_intxns_tbl), nrow(aucs.tmp))
            all.aucs <- rbind(all.aucs, aucs.tmp)
        }
   }
    p3 <- rocr.sum.ggplot2(aucs, "AUC", n, NULL, 12)

    ## Precision-Recall curve
    prs <- build.perf2d.df(preds, "rec", "prec",
                           colnames(pred_score_full)[3:ncol(pred_score_full)],
                           n, TRUE)
    prs$feature <- factor(prs$feature,
                          levels=levels(prs$feature),
                          labels=feature_names[levels(prs$feature)])
    p4 <- rocr2d.ggplot2(prs, "recall", "precision", n, c(0.5, 0), 12)

    grid.arrange(p1, p2, p3, p4,
                 top=paste0(kinase, " (n=", nrow(kin_true_intxns_tbl), ")"))

    ## ## All-vs-all precision-recall
    ## prs <- NULL
    ## for (i in 1:length(preds)){
    ##     if (is.na(preds[[i]]))
    ##         next
    ##     pred_score <- pred_score_full[,c(1, 2, (2+i))]
    ##     pred_score <- pred_score[complete.cases(pred_score),]
    ##     if (!all(is.numeric(pred_score[,3]))){
    ##         next
    ##     }
    ##     possible_true_intxns <- intersect(possible_true_intxns_full,
    ##                                       rownames(pred_score))
    ##     possible_false_intxns <- setdiff(rownames(pred_score), possible_true_intxns)
    ##     intxns <- c(possible_true_intxns, possible_false_intxns)
    ##     pred_score <- pred_score[intxns,]
    ##     labels <- c(rep(TRUE, length(possible_true_intxns)),
    ##                 rep(FALSE, length(possible_false_intxns)))
    ##     pred_full <- prediction(pred_score[,3], labels)
    ##     perf_pr <- performance(pred_full, measure = "prec", x.measure = "rec")
    ##     if (is.null(prs)){
    ##         prs <- data.frame(x=perf_pr@x.values[[1]],
    ##                           y=perf_pr@y.values[[1]],
    ##                           feature=rep(colnames(pred_score_full)[i+2],
    ##                                       length(perf_pr@x.values[[1]])))
    ##     }else{
    ##         new_pr <- data.frame(x=perf_pr@x.values[[1]],
    ##                              y=perf_pr@y.values[[1]],
    ##                              feature=rep(colnames(pred_score_full)[i+2],
    ##                                          length(perf_pr@x.values[[1]])))
    ##         prs <- rbind(prs, new_pr)
    ##     }
    ## }

    ## prs$feature <- factor(prs$feature,
    ##                       levels=levels(prs$feature),
    ##                       labels=feature_names[levels(prs$feature)])

    ## rocr2d.ggplot2(prs, "recall", "precision", 1, NULL)
}

## median.aucs <- all.aucs %>%
##     group_by(kinase) %>%
##     summarise(median.auc=median(y),
##               n=unique(n))

## all.aucs$kinase <- factor(all.aucs$kinase,
##                           levels=median.aucs$kinase[order(median.aucs$median.auc, decreasing=TRUE)])
## median.aucs$kinase <- factor(median.aucs$kinase,
##                              levels=median.aucs$kinase[order(median.aucs$median.auc, decreasing=TRUE)])

if (!is.null(all.aucs)){
    ggplot(data=all.aucs, aes(reorder(kinase, y, FUN=median), y=y)) +
        geom_boxplot(aes(fill=n)) +
        scale_fill_gradient(low="white", high="#045a8d") +
        geom_abline(intercept=0.5, slope=0.0, size=1,
                             linetype=2, colour="grey") +
        ylim(0, 1) +
        ylab("AUC") +
        xlab("kinase") +
        theme_bw(base_size=12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

dev.off()
