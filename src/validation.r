suppressMessages(library(ROCR))
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
source("src/rocr-ggplot2.r")
library(cowplot)

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
possible_true_intxns_full <- paste(true_intxns_tbl[, 1], true_intxns_tbl[, 2],
                              sep = "-")

if (is.null(pred_col)){
    pred_cols <- 3:ncol(pred_score_full)
}else{
    pred_cols <- pred_col
}
preds <- vector("list", length(pred_cols))
## This procedure is sufficient for training-free predictions.
for (i in pred_cols){
    pred_data <- NULL
    label_data <- NULL
    pred_score <- pred_score_full[,c(1, 2, i)]
    pred_score <- pred_score[complete.cases(pred_score),]
    if (!all(is.numeric(pred_score[,3]))){
        preds[[i-2]] <- NA
        next
    }
    possible_true_intxns <- intersect(possible_true_intxns_full,
                                      rownames(pred_score))
    if (use_rand_negs){
        possible_false_intxns <- rownames(pred_score)
        n <- 100
    }else if (directed){
        rev_true_intxns <- paste(true_intxns_tbl[, 2], true_intxns_tbl[, 1], sep = "-")
        possible_false_intxns <- rev_true_intxns
        possible_false_intxns <- intersect(possible_false_intxns, rownames(pred_score))
        n <- 1
    }else{
        false_intxns_tbl <- read.table(neg_val_set_file, as.is = TRUE)
        possible_false_intxns <- paste(false_intxns_tbl[, 1], false_intxns_tbl[, 2],
                                       sep = "-")
        possible_false_intxns <- intersect(possible_false_intxns, rownames(pred_score))
        n <- 1
    }
    possible_false_intxns <- setdiff(possible_false_intxns, possible_true_intxns)
    ## sample_size <- round(0.8 * min(length(possible_false_intxns),
    ##                                length(possible_true_intxns)))
    num_negs <- min(length(possible_false_intxns),
                    8*(2/3)*length(possible_true_intxns))
    for (j in 1:n){
        ## false_intxns <- sample(possible_false_intxns, size=sample_size)
        ## true_intxns <- sample(possible_true_intxns, size=sample_size)
        false_intxns <- sample(possible_false_intxns, num_negs)
        if (num_negs <= length(possible_true_intxns)){
            true_intxns <- sample(possible_true_intxns, num_negs)
        }else{
            true_intxns <- sample(possible_true_intxns, (2/3)*length(possible_true_intxns))
        }
        ## Put together the tables of predictions and TRUE/FALSE labels
        intxns <- c(true_intxns, false_intxns)
        pred_vals <- pred_score[intxns, 3]
        labels <- c(rep(TRUE, length(true_intxns)), rep(FALSE, length(false_intxns)))
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

if (all(is.na(preds))){
    stop("No numeric predictors present")
}

data_basename <- strsplit(basename(pred_score_file), split = "\\.")[[1]][1]
data_basename <- paste(data_basename, "val", sep="-")
if (use_rand_negs){
    data_basename <- paste(data_basename, "rand-negs", sep="-")
}
if (directed){
    data_basename <- paste(data_basename, "direct", sep="-")
}
out_img <- paste0("img/", data_basename, ".pdf")
rds_basename <- paste0("img/rds/", data_basename)

pred_score_melt <- melt(pred_score_full[,c(1, 2, pred_cols)])[c("prot1", "prot2", "variable", "value")]
names(pred_score_melt)[3] <- "feature"

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
                      bart.pred.mean.hinegs="Mean BART prob. - full, hinegs",
                      assoc.pred.mean="Mean BART prob. - assoc.",
                      direct.pred.mean="Mean BART prob. - direct.")

pred_score_melt$feature <- factor(pred_score_melt$feature,
                                  levels=levels(pred_score_melt$feature),
                                  labels=feature_names[levels(pred_score_melt$feature)])

if (length(unique(levels(pred_score_melt$feature)))==1){
    fig_width=7
}else{
    fig_width=10.5
}

pdf(out_img, title=data_basename, width=fig_width, height=7)

## Density plot of prediction scores
dens_plot <- ggplot(pred_score_melt, aes(value, colour=feature, fill=feature)) +
    geom_density(alpha=0.1) +
    xlim(0, 1)
if (length(unique(levels(pred_score_melt$feature)))==1)
    dens_plot <- dens_plot + guides(colour=FALSE, fill=FALSE)
saveRDS(dens_plot, file=paste0(rds_basename, "-dens.rds"))
print(dens_plot)

## ROC curve w/ AUC info
rocs <- build.perf2d.df(preds, "fpr", "tpr",
                      colnames(pred_score_full)[3:ncol(pred_score_full)],
                      n, TRUE)
rocs$feature <- factor(rocs$feature,
                       levels=levels(rocs$feature),
                       labels=feature_names[levels(rocs$feature)])
roc_plot <- rocr2d.ggplot2(rocs, "false postive rate", "true positive rate", n, c(0, 1))
saveRDS(roc_plot, file=paste0(rds_basename, "-roc.rds"))
saveRDS(rocs, file=paste0(rds_basename, "-roc-data.rds"))
print(roc_plot)

aucs <- build.perf1d.df(preds, "auc",
                        colnames(pred_score_full)[3:ncol(pred_score_full)],
                        n, FALSE)
aucs$feature <- factor(aucs$feature,
                       levels=levels(aucs$feature),
                       labels=feature_names[levels(aucs$feature)])
auc_plot <- rocr.sum.ggplot2(aucs, "AUC", n, c(0.5, 0))
saveRDS(auc_plot, file=paste0(rds_basename, "-auc.rds"))
saveRDS(aucs, file=paste0(rds_basename, "-auc-data.rds"))
print(auc_plot)

## Precision-Recall curve
prs <- build.perf2d.df(preds, "rec", "prec",
                      colnames(pred_score_full)[3:ncol(pred_score_full)],
                      n, TRUE)

prs$feature <- factor(prs$feature,
                      levels=levels(prs$feature),
                      labels=feature_names[levels(prs$feature)])

prec_plot <- rocr2d.ggplot2(prs, "recall", "precision", n, c(1/9, 0))
saveRDS(prec_plot, file=paste0(rds_basename, "-prec.rds"))
saveRDS(prs, file=paste0(rds_basename, "-prec-data.rds"))
print(prec_plot)

## All-vs-all precision-recall
prs_full <- NULL
for (i in 1:length(preds)){
    if (is.na(preds[[i]]))
        next
    pred_score <- pred_score_full[,c(1, 2, (2+i))]
    pred_score <- pred_score[complete.cases(pred_score),]
    if (!all(is.numeric(pred_score[,3]))){
        next
    }
    possible_true_intxns <- intersect(possible_true_intxns_full,
                                      rownames(pred_score))
    possible_false_intxns <- setdiff(rownames(pred_score), possible_true_intxns)
    intxns <- c(possible_true_intxns, possible_false_intxns)
    pred_score <- pred_score[intxns,]
    labels <- c(rep(TRUE, length(possible_true_intxns)),
                rep(FALSE, length(possible_false_intxns)))
    pred_full <- prediction(pred_score[,3], labels)
    perf_pr <- performance(pred_full, measure = "prec", x.measure = "rec")
    if (is.null(prs_full)){
        prs_full <- data.frame(x=perf_pr@x.values[[1]],
                          y=perf_pr@y.values[[1]],
                          feature=rep(colnames(pred_score_full)[i+2],
                                      length(perf_pr@x.values[[1]])))
    }else{
        new_pr <- data.frame(x=perf_pr@x.values[[1]],
                             y=perf_pr@y.values[[1]],
                             feature=rep(colnames(pred_score_full)[i+2],
                                         length(perf_pr@x.values[[1]])))
        prs_full <- rbind(prs_full, new_pr)
    }
}

prs_full$feature <- factor(prs_full$feature,
                      levels=levels(prs_full$feature),
                      labels=feature_names[levels(prs_full$feature)])

prs_full$y[which(is.nan(prs_full$y))] <- NA
prs_full$y[which(is.infinite(prs_full$y))] <- NA
prs_full$x[which(is.nan(prs_full$x))] <- NA
prs_full$x[which(is.infinite(prs_full$x))] <- NA

prec_full_plot <- rocr2d.ggplot2(prs_full, "recall", "precision", 1, NULL)
saveRDS(prec_full_plot, file=paste0(rds_basename, "-prec-full.rds"))
saveRDS(prs_full, file=paste0(rds_basename, "-prec-full-data.rds"))
print(prec_full_plot)

prs_full$y.ci <- rep(NA, nrow(prs_full))
prs_full$neg.set <- rep("full", nrow(prs_full))
prs$neg.set <- rep("downsampled", nrow(prs))

prs_merged <- rbind(prs, prs_full)
prs_merged$neg.set <- as.factor(prs_merged$neg.set)
prs_merged$y.min <- prs_merged$y - prs_merged$y.ci
prs_merged$y.max <- prs_merged$y + prs_merged$y.ci
errorbar_rows <- c()
for (f in levels(prs_merged$feature)){
    for (ns in levels(prs_merged$neg.set)){
        prs_merged_sub <- which(prs_merged$feature==f & prs_merged$neg.set==ns)
        rows <- prs_merged_sub[seq(1, length(prs_merged_sub), length=6)[2:5]]
        errorbar_rows <- c(errorbar_rows, rows)
    }
}
errorbar_df <- data.frame(prs_merged[errorbar_rows,])

prs_poly <- rbind(prs_full, prs[nrow(prs):1,])[,1:2]

prec_merged_plot <- ggplot(prs_merged, aes(x=x, y=y, colour=neg.set)) +
    geom_polygon(data=prs_poly, aes(x=x, y=y), fill="#E9E9F7", colour="#E9E9F7") +
    geom_line(size=1.5) +
    scale_fill_identity() +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(x="recall", y="precision") +
    geom_errorbar(mapping=aes(x=x, ymin=y.min, ymax=y.max),
                  data=errorbar_df,
                  inherit.aes=FALSE,
                  width=0.025)
saveRDS(prec_merged_plot, file=paste0(rds_basename, "-prec-merged.rds"))
print(prec_merged_plot)

dev.off()
