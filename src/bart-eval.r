suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx512g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))
library(ggplot2)
source("src/rocr-ggplot2.r")

argv <- commandArgs(TRUE)
if (length(argv) != 5){
    stop("USAGE: <script> MERGED_PREDS VAL_SET NEG_VAL_SET RAND_NEGS DIRECTED")
}

merged_pred_file <- argv[1]
val_set_file <- argv[2]
neg_val_set_file <- argv[3]
use_rand_negs <- as.logical(argv[4])
directed <- as.logical(argv[5])

if (use_rand_negs && directed)
    stop("Select either RAND_NEGS or DIRECTED")

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")
## Remove missing data
if (directed){
    merged_pred <- merged_pred[complete.cases(merged_pred),]
}else{
    good_rows <- which(apply(merged_pred[,3:ncol(merged_pred)], 1,
                             function(row){
                                 any(!is.na(row))
                             }))
    merged_pred <- merged_pred[good_rows,]
}
merged_pred <- subset(merged_pred, prot1 != prot2)
rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2, sep="-")

true_intxns_tbl <- read.table(val_set_file, as.is = TRUE)
possible_true_intxns <- paste(true_intxns_tbl[, 1], true_intxns_tbl[, 2],
                              sep = "-")
possible_true_intxns <- intersect(possible_true_intxns, rownames(merged_pred))

## Generate some reverse true interactions to be mixed in as negatives
if (use_rand_negs){
    possible_false_intxns <- rownames(merged_pred)
}else if (directed){
    rev_true_intxns <- paste(true_intxns_tbl[, 2], true_intxns_tbl[, 1], sep = "-")
    possible_false_intxns <- rev_true_intxns
    possible_false_intxns <- intersect(possible_false_intxns, rownames(merged_pred))
}else{
    false_intxns_tbl <- read.table(neg_val_set_file, as.is = TRUE)
    possible_false_intxns <- paste(false_intxns_tbl[, 1], false_intxns_tbl[, 2],
                                   sep = "-")
    possible_false_intxns <- intersect(possible_false_intxns, rownames(merged_pred))
}
possible_false_intxns <- setdiff(possible_false_intxns, possible_true_intxns)

## sample_size <- round(0.8 * min(length(possible_false_intxns),
##                                length(possible_true_intxns)))
num_negs <- min(length(possible_false_intxns),
                8*length(possible_true_intxns))

num_reps <- 20
k <- 3
all_probs <- NULL
all_labels <- NULL

for (n in 1:num_reps){
    message(paste("Rep.", n))
    ## false_intxns <- sample(possible_false_intxns, size=sample_size)
    ## true_intxns <- sample(possible_true_intxns, size=sample_size)
    false_intxns <- sample(possible_false_intxns, num_negs)
    if (num_negs <= length(possible_true_intxns)){
        true_intxns <- sample(possible_true_intxns, num_negs)
    }else{
        true_intxns <- sample(possible_true_intxns)
    }
    ## Put together the tables of predictions and TRUE/FALSE labels
    intxns <- c(true_intxns, false_intxns)
    preds <- merged_pred[intxns, 3:ncol(merged_pred)]
    labels <- c(rep(TRUE, length(true_intxns)),
                rep(FALSE, length(false_intxns)))
    rand.rows <- sample(1:nrow(preds))
    preds <- preds[rand.rows,]
    labels <- as.factor(labels[rand.rows])
    val_n <- ceiling(nrow(preds)/k)
    train_n <- nrow(preds) - val_n
    rep_probs <- NULL
    rep_labels <- NULL
    for (i in 1:k){
        val_i <- (val_n*(i-1))+1
        if (i < k){
            val_rows <- val_i:(val_i + val_n - 1)
        }else{
            val_rows <- val_i:nrow(preds)
        }
        train_rows <- setdiff(1:nrow(preds), val_rows)
        val_set <- preds[val_rows,]
        val_labels <- labels[val_rows]
        train_set <- preds[train_rows,]
        train_labels <- labels[train_rows]
        bart <- bartMachineCV(train_set,
                              train_labels,
                              use_missing_data=TRUE,
                              use_missing_data_dummies_as_covars=(!directed),
                              replace_missing_data_with_x_j_bar=FALSE,
                              impute_missingness_with_x_j_bar_for_lm=FALSE,
                              prob_rule_class=0.5,
                              verbose=FALSE,
                              serialize=FALSE)
        ## val_preds <- bart_predict_for_test_data(bart, val_set, val_labels,
        ##                                         prob_rule_class=0.0)
        val_preds <- 1.0-bart_machine_get_posterior(bart, val_set)$y_hat
        if (is.null(all_probs)){
            all_probs <- val_preds
            all_labels <- val_labels
        }else{
            if (length(val_rows) < val_n){
                num.na <- val_n - length(val_rows)
                val_preds <- c(val_preds, rep(NA, num.na))
                val_labels <- c(val_labels, rep(NA, num.na))
                val_labels <- as.factor(val_labels)
                levels(val_labels) <- levels(labels)
            }
            all_probs <- cbind.data.frame(all_probs, val_preds)
            all_labels <- cbind.data.frame(all_labels, val_labels)
        }
    }
}

rocr_pred <- prediction(all_probs, all_labels)
rocr_roc <- build.cvperf2d.df(rocr_pred, "fpr", "tpr", c("BART prob."),
                              ncol(all_probs), TRUE)
rocr_prec <- build.cvperf2d.df(rocr_pred, "rec", "prec", c("BART prob."),
                              ncol(all_probs), TRUE)
rocr_auc <- build.cvperf.sum.df(rocr_pred, "auc", c("BART prob."),
                                ncol(all_probs), FALSE)

file_base <- strsplit(basename(merged_pred_file), split="\\.")[[1]][1]
img_file <- paste0(file_base, "-bart")
if (use_rand_negs){
    img_file <- paste(img_file, "rand-negs", sep="-")
}
if (directed){
    img_file <- paste(img_file, "direct", sep="-")
}
rds_basename <- paste0("img/rds/", img_file)
img_file <- paste0("img/", img_file, "-roc.pdf")

mean_auc <- mean(rocr_auc$y, na.rm=TRUE)
roc_plot <- rocr2d.ggplot2(rocr_roc, "false positive rate", "true positive rate",
                           ncol(all_probs), c(0, 1))##  +
    ## geom_text(x=1.0, y=0.1, hjust=1.0,
##           label=paste("Mean AUC =", format(mean_auc, digits=2)))
print(paste("Mean AUC =", format(mean_auc, digits=2)))
if (num_negs <= length(possible_true_intxns)){
    line_intcpt <- 0.5
}else{
    num_true <- length(possible_true_intxns)
    line_intcpt <- num_true / (num_true + num_negs)
}
prec_plot <- rocr2d.ggplot2(rocr_prec, "recall", "precision",
                            ncol(all_probs), c(line_intcpt, 0))
auc_plot <- rocr.sum.ggplot2(rocr_auc, "AUROC", ncol(all_probs), c(0.5,  0)) +
    geom_text(x=0.1, y=1, label=paste("Mean AUROC =", format(mean_auc, digits=2)))
saveRDS(roc_plot, file=paste0(rds_basename, "-cv-roc.rds"))
saveRDS(rocr_roc, file=paste0(rds_basename, "-cv-roc-data.rds"))
saveRDS(prec_plot, file=paste0(rds_basename, "-cv-prec.rds"))
saveRDS(rocr_prec, file=paste0(rds_basename, "-cv-prec-data.rds"))
saveRDS(auc_plot, file=paste0(rds_basename, "-cv-auc.rds"))
saveRDS(rocr_auc, file=paste0(rds_basename, "-cv-auc-data.rds"))
pdf(img_file)
print(roc_plot)
print(prec_plot)
print(auc_plot)
dev.off()

rocr_roc <- performance(rocr_pred, measure="tpr", x.measure="fpr")
rocr_prec <- performance(rocr_pred, measure="prec", x.measure="rec")
rocr_auc <- performance(rocr_pred, "auc")
mean_auc <- mean(unlist(rocr_auc@y.values), na.rm=TRUE)
se_auc <- sd(unlist(rocr_auc@y.values), na.rm=TRUE)/sqrt(num_reps*k)

fprs <- rocr_roc@x.values
fpr0.1s <- unlist(lapply(1:length(fprs),
                        function(f){
                            fpr0.1 <-max(which(fprs[[f]]<0.1))
                            rocr_roc@alpha.values[[f]][fpr0.1]
                        }))
fpr.cutoff <- mean(fpr0.1s)
message(paste("FPR 0.1 cutoff: ", fpr.cutoff))
precs <- rocr_prec@y.values
prec0.9s <- unlist(lapply(1:length(precs),
                        function(f){
                            prec0.9 <-max(which(precs[[f]]>0.9))
                            rocr_prec@alpha.values[[f]][prec0.9]
                        }))
prec.cutoff <- mean(prec0.9s)
message(paste("Prec 0.9 cutoff: ", prec.cutoff))
