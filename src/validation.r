suppressMessages(library(ROCR))
source("src/rocr-helpers.r")

argv <- commandArgs(TRUE)
if (length(argv) != 5){
    stop("USAGE: <script> DATA_FILE VAL_SET NEG_VAL_SET RAND_NEGS DIRECTED")
}

pred_score_file <- argv[1]
val_set_file <- argv[2]
neg_val_set_file <- argv[3]
use_rand_negs <- as.logical(argv[4])
directed <- as.logical(argv[5])

if (use_rand_negs && directed)
    stop("Select either RAND_NEGS or DIRECTED")

## Get information from the prediction file name.  Build up a name of
## the prediction method for adding a title to plots, and figure out
## what kind of method was in use in order to determine any kind of
## normalisation that should be done before making predictions.
pred_name <- strsplit(basename(pred_score_file), split = "\\.")[[1]][1]
pred_name <- sub("kinact-", "", pred_name)
assoc_method <- "[undefined]"
method <- "other"
if (grepl("pcor2", pred_score_file)){
    assoc_method <- "Pearson's Correlation^2"
    method <- "cor"
}else if (grepl("pcor", pred_score_file)){
    assoc_method <- "Pearson's Correlation"
    method <- "cor"
}else if (grepl("scor2", pred_score_file)){
    assoc_method <- "Spearman's Correlation^2"
    method <- "cor"
}else if (grepl("scor", pred_score_file)){
    assoc_method <- "Spearman's Correlation"
    method <- "cor"
}else if (grepl("paircor", pred_score_file)){
    assoc_method <- "Pairwise Spearman's Correlation"
    method <- "cor"
}else if (grepl("nfchisq", pred_score_file)){
    assoc_method <- "FunChiSq"
    method <- "nfchisq"
}else if (grepl("fnn_mut_info", pred_score_file)){
    assoc_method <- "FNN Mutual Information"
    method <- "mutinfo"
}else if (grepl("mut_info", pred_score_file)){
    assoc_method <- "Mutual Information"
    method <- "mutinfo"
}else if (grepl("partcor", pred_score_file)){
    assoc_method <- "Partial Correlation"
    method <- "cor"
}else if (grepl("fvalue", pred_score_file)){
    assoc_method <- "ANOVA F value_"
    method <- "meansq"
}else if (grepl("pssm", pred_score_file)){
    assoc_method <- "PSSM Score"
    method <- "pssm"
}else if (grepl("string-coexp2", pred_score_file)){
    assoc_method <- "STRING Coexpression Score (cor^2)"
}else if (grepl("string-coexp", pred_score_file)){
    assoc_method <- "STRING Coexpression Score"
}else if (grepl("string-exper", pred_score_file)){
    assoc_method <- "STRING Experimental Score"
}else if (grepl("all", pred_score_file)){
    assoc_method <- "Mean Score"
}else{
    assoc_method <- pred_score_file
}
## The next two append info to the plot title and override for
## normalisation purposes the method used, in order of increasing
## precedence.
if (grepl("filter", pred_score_file)){
    assoc_method <- paste0(assoc_method, " (filtered)")
    method <- "filtered"
}
if (grepl("final-predictor", pred_score_file)){
    assoc_method <- paste0(assoc_method, " (merged)")
    method <- "predictor"
}
if (grepl("balanced", pred_score_file)){
    table_method <- "balanced table"
}else if (grepl("max-rows", pred_score_file)){
    table_method <- "max kinases"
}else if (grepl("max-cols", pred_score_file)){
    table_method <- "max conditions"
}else{
    table_method <- ""
}
pred_score <- read.delim(pred_score_file, as.is = TRUE)
names(pred_score) <- c("prot1", "prot2", "pred_score")
pred_score <- subset(pred_score, prot1 != prot2)
pred_score <- pred_score[complete.cases(pred_score),]
rownames(pred_score) <- paste(pred_score$prot1, pred_score$prot2, sep = "-")

if (method %in% c("cor", "nfchisq")){
    min_pred <- min(pred_score$pred_score)
    max_pred <- max(pred_score$pred_score)
    pred_score$pred_score <- ((pred_score$pred_score - min_pred) /
                              (max_pred - min_pred))
}else if (method %in% c("mutinfo", "fvalue")){
    pred_score$pred_score <- pred_score$pred_score / max(pred_score$pred_score)
}


true_intxns_tbl <- read.table(val_set_file, as.is = TRUE)
possible_true_intxns <- paste(true_intxns_tbl[, 1], true_intxns_tbl[, 2],
                              sep = "-")
possible_true_intxns <- intersect(possible_true_intxns, rownames(pred_score))

if (use_rand_negs){
    possible_false_intxns <- rownames(pred_score)
    ## foo <- expand.grid(true_intxns_tbl[,1], true_intxns_tbl[,2], stringsAsFactors=FALSE)
    ## possible_false_intxns <- paste(foo[,1], foo[,2], sep="-")
    ## possible_false_intxns <- intersect(possible_false_intxns, rownames(pred_score))
}else if (directed){
    rev_true_intxns <- paste(true_intxns_tbl[, 2], true_intxns_tbl[, 1], sep = "-")
    possible_false_intxns <- rev_true_intxns
}else{
    false_intxns_tbl <- read.table(neg_val_set_file, as.is = TRUE)
    possible_false_intxns <- paste(false_intxns_tbl[, 1], false_intxns_tbl[, 2],
                                   sep = "-")
    possible_false_intxns <- intersect(possible_false_intxns, rownames(pred_score))
}

possible_false_intxns <- setdiff(possible_false_intxns, possible_true_intxns)

pred_data <- NULL
label_data <- NULL

print(c(length(possible_false_intxns), length(possible_true_intxns)))

n <- 100

sample_size <- round(0.8 * min(length(possible_false_intxns),
                               length(possible_true_intxns)))

## This procedure is sufficient for training-free predictions.
for (i in 1:n){
    ## false_intxns <- sample(possible_false_intxns, size=sample_size)
    false_intxns <- sample(possible_false_intxns, size=sample_size)
    true_intxns <- sample(possible_true_intxns, size=sample_size)
    ## Put together the tables of predictions and TRUE/FALSE labels
    intxns <- c(true_intxns, false_intxns)
    preds <- pred_score[intxns, "pred_score"]
    labels <- c(rep(TRUE, length(true_intxns)), rep(FALSE, length(false_intxns)))
    if (is.null(pred_data)){
        pred_data <- preds
        label_data <- labels
    }else{
        pred_data <- cbind(pred_data, preds)
        label_data <- cbind(label_data, labels)
    }
}

pred <- prediction(pred_data, label_data)

data_basename <- strsplit(basename(pred_score_file), split = "\\.")[[1]][1]
data_basename <- paste(data_basename, "val", sep="-")
if (use_rand_negs){
    data_basename <- paste(data_basename, "rand-negs", sep="-")
}
if (directed){
    data_basename <- paste(data_basename, "direct", sep="-")
}
out_img <- paste0("img/", data_basename, ".pdf")

pdf(out_img)
par(cex = 1.25, cex.main = 0.8)

## ROC curve w/ AUC info
alpha <- 0.9

perf_roc <- performance(pred, measure = "tpr", x.measure = "fpr")
perf_auc <- performance(pred, "auc")
mean_auc <- mean(unlist(perf_auc@y.values))
se_auc <- sd(unlist(perf_auc@y.values)) / sqrt(n)
perf_f <- performance(pred, "f", alpha=alpha)
min_prec <- 0.5
min_prec_f_val <- 1/(alpha*(1/min_prec)+(1-alpha))
## perf_f@y.values <- lapply(perf_f@y.values,
##                                function(vals) vals/min_prec_f_val - 1.0)
avg_fs <- perf_vert_avg(perf_f)
max_f <- max(avg_fs@y.values[[1]])
recs <- seq(0.05, 1.0, length.out=100)
precs <- (alpha*min_prec_f_val)/(1.0-((1.0-alpha)*min_prec_f_val)/recs)
val_precs <- which(precs <= 1.0 & precs >= min_prec)
precs <- precs[val_precs]
recs <- recs[val_precs]
message(assoc_method)
message(paste("Mean AUC =", format(mean_auc, digits = 2)))
message(paste("S.E.M. =", format(se_auc, digits = 2)))
message(paste("Max F-score =", format(max_f, digits = 2)))
message(paste("n =", n))
message(paste("sample size =", sample_size))


possible_false_intxns <- setdiff(rownames(pred_score), possible_true_intxns)
intxns <- c(possible_true_intxns, possible_false_intxns)
preds <- pred_score[intxns, "pred_score"]
labels <- c(rep(TRUE, length(possible_true_intxns)),
            rep(FALSE, length(possible_false_intxns)))
pred_full <- prediction(preds, labels)
perf_f_full <- performance(pred_full, "f", alpha=alpha)
perf_prec_full <- performance(pred_full, "prec", x.measure="rec")
min_prec <- length(possible_true_intxns)/(length(possible_true_intxns)+length(possible_false_intxns))
min_prec_f_val <- 1/(alpha*(1/min_prec)+(1-alpha))
print(c(min_prec, min_prec_f_val))
## perf_f_full@y.values[[1]] <- perf_f_full@y.values[[1]]/min_prec_f_val - 1.0
max_f_full <- max(perf_f_full@y.values[[1]], na.rm=TRUE)
max_f_full_cutoff <- perf_f_full@x.values[[1]][which.max(perf_f_full@y.values[[1]])]
recs_full <- seq(0.05, 1.0, length.out=100)
precs_full <- (alpha*min_prec_f_val)/(1.0-((1.0-alpha)*min_prec_f_val)/recs)
val_precs_full <- which(precs_full <= 1.0 & precs_full >= min_prec)
precs_full <- precs_full[val_precs_full]
recs_full <- recs_full[val_precs_full]
message(paste("Full predictor max F-score =",
              format(max_f_full, digits = 2),
              "at",
              format(max_f_full_cutoff, digits = 2)))
plot(density(preds))
plot(perf_roc, avg = "vertical", spread.estimate = "boxplot",
     main = paste(assoc_method,
                  table_method,
                  paste0("Mean AUC=", format(mean_auc, digits = 2)),
                  paste0("S.E.M.=", format(se_auc, digits = 2)),
                  paste0("Sample size=", sample_size),
                  sep = "\n"))
## Precision-recall curve
perf_pr <- performance(pred, measure="prec", x.measure="rec")
plot(perf_pr, avg="threshold", spread.estimate="stderror", colorize=TRUE,
     ylim=c(0.0, 1.0), lwd=2)
lines(recs, precs, lty=2, col="grey", lwd=1.5)
plot(perf_f, avg="vertical", spread.estimate="boxplot", lwd=2)
plot(perf_prec_full, avg="threshold", spread.estimate="stderror", colorize=TRUE,
     lwd=2, ylim=c(0.0, 1.0))
lines(recs_full, precs_full, lty=2, col="grey", lwd=2)
plot(perf_f_full, lwd=2)
dev.off()
