suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx16g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))

argv <- commandArgs(TRUE)
if (length(argv) != 5){
    stop("USAGE: <script> MERGED_PREDS VAL_SET NEG_VAL_SET RAND_NEGS DIRECTED")
}

merged_pred_file <- argv[1]
val_set_file <- argv[2]
neg_val_set_file <- argv[3]
use_rand_negs <- as.logical(argv[4])
directed <- as.logical(argv[5])

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")
merged_pred <- subset(merged_pred, prot1 != prot2)
## merged_pred <- merged_pred[complete.cases(merged_pred),]
rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2, sep="-")

for (i in 3:ncol(merged_pred)){
    min_pred <- min(merged_pred[,i], na.rm=TRUE)
    max_pred <- max(merged_pred[,i], na.rm=TRUE)
    merged_pred[,i] <- (merged_pred[,i] - min_pred)/(max_pred - min_pred)
}

true_intxns_tbl <- read.table(val_set_file, as.is = TRUE)
possible_true_intxns <- paste(true_intxns_tbl[, 1], true_intxns_tbl[, 2],
                              sep = "-")
true_intxns <- intersect(possible_true_intxns, rownames(merged_pred))

## Generate some reverse true interactions to be mixed in as negatives
if (use_rand_negs){
    rand_neg_intxns <- rownames(merged_pred)
    possible_false_intxns <- setdiff(rand_neg_intxns, possible_true_intxns)
}else{
    false_intxns_tbl <- read.table(neg_val_set_file, as.is = TRUE)
    possible_false_intxns <- paste(false_intxns_tbl[, 1], false_intxns_tbl[, 2],
                                   sep = "-")
    possible_false_intxns <- intersect(possible_false_intxns, rownames(merged_pred))
}

rev_true_intxns <- paste(true_intxns_tbl[, 2], true_intxns_tbl[, 1], sep = "-")
rev_true_intxns <- setdiff(possible_true_intxns, rev_true_intxns)
if (directed){
    ## The predictor is directed, so include reverse-true interactions
    ## as possible false ones since we take bi-directional
    ## relationships to be statistically unlikely
    possible_false_intxns <- union(possible_false_intxns, rev_true_intxns)
}else{
    ## Since we're testing for symmetrical associations, remove
    ## reverse-true interactions from the possible false interactions
    possible_false_intxns <- setdiff(possible_false_intxns, rev_true_intxns)
}

num_negs <- min(length(possible_false_intxns),
                3*length(true_intxns))
false_intxns <- sample(possible_false_intxns, num_negs)
true_intxns <- possible_true_intxns

## Put together the tables of predictions and TRUE/FALSE labels
intxns <- c(true_intxns, false_intxns)
preds <- merged_pred[intxns, 3:ncol(merged_pred)]
labels <- c(rep(TRUE, length(true_intxns)),
            rep(FALSE, length(false_intxns)))
bart <- bartMachineCV(preds,
                      as.integer(labels),
                      use_missing_data=TRUE,
                      replace_missing_data_with_x_j_bar=FALSE,
                      impute_missingness_with_x_j_bar_for_lm=FALSE,
                      verbose=FALSE,
                      serialize=TRUE)

file_base <- strsplit(merged_pred_file, split="\\.")[[1]][1]
if (use_rand_negs){
    file_base <- paste0(file_base, "rand-negs", sep="-")
}
if (directed){
    file_base <- paste0(file_base, "directed", sep="-")
}
save(bart, file=paste0(file_base, "-bart.Rdata"))

file_base <- basename(file_base)
img_file <- paste0("img/", file_base, "-info.pdf")

pdf(img_file)
## check_bart_error_assumptions(bart)
investigate_var_importance(bart, type="splits")
investigate_var_importance(bart, type="trees")
cov_importance_test(bart, plot=TRUE)
dev.off()
