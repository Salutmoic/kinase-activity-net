suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx256g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))

argv <- commandArgs(TRUE)
if (length(argv) != 6){
    stop("USAGE: <script> MERGED_PREDS VAL_SET NEG_VAL_SET RAND_NEGS DIRECTED N")
}

merged_pred_file <- argv[1]
val_set_file <- argv[2]
neg_val_set_file <- argv[3]
use_rand_negs <- as.logical(argv[4])
directed <- as.logical(argv[5])
n <- as.integer(sub("^\\\\", "", argv[6]))

set.seed(n)

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
true_intxns <- intersect(possible_true_intxns, rownames(merged_pred))

## Generate some reverse true interactions to be mixed in as negatives
if (use_rand_negs){
    possible_false_intxns <- rownames(merged_pred)
    ## false_pairs <- expand.grid(unique(true_intxns_tbl[,1]),
    ##                            unique(true_intxns_tbl[,2]),
    ##                            stringsAsFactors=FALSE)
    ## false_pairs <- false_pairs[which(false_pairs[,1]!=false_pairs[,2]),]
    ## possible_false_intxns <- paste(false_pairs[, 1], false_pairs[, 2],
    ##                               sep = "-")
}else if (directed){
    rev_true_intxns <- paste(true_intxns_tbl[, 2], true_intxns_tbl[, 1], sep = "-")
    possible_false_intxns <- rev_true_intxns
}else{
    false_intxns_tbl <- read.table(neg_val_set_file, as.is = TRUE)
    possible_false_intxns <- paste(false_intxns_tbl[, 1], false_intxns_tbl[, 2],
                                   sep = "-")
}

possible_false_intxns <- setdiff(possible_false_intxns, true_intxns)
possible_false_intxns <- intersect(possible_false_intxns, rownames(merged_pred))

num_negs <- min(length(possible_false_intxns),
                8*length(true_intxns))

file_base <- strsplit(merged_pred_file, split="\\.")[[1]][1]

if (!dir.exists(file_base)){
    dir.create(file_base)
}

false_intxns <- sample(possible_false_intxns, num_negs)
if (num_negs <= length(true_intxns)){
    true_intxns <- sample(true_intxns, num_negs)
}else{
    true_intxns <- sample(true_intxns)
}

## Put together the tables of predictions and TRUE/FALSE labels
intxns <- c(true_intxns, false_intxns)
preds <- merged_pred[intxns, 3:ncol(merged_pred)]
labels <- c(rep(TRUE, length(true_intxns)),
            rep(FALSE, length(false_intxns)))
bart <- bartMachineCV(preds,
                      as.factor(labels),
                      use_missing_data=(!directed),
                      use_missing_data_dummies_as_covars=(!directed),
                      replace_missing_data_with_x_j_bar=FALSE,
                      impute_missingness_with_x_j_bar_for_lm=FALSE,
                      prob_rule_class=0.5,
                      verbose=FALSE,
                      serialize=TRUE)
save(bart, file=paste0(file_base, "/", sub("out/", "", file_base), "-bart-", n, ".Rdata"))
## check_bart_error_assumptions(bart)
img_file_base <- basename(file_base)
img_file <- paste0("img/", img_file_base, "-info.pdf")
pdf(img_file)
investigate_var_importance(bart, type="splits")
investigate_var_importance(bart, type="trees")
cov_importance_test(bart, plot=TRUE)
dev.off()
