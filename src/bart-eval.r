suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx32g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))

argv <- commandArgs(TRUE)
if (length(argv) != 6){
    stop("USAGE: <script> MERGED_PREDS VAL_SET NEG_VAL_SET RAND_NEGS CLASSIFIER DIRECTED")
}

merged_pred_file <- argv[1]
val_set_file <- argv[2]
neg_val_set_file <- argv[3]
use_rand_negs <- as.logical(argv[4])
classifier <- as.logical(argv[5])
directed <- as.logical(argv[6])

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")
merged_pred <- subset(merged_pred, prot1 != prot2)
merged_pred <- merged_pred[complete.cases(merged_pred),]
rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2, sep="-")

true_intxns_tbl <- read.table(val_set_file, as.is = TRUE)
possible_true_intxns <- paste(true_intxns_tbl[, 1], true_intxns_tbl[, 2],
                              sep = "-")
possible_true_intxns <- intersect(possible_true_intxns, rownames(merged_pred))

## Generate some reverse true interactions to be mixed in as negatives
if (use_rand_negs){
    possible_false_intxns <- rownames(merged_pred)
    possible_false_intxns <- setdiff(possible_false_intxns, possible_true_intxns)
}else{
    false_intxns_tbl <- read.table(neg_val_set_file, as.is = TRUE)
    possible_false_intxns <- paste(false_intxns_tbl[, 1], false_intxns_tbl[, 2],
                                   sep = "-")
    possible_false_intxns <- intersect(possible_false_intxns, rownames(merged_pred))
}
possible_false_intxns <- setdiff(possible_false_intxns, possible_true_intxns)

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

sample_size <- round(0.8 * min(length(possible_false_intxns),
                               length(possible_true_intxns)))

num_reps <- 50
k <- 5
all_probs <- NULL
all_labels <- NULL

for (n in 1:num_reps){
    message(paste("Rep.", n))
    false_intxns <- sample(possible_false_intxns, size=sample_size)
    true_intxns <- sample(possible_true_intxns, size=sample_size)
    ## Put together the tables of predictions and TRUE/FALSE labels
    intxns <- c(true_intxns, false_intxns)
    preds <- merged_pred[intxns, 3:ncol(merged_pred)]
    labels <- c(rep(TRUE, length(true_intxns)),
                rep(FALSE, length(false_intxns)))
    rand.rows <- sample(1:nrow(preds))
    preds <- preds[rand.rows,]
    labels <- labels[rand.rows]
    val_n <- floor(nrow(preds)/k)
    train_n <- nrow(preds) - val_n
    rep_probs <- NULL
    rep_labels <- NULL
    for (i in 1:k){
        val_i <- ((i-1)*k + 1)
        if (i < k){
            val_rows <- val_i:(val_i + val_n)
        }else{
            val_rows <- val_i:nrow(preds)
        }
        train_rows <- setdiff(1:nrow(preds), val_rows)
        val_set <- preds[val_rows,]
        val_labels <- labels[val_rows]
        train_set <- preds[train_rows,]
        train_labels <- labels[train_rows]
        if (classifier){
            val_labels <- as.factor(val_labels)
            train_labels <- as.factor(train_labels)
        }else{
            val_labels <- as.integer(val_labels)
            train_labels <- as.integer(train_labels)
        }
        print(head(cbind(train_set, train_labels)))
        bart <- bartMachineCV(train_set,
                              train_labels,
                              use_missing_data=TRUE,
                              replace_missing_data_with_x_j_bar=FALSE,
                              impute_missingness_with_x_j_bar_for_lm=FALSE,
                              verbose=FALSE,
                              serialize=FALSE)
        val_preds <- bart_predict_for_test_data(bart, val_set, val_labels)
        if (classifier){
            neg_val_preds <- which(val_preds$y_hat == FALSE)
            val_preds$p_hat[neg_val_preds] <- 1.0 - val_preds$p_hat[neg_val_preds]
            val_preds <- val_preds$p_hat
        }else{
            val_preds <- val_preds$y_hat
        }
        if (is.null(all_probs)){
            all_probs <- val_preds
            all_labels <- val_labels
        }else{
            all_probs <- cbind(all_probs, val_preds)
            all_labels <- cbind(all_labels, val_labels)
        }
    }
}

rocr_pred <- prediction(all_probs, all_labels)
rocr_roc <- performance(rocr_pred, measure="tpr", x.measure="fpr")
rocr_auc <- performance(rocr_pred, "auc")
mean_auc <- mean(unlist(rocr_auc@y.values), na.rm=TRUE)
se_auc <- sd(unlist(rocr_auc@y.values), na.rm=TRUE)/sqrt(num_reps*k)

file_base <- strsplit(basename(merged_pred_file), split="\\.")[[1]][1]
img_file <- paste0(file_base, "-bart")
if (classifier){
    img_file <- paste(img_file, "classifier",  sep="-")
}
if (use_rand_negs){
    img_file <- paste(img_file, "rand-negs", sep="-")
}
if (directed){
    img_file <- paste(img_file, "direct", sep="-")
}
img_file <- paste0("img/", img_file, "-roc.pdf")

pdf(img_file)
plot(rocr_roc, spread.estimate="boxplot", avg="vertical",
     main=paste0("BART\n(mean AUC=", format(mean_auc, digits=2),
                 ", S.E.M=", format(se_auc, digits=2), ")"))
plot(rocr_roc, 
     main=paste0("BART\n(mean AUC=", format(mean_auc, digits=2),
                 ", S.E.M=", format(se_auc, digits=2), ")"))
dev.off()
