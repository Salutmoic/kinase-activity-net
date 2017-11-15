suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx24g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> MERGED_PREDS VAL_SET NEG_VAL_SET")
}

merged_pred_file <- argv[1]
val_set_file <- argv[2]
neg_val_set_file <- argv[3]

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")
merged_pred <- subset(merged_pred, prot1 != prot2)
rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2, sep="-")

true_intxns_tbl <- read.table(val_set_file, as.is = TRUE)
possible_true_intxns <- paste(true_intxns_tbl[, 1], true_intxns_tbl[, 2],
                              sep = "-")
possible_true_intxns <- intersect(possible_true_intxns, rownames(merged_pred))

false_intxns_tbl <- read.table(neg_val_set_file, as.is = TRUE)
possible_false_intxns <- paste(false_intxns_tbl[, 1], false_intxns_tbl[, 2],
                               sep = "-")
possible_false_intxns <- intersect(possible_false_intxns, rownames(merged_pred))

## Generate some reverse true interactions to be mixed in as negatives
rev_true_intxns <- paste(true_intxns_tbl[, 2], true_intxns_tbl[, 1], sep = "-")
rev_true_intxns <- intersect(rev_true_intxns, rownames(merged_pred))
rev_true_intxns <- setdiff(rev_true_intxns, possible_true_intxns)
rev_true_intxns <- setdiff(rev_true_intxns, possible_false_intxns)

num_rev_true <- 0.5 * length(rev_true_intxns)
sample_size <- round(0.8 * min(length(possible_false_intxns) + num_rev_true,
                               length(possible_true_intxns)))

num_reps <- 5
k <- 5
all_probs <- NULL
all_labels <- NULL

for (n in 1:num_reps){
    message(paste("Rep.", n))
    false_intxns <- sample(union(possible_false_intxns,
                                 sample(rev_true_intxns,
                                        size=as.integer(num_rev_true))),
                                 size=sample_size)
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
        bart <- bartMachineCV(train_set,
                              as.integer(train_labels),
                              use_missing_data=FALSE,
                              replace_missing_data_with_x_j_bar=FALSE,
                              impute_missingness_with_x_j_bar_for_lm=FALSE,
                              verbose=FALSE,
                              serialize=FALSE)
        val_preds <- bart_predict_for_test_data(bart, val_set, as.integer(val_labels))
        ## neg_val_preds <- which(val_preds$y_hat == FALSE)
        ## val_preds$p_hat[neg_val_preds] <- 1.0 - val_preds$p_hat[neg_val_preds]
        if (is.null(all_probs)){
            ## all_probs <- val_preds$p_hat
            all_probs <- val_preds$y_hat
            all_labels <- val_labels
        }else{
            ## all_probs <- cbind(all_probs, val_preds$p_hat)
            all_probs <- cbind(all_probs, val_preds$y_hat)
            all_labels <- cbind(all_labels, val_labels)
        }
    }
}

rocr_pred <- prediction(all_probs, all_labels)
rocr_roc <- performance(rocr_pred, measure="tpr", x.measure="fpr")
rocr_auc <- performance(rocr_pred, "auc")
mean_auc <- mean(unlist(rocr_auc@y.values), na.rm=TRUE)
se_auc <- sd(unlist(rocr_auc@y.values), na.rm=TRUE)/sqrt(num_reps*k)

pdf("img/bartmachine.pdf")
plot(rocr_roc, spread.estimate="boxplot", avg="vertical",
     main=paste0("BART\n(mean AUC=", format(mean_auc, digits=2),
                 ", S.E.M=", format(se_auc, digits=2), ")"))
plot(rocr_roc, 
     main=paste0("BART\n(mean AUC=", format(mean_auc, digits=2),
                 ", S.E.M=", format(se_auc, digits=2), ")"))
dev.off()
