suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx256g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> MERGED_FEATURES DIRECT_VAL_SET N")
}

feat_file <- argv[1]
val_set_file <- argv[2]
n <- as.integer(sub("^\\\\", "", argv[3]))

features <- read.delim(feat_file, as.is=TRUE)
true_intxns_tbl <- read.table(val_set_file, as.is=TRUE)
names(true_intxns_tbl) <- c("prot1", "prot2")

names(features)[1:2] <- c("prot1", "prot2")
good_rows <- which(apply(features[,3:ncol(features)], 1,
                         function(row){
                             any(!is.na(row))
                         }))
features <- features[good_rows,]
features <- subset(features, prot1 != prot2)
rownames(features) <- paste(features$prot1, features$prot2, sep="-")

citations <- read.delim("out/kinase-citations.tsv", as.is=TRUE)
citations <- subset(citations, count.type=="num.subs")

kins <- sort(unique(c(true_intxns_tbl$prot1, true_intxns_tbl$prot2)))

kin <- kins[n]

all_true_intxns <- paste(true_intxns_tbl$prot1, true_intxns_tbl$prot2,
                         sep = "-")
true_intxns_sub <- subset(true_intxns_tbl, prot1 != kin & prot2 != kin)
possible_true_intxns <- paste(true_intxns_sub$prot1, true_intxns_sub$prot2,
                              sep = "-")
possible_true_intxns <- intersect(possible_true_intxns, rownames(features))
possible_false_intxns <- setdiff(rownames(features), possible_true_intxns)
## sample_size <- round(0.8 * min(length(possible_false_intxns),
##                                length(possible_true_intxns)))
num_negs <- min(length(possible_false_intxns),
                8*length(possible_true_intxns))

kin_feats <- subset(features, prot1==kin | prot2==kin)
kin_labels <- rownames(kin_feats) %in% all_true_intxns

num.reps <- 3

all_probs <- NULL
all_labels <- NULL
reps <- NULL
for (i in 1:num.reps){
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
    feats <- features[intxns, 3:ncol(features)]
    labels <- c(rep(TRUE, length(true_intxns)),
                rep(FALSE, length(false_intxns)))
    rand.rows <- sample(1:nrow(feats))
    feats <- feats[rand.rows,]
    labels <- as.factor(labels[rand.rows])
    bart <- bartMachineCV(feats,
                          labels,
                          use_missing_data=TRUE,
                          use_missing_data_dummies_as_covars=TRUE,
                          replace_missing_data_with_x_j_bar=FALSE,
                          impute_missingness_with_x_j_bar_for_lm=FALSE,
                          prob_rule_class=0.5,
                          verbose=FALSE,
                          serialize=FALSE)
    kin_preds <- 1.0-bart_machine_get_posterior(bart, kin_feats[,3:ncol(kin_feats)])$y_hat
    if (is.null(all_probs)){
        all_probs <- kin_preds
        reps <- rep(i, length(kin_preds))
    }else{
        all_probs <- c(all_probs, kin_preds)
        reps <- c(reps, rep(i, length(kin_preds)))
    }
}

kin_preds_df <- data.frame(prot1=kin_feats$prot1,
                           prot2=kin_feats$prot2,
                           rep=reps,
                           lo.kin=rep(kin, length(all_probs)),
                           bart.pred=all_probs)
kin_preds_df$id <- paste(kin_preds_df$prot1, kin_preds_df$prot2, sep="-")

kin_preds_mean <- kin_preds_df %>%
    group_by(id, prot1, prot2, lo.kin) %>%
    summarise(bart.pred.mean=mean(bart.pred, na.rm=TRUE))

kin_preds_mean <- kin_preds_mean[,c("prot1", "prot2", "lo.kin", "bart.pred.mean")]

file_base <- strsplit(feat_file, split="\\.")[[1]][1]

if (!dir.exists(paste0(file_base, "-lo-kin"))){
    dir.create(paste0(file_base, "-lo-kin"))
}

kin_preds_out_df <- subset(kin_preds_mean, prot1==kin)
kin_preds_in_df <- subset(kin_preds_mean, prot2==kin)

write.table(kin_preds_out_df,
            paste0(file_base, "-lo-kin/", sub("out/", "", file_base),
                   "-no-", kin, "-out.tsv"),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

write.table(kin_preds_in_df,
            paste0(file_base, "-lo-kin/", sub("out/", "", file_base),
                   "-no-", kin, "-in.tsv"),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
