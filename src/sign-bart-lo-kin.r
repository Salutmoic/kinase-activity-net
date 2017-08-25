suppressPackageStartupMessages(library(data.table))
options(java.parameters = "-Xmx256g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

argv <- commandArgs(TRUE)
if (length(argv) != 3){
    stop("USAGE: <script> MERGED_SIGN_FEATURES VAL_SET KINASE")
}

feat_file <- argv[1]
val_set_file <- argv[2]
kin <- argv[3]

features <- read.delim(feat_file, as.is=TRUE)
val_set_tbl <- read.table(val_set_file, as.is=TRUE)
names(val_set_tbl) <- c("prot1", "prot2", "sign")
val_set_tbl$sign <- as.factor(val_set_tbl$sign)

names(features)[1:2] <- c("prot1", "prot2")
good_rows <- which(apply(features[,3:ncol(features)], 1,
                         function(row){
                             any(!is.na(row))
                         }))
features <- features[good_rows,]
features <- subset(features, prot1 != prot2)

rownames(features) <- paste(features$prot1, features$prot2, sep="-")
rownames(val_set_tbl) <- paste(val_set_tbl$prot1, val_set_tbl$prot2, sep="-")

final_feats <- c("signed.dcg", "signed.func.score",
                 "kinase.type", "sub.kinase.type", "cor.est.293", "cor.est.272")

good_pairs <- intersect(rownames(features), rownames(val_set_tbl))
feats_sub <- features[good_pairs, ]
feats_sub <- subset(feats_sub, prot1 != kin & prot2 != kin)
val_set_sub <- val_set_tbl[good_pairs, ]
val_set_sub <- subset(val_set_sub, prot1 != kin & prot2 != kin)

kin_feats <- subset(features, prot1==kin | prot2==kin)

num.reps <- 3
samp_size = min(nrow(subset(val_set_sub, sign=="inhibits")),
                nrow(subset(val_set_sub, sign=="activates")))

all_probs <- NULL
reps <- NULL
for (i in 1:num.reps){
    random_pairs <- c(sample(rownames(subset(val_set_sub, sign=="inhibits")),
                             samp_size),
                      sample(rownames(subset(val_set_sub, sign=="activates")),
                             samp_size))
    random_pairs <- sample(random_pairs)
    preds <- feats_sub[random_pairs, final_feats]
    labels <- val_set_sub[random_pairs, "sign"]
    bart <- bartMachineCV(preds,
                          labels,
                          use_missing_data=TRUE,
                          use_missing_data_dummies_as_covars=TRUE,
                          replace_missing_data_with_x_j_bar=FALSE,
                          impute_missingness_with_x_j_bar_for_lm=FALSE,
                          prob_rule_class=0.5,
                          verbose=FALSE,
                          serialize=FALSE)
    kin_preds <- bart_machine_get_posterior(bart, kin_feats[,3:ncol(kin_feats)])$y_hat
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
                           bart.act=all_probs)
kin_preds_df$id <- paste(kin_preds_df$prot1, kin_preds_df$prot2, sep="-")

kin_preds_mean <- kin_preds_df %>%
    group_by(id, prot1, prot2, lo.kin) %>%
    summarise(bart.act.mean=mean(bart.act, na.rm=TRUE))

kin_preds_mean <- kin_preds_mean[,c("prot1", "prot2", "lo.kin", "bart.act.mean")]

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
