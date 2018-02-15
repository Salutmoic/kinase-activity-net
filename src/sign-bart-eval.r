options(java.parameters = "-Xmx16g")
suppressPackageStartupMessages(library(bartMachine))
set_bart_machine_num_cores(24)
suppressPackageStartupMessages(library(ROCR))
library(reshape2)

source("src/rocr-helpers.r")

argv <- commandArgs(TRUE)
if (length(argv) != 4){
    stop("USAGE: <script> MERGED_SIGN_PREDS EDGE_PREDS SIGN_GUESSES VAL_SET")
}

merged_pred_file <- argv[1]
edge_pred_file <- argv[2]
sign_guess_file <- argv[3]
val_set_file <- argv[4]

merged_pred <- read.delim(merged_pred_file, as.is=TRUE)
names(merged_pred)[1:2] <- c("prot1", "prot2")

edge_pred <- read.delim(edge_pred_file, as.is=TRUE)

sign_guess <- read.delim(sign_guess_file, as.is=TRUE)
rownames(sign_guess) <- sign_guess$kinase

rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2, sep="-")

## Temporary fix:
if (any(is.infinite(merged_pred$reg.site.cor.est))){
    inf.reg.site.cor <- which(is.infinite(merged_pred$reg.site.cor.est))
    max.noninf.reg.site.cor <- max(subset(merged_pred, !is.infinite(reg.site.cor.est))$reg.site.cor.est, na.rm=TRUE)
    merged_pred$reg.site.cor.est[inf.reg.site.cor] <- max.noninf.reg.site.cor*sign(merged_pred$reg.site.cor.est[inf.reg.site.cor])
}
if (any(is.infinite(merged_pred$kinact.cor.est))){
    inf.kinact.cor <- which(is.infinite(merged_pred$kinact.cor.est))
    max.noninf.kinact.cor <- max(subset(merged_pred, !is.infinite(kinact.cor.est))$kinact.cor.est, na.rm=TRUE)
    merged_pred$kinact.cor.est[inf.kinact.cor] <- max.noninf.kinact.cor*sign(merged_pred$kinact.cor.est[inf.kinact.cor])
}
if (any(is.nan(merged_pred$kinact.cor.est))){
    nan.kinact.cor <- which(is.nan(merged_pred$kinact.cor.est))
    merged_pred$kinact.cor.est[nan.kinact.cor] <- NA
}

edge_pred_m <- acast(edge_pred, prot1 ~ prot2, value.var="bart.pred")
edge_pred_m <- edge_pred_m[,rownames(edge_pred_m)]
reg_site_cor_m <- acast(merged_pred, prot1 ~ prot2, value.var="reg.site.cor.est")
good_kins <- rownames(edge_pred_m)
reg_site_cor_m <- reg_site_cor_m[good_kins, good_kins]

reg_site_cor_med <- sapply(rownames(reg_site_cor_m),
                           function(kinase){
                               if (!(kinase %in% rownames(edge_pred_m)))
                                   return(NA)
                               kinase_cors <- reg_site_cor_m[kinase,]
                               kinase_edges <- edge_pred_m[kinase,]
                               return(weighted.mean(x=kinase_cors,
                                                    w=kinase_edges,
                                                    na.rm=TRUE))
                           })

reg_site_cor_med[is.nan(reg_site_cor_med)] <- NA

## reg_site_new_pred_m <- t(apply(edge_pred_m, 1,
##                                function(row) return(row*reg_site_cor_med)))
## reg_site_new_pred <- melt(reg_site_new_pred_m)
## names(reg_site_new_pred) <- c("prot1", "prot2", "reg.site.sign")
## merged_pred <- merge(merged_pred, reg_site_new_pred)
merged_pred$two.step.cors <- reg_site_cor_med

## merged_pred$sign.guess <- sapply(rownames(merged_pred),
##                                      function(pair){
##                                          prot2 <- new_merged_pred[pair, "prot2"]
##                                          guess <- sign_guess[prot2, "sign.guess"]
##                                          edge <- edge_pred[pair, "bart.pred"]
##                                          return(guess*edge)
##                                      })
merged_pred$sign.guess <- sign_guess[merged_pred$prot2, "sign.guess"]
rownames(merged_pred) <- paste(merged_pred$prot1, merged_pred$prot2,
                                   sep="-")

val_set_tbl <- read.table(val_set_file, as.is = TRUE)
val_set_tbl[,3] <- as.factor(val_set_tbl[,3])
rownames(val_set_tbl) <- paste(val_set_tbl[,1], val_set_tbl[,2], sep="-")

good_pairs <- intersect(rownames(merged_pred), rownames(val_set_tbl))

merged_pred <- merged_pred[good_pairs,]
val_set_tbl <- val_set_tbl[good_pairs,]

## ## Remove missing data
## new_merged_pred <- new_merged_pred[complete.cases(new_merged_pred),]

random_pairs <- c(sample(rownames(subset(val_set_tbl, V3=="inhibits")), 100),
                  sample(rownames(subset(val_set_tbl, V3=="activates")), 100))
random_pairs <- sample(random_pairs)
preds <- merged_pred[random_pairs, 3:ncol(merged_pred)]
labels <- val_set_tbl[random_pairs, 3]

sink("/dev/null")
bart_init <- bartMachineCV(preds, labels,
                           use_missing_data=TRUE,
                           use_missing_data_dummies_as_covars=TRUE,
                           replace_missing_data_with_x_j_bar=FALSE,
                           impute_missingness_with_x_j_bar_for_lm=FALSE,
                           prob_rule_class=0.5,
                           verbose=FALSE,
                           serialize=FALSE)
sink()

var_select <- var_selection_by_permute(bart_init)

var_select$important_vars_local_names

final_feats <- c("pssm.score", "sign.guess") ## , "reg.site.sign")
## final_feats <- c("pssm.score", "sign.guess", "two.step.cors")

k <- 5
n <- 2
all_probs <- NULL
all_labels <- NULL
for (j in 1:n){
    message(paste("Run:", j))
    random_pairs <- c(sample(rownames(subset(val_set_tbl, V3=="inhibits")), 100),
                      sample(rownames(subset(val_set_tbl, V3=="activates")), 100))
    random_pairs <- sample(random_pairs)
    preds <- merged_pred[random_pairs, final_feats]
    labels <- val_set_tbl[random_pairs, 3]
    val_n <- ceiling(nrow(preds)/k)
    train_n <- nrow(preds) - val_n
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
        sink("/dev/null")
        bart <- bartMachineCV(train_set,
                              train_labels,
                              use_missing_data=TRUE,
                              use_missing_data_dummies_as_covars=TRUE,
                              replace_missing_data_with_x_j_bar=FALSE,
                              impute_missingness_with_x_j_bar_for_lm=FALSE,
                              prob_rule_class=0.0,
                              verbose=FALSE,
                              serialize=FALSE)
        val_preds <- 1.0-bart_machine_get_posterior(bart, val_set)$y_hat
        sink()
        if (is.null(all_probs)){
            all_probs <- val_preds
            all_labels <- val_labels
        }else{
            if (length(val_preds) < val_n){
                num.na <- val_n - length(val_preds)
                val_preds <- c(val_preds, rep(NA, num.na))
                val_labels <- c(val_labels, rep(NA, num.na))
            }
            all_probs <- cbind(all_probs, val_preds)
            all_labels <- cbind(all_labels, val_labels)
        }
    }
}

rocr_pred <- prediction(all_probs, all_labels)
rocr_roc <- performance(rocr_pred, measure="tpr", x.measure="fpr")
rocr_roc2 <- performance(rocr_pred, measure="tnr", x.measure="fnr")
rocr_prec <- performance(rocr_pred, measure="prec", x.measure="rec")
rocr_auc <- performance(rocr_pred, "auc")
rocr_mcc <- performance(rocr_pred, "mat")
mean_auc <- mean(unlist(rocr_auc@y.values), na.rm=TRUE)
se_auc <- sd(unlist(rocr_auc@y.values), na.rm=TRUE)/sqrt(k)

rocr_balance <- performance(rocr_pred, measure="tpr", x.measure="tnr")

avg_mccs <- perf_vert_avg(rocr_mcc)
max_mcc_cutoff <- avg_mccs@x.values[[1]][which.max(avg_mccs@y.values[[1]])]
avg_balance <- perf_thresh_avg(rocr_balance)
approx_max_mcc_cutoff_i <- which.min(abs(avg_balance@alpha.values[[1]]-max_mcc_cutoff))
max_mcc_tnr <- avg_balance@x.values[[1]][approx_max_mcc_cutoff_i]
max_mcc_tpr <- avg_balance@y.values[[1]][approx_max_mcc_cutoff_i]

file_base <- strsplit(basename(merged_pred_file), split="\\.")[[1]][1]
img_file <- paste0(file_base, "-bart-sign")
img_file <- paste0("img/", img_file, "-roc.pdf")

pdf(img_file)
plot(rocr_roc, spread.estimate="stderror", avg="threshold", colorize=TRUE,
     main=paste0("BART\n(mean AUC=", format(mean_auc, digits=2),
                 ", S.E.M=", format(se_auc, digits=2), ")"),
     xlab="False activating rate", ylab="True activating rate")
abline(v=0.05, lty=2, col="blue")
plot(rocr_roc2, spread.estimate="stderror", avg="threshold", colorize=TRUE,
     main=paste0("BART\n(mean AUC=", format(mean_auc, digits=2),
                 ", S.E.M=", format(se_auc, digits=2), ")"),
     xlab="False inhibiting rate", ylab="True inhibiting rate")
abline(v=0.05, lty=2, col="blue")
plot(rocr_mcc, spread.estimate="stderror", avg="vertical")
abline(v=0.5, col="blue", lty=2)
abline(v=max_mcc_cutoff, col="red", lty=2)
plot(rocr_balance, spread.estimate="stderror", avg="threshold", colorize=TRUE,
     main=paste0("BART"), xlab="Average true inhibitory rate",
     ylab="Average true activating rate")
lines(c(max_mcc_tnr, max_mcc_tnr),
      c(-0.03, max_mcc_tpr),
      lty=2, col="blue")
lines(c(-0.03, max_mcc_tnr),
       c(max_mcc_tpr, max_mcc_tpr),
       lty=2, col="blue")
points(max_mcc_tnr, max_mcc_tpr, col="blue", pch=19)
legend("bottomleft",
       legend=c(paste0("max-MCC cutoff: ", format(max_mcc_cutoff, digits=2))),
       lty=2, col=c("blue", "red"), pch=19)
investigate_var_importance(bart)
dev.off()

bart <- bartMachineCV(preds, labels,
                      use_missing_data=TRUE,
                      use_missing_data_dummies_as_covars=TRUE,
                      replace_missing_data_with_x_j_bar=FALSE,
                      impute_missingness_with_x_j_bar_for_lm=FALSE,
                      prob_rule_class=(1-max_mcc_cutoff),
                      verbose=FALSE,
                      serialize=FALSE)

new_preds <- predict(bart, merged_pred[colnames(preds)],
                     type="class", prob_rule_class=0.5)
new_probs <- 1.0-bart_machine_get_posterior(bart, merged_pred[colnames(preds)])$y_hat

out_tbl <- data.frame(prot1=merged_pred$prot1,
                      prot2=merged_pred$prot2,
                      prob.act=new_probs)
out_tbl <- out_tbl[order(out_tbl$prot1, out_tbl$prot2),]

out_file <- paste0(file_base, "-bart-sign-preds")
out_file <- paste0("out/", out_file, ".tsv")

write.table(out_tbl, out_file, quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")
